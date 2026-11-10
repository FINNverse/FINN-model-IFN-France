#' Build the model structure and set its fitted parameters
#'
#' The model definition of this package: the same `createProcess()` and
#' `finn()` calls as the calibration fit script, followed by copying the fitted
#' parameters into the new object. Rebuilding (rather than `torch_load()` of a
#' whole model object) means the model always runs with the methods of the
#' installed FINN.
#'
#' Variant `"process"` (v0.2, `calibration/v02/fit.R`): FINN process model for
#' the 30 species of [species_table()] (no "other" class; valid for stands
#' where species outside the set hold < 2 % of the basal area), annual timestep,
#' covariates `tavg_c` (WorldClim 2.1 bio1, degC), `prec_mm` (bio12, mm) and
#' `gest` (IGN management trace at the start of the interval: 0 none, 1 weak or
#' old, 2 clear signs of exploitation). Each process (growth, mortality =
#' total removal incl. cuts, saturating regeneration) has per-species
#' mechanistic parameters and a two-layer SELU network (20, 20) on the
#' covariates; competition has none. The recruitment dispersion and the
#' per-species Beverton-Holt cap are fixed quantities from the calibration data
#' and travel in `parameters.pt` (`par_theta_recruits_raw`, `buffer.reg_logK`).
#' Calibrated on one patch of 0.0707 ha per site (the IGN plot).
#'
#' The parameter list must carry, besides the tensors of `model$parameters`,
#' the buffers (`buffer.<name>`) and the standardisation constants
#' `env_scaling.<variable>.center` / `.scale`: FINN z-scales the covariates
#' internally and reuses these constants at prediction. Matching is by name and
#' fails loudly on anything missing, unexpected or of another shape.
#'
#' @param variant Name of the model variant (see [model_info()]).
#' @param parameters Named list of tensors, as saved by the calibration script
#'   (`parameters.pt`).
#' @param version Structure of the parameter file: `"v0.1"` (the structure
#'   written out below, the released v0.1 files), `"v0.2-recipe"` (the zoo v0.2
#'   recipe, [build_model_v02()]) or `"auto"` (the default: the recipe when the
#'   file's `settings$recipe` tag says so, v0.1 otherwise).
#' @return A fitted `FINN` model.
#' @keywords internal
build_model = function(variant = "process", parameters, version = c("auto", "v0.1", "v0.2-recipe")) {
  version = match.arg(version)
  ## the structure switch: a zoo v0.2 recipe file describes itself through its
  ## settings$recipe tag (R/recipe_v02.R); everything else is the v0.1 structure below
  if (version == "auto") version = if (is_recipe_v02(parameters)) "v0.2-recipe" else "v0.1"
  if (version == "v0.2-recipe") return(build_model_v02(parameters, variant = if (missing(variant)) NULL else variant))
  if (!identical(variant, "process"))
    stop("Unknown variant '", variant, "'; this package has the variant 'process'.", call. = FALSE)
  n_species = nrow(species_table())                       # 30
  env_vars  = c("tavg_c", "prec_mm", "gest")
  hidden    = c(20L, 20L)
  scal_names = paste0("env_scaling.", rep(env_vars, each = 2), c(".center", ".scale"))
  for (nm in c(scal_names, "buffer.reg_logK", "par_theta_recruits_raw"))
    if (is.null(parameters[[nm]])) stop("parameters.pt lacks ", nm, call. = FALSE)
  ## the fixed cap K comes from the parameter file itself (a buffer, not trainable);
  ## the frozen dispersion travels as the parameter par_theta_recruits_raw and is copied below
  K = as.numeric(torch::torch_exp(parameters[["buffer.reg_logK"]]$to(dtype = torch::torch_float64())))

  ## createProcess() evaluates its formula in the calling frame's parent, so the
  ## formula is written out literally in every call
  model = FINN::finn(
    N_species    = n_species,
    recruits_dbh = 8.3,                                   # median dbh of the observed recruits (tally threshold 7.5 cm)
    reg_floor    = 1e-3,                                  # stems/ha/yr: minimum recruitment
    recruit_obs_weight = pi * 6^2 / 1e4,                  # NB on tallied recruits (the 6 m circle); only used in fitting
    competition_process  = FINN::createProcess(~0, FINN::competition, optimizeSpecies = TRUE),
    growth_process       = FINN::createProcess(~ tavg_c + prec_mm + gest, FINN::growth,
                                               hidden = hidden, dropout = 0.1, optimizeSpecies = TRUE, optimizeEnv = TRUE),
    mortality_process    = FINN::createProcess(~ tavg_c + prec_mm + gest, FINN::mortality,
                                               hidden = hidden, dropout = 0.1, optimizeSpecies = TRUE, optimizeEnv = TRUE),
    regeneration_process = FINN::createProcess(~ tavg_c + prec_mm + gest, FINN::regeneration_saturation,
                                               hidden = hidden, dropout = 0.1, optimizeSpecies = TRUE, optimizeEnv = TRUE,
                                               custom_parameters = list(reg_logK = list(init = log(K), trainable = FALSE)))
  )
  ## the env standardisation learned at fitting (natural units in, z-scores inside)
  model$env_scaling = data.frame(
    variable = env_vars,
    center   = vapply(env_vars, function(v) as.numeric(parameters[[paste0("env_scaling.", v, ".center")]]), numeric(1)),
    scale    = vapply(env_vars, function(v) as.numeric(parameters[[paste0("env_scaling.", v, ".scale")]]), numeric(1)),
    stringsAsFactors = FALSE)
  ## FINN builds the env networks lazily on the first simulation (their input
  ## width comes from the env table): one minimal bare-ground step creates them
  dummy_env = data.frame(siteID = 1L, year = 1L, tavg_c = model$env_scaling$center[1], prec_mm = model$env_scaling$center[2], gest = model$env_scaling$center[3])
  invisible(utils::capture.output(suppressMessages(
    FINN::simulateForest(model, env = dummy_env, init_cohort = NULL, patches = 1L, patch_size = pi * 15^2 / 1e4, device = "cpu"))))

  ## the fitted tensors, matched by NAME; anything missing or unexpected stops here
  own = model$parameters
  missing = setdiff(names(own), names(parameters))
  if (length(missing)) stop("parameters.pt lacks: ", paste(missing, collapse = ", "), call. = FALSE)
  buf_names = paste0("buffer.", names(model$buffers))
  unexpected = setdiff(names(parameters), c(names(own), buf_names, scal_names))
  if (length(unexpected))
    stop("parameters.pt holds tensors this model structure does not have: ", paste(unexpected, collapse = ", "), " (another variant or FINN version?)", call. = FALSE)
  torch::with_no_grad({
    for (nm in names(own)) {
      src = parameters[[nm]]
      if (!identical(as.integer(dim(src)), as.integer(dim(own[[nm]]))))
        stop("parameter '", nm, "' has dims ", paste(dim(src), collapse = "x"), ", the model expects ", paste(dim(own[[nm]]), collapse = "x"), call. = FALSE)
      own[[nm]]$copy_(src$to(device = own[[nm]]$device, dtype = own[[nm]]$dtype))
    }
    for (nm in names(model$buffers)) {
      src = parameters[[paste0("buffer.", nm)]]
      if (is.null(src)) next
      if (!identical(as.integer(dim(src)), as.integer(dim(model$buffers[[nm]]))))
        stop("buffer '", nm, "' has dims ", paste(dim(src), collapse = "x"), ", the model expects ", paste(dim(model$buffers[[nm]]), collapse = "x"), call. = FALSE)
      model$buffers[[nm]]$copy_(src$to(device = model$buffers[[nm]]$device, dtype = model$buffers[[nm]]$dtype))
    }
  })
  model$eval()
  model
}
