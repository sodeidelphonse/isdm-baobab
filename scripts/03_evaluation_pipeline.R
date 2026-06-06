# Copyright (c) 2026 SODE Akoeugnigan Idelphonse
# Licensed under the MIT License (see LICENSE file for details)
 

#----------------------------------------------
#--- 1) Model fitting and evaluation wrapper 
#----------------------------------------------

#' Fit and evaluate ISDM using a spatial block cross-validation strategy
#' 
#' @param fold integer. The current spatial fold to evaluate
#' @param model_name character. The name of the actual candidate model
#' @param response_type character. The type of response to fit ("joint.po", "count" or "po")
#' @param dataset A named list of spatial datasets to integrate
#' @param mesh A two-dimensional mesh to be created with `fmesher` for the spatial latent process
#' @param proj character. It specify the spatial data projection system 
#' @param boundary An sf polygon object defining the study region
#' @param xy_excluded An optional sf object to specify locations to exclude during the background sampling
#' @param cov_path The path to the serialised raster covariates (useful for parallel process)
#' @param bias_names character. The vector of names of datasets for which a bias field must be define
#' @param metrics character. The vector of metrics to request during the evaluation process
#' @param roc_composite character. The vector of ROC-based metrics to request for the overall score
#' @param prior.sigma The vector of PC prior for the marginal variance
#' @param prior.range The vector of PC prior parameters for the spatial range 
#' @param Offset character. The name of the offset variable in the abundance dataset
#' @param responseCounts character. The name of the count response variable in the count dataset
#' @param has_offset logical. Indicates if the count dataset includes an offset variable.
#' @param int.strategy character. Indicates the integration strategy to use by the INLA engine
#' Defaulted to empirical Bayes ("eb")
#' @param diagonal numeric. A value added to the diagonal elements to stabilize the precision matrix
#' @param seed integer. The seed generator for the reproducibility of the background sample
#' @param verbose logical. Indicates if the details on each model fitting will be shown in the R session
#' @param ... additional arguments passed on to internal functions
#' 
evaluate_model <- function(fold, model_name, response_type = c("joint.po", "count", "po"), dataset,
                           mesh, proj, boundary, xy_excluded, cov_path, bias_names, metrics = NULL, 
                           roc_composite = NULL,  prior.sigma = c(0.1, 0.01), prior.range = c(10, 0.01),  
                           Offset = 'area', responseCounts = 'counts', has_offset = TRUE, int.strategy ='eb',  
                           diagonal = 0.1, seed = 23, verbose = FALSE, ...) {  
  
  response_type <- match.arg(response_type)
  
  #--- Step 1: Data preparation ---
  covariates_pc <- terra::unwrap(readRDS(cov_path)) 
  covariates <- c(covariates_pc, covariates_pc)
  cov_names <- names(covariates_pc)
  cov_pp <- paste0(cov_names, "_pp")
  names(covariates) <- c(cov_names, cov_pp)
  
  cat("========== Processing Fold:", fold, "for Model:", as.character(model_name), "==========\n")  
  folds_splits <- create_folds(dataset, boundary, seed = seed, ...)
  train_data <- extract_fold(folds_splits, fold = fold)$train
  test_data  <- extract_fold(folds_splits, fold = fold)$test
  
  #--- Step 2: Define the models -----
  points_spatial_type <- set_points_spatial(model_name)
  
  if(model_name == "count") {      # separate count model
    model <- startISDM(data = list(Count = train_data$Count), Mesh = mesh, Projection = proj, 
                       responseCounts = responseCounts, Offset = Offset, spatialCovariates = covariates_pc, 
                       pointsSpatial = "shared", Boundary = boundary)
    model$changeComponents("bio14_wc30s(main = bio14_wc30s, model = spde_bio14)")
  } else if (model_name == "po") { # separate PO model
    model <- startISDM(data = list(Presence = train_data$Presence), Mesh = mesh, Projection = proj, 
                       spatialCovariates = covariates_pc, pointsSpatial = "shared", Boundary = boundary)
    model$changeComponents("bio14_wc30s(main = bio14_wc30s, model = spde_bio14)")
  } else if (grepl("_sc", model_name)) {  # dataset-specific covariates' effects
    model <- startISDM(data = train_data, Mesh = mesh, Projection = proj, responseCounts = responseCounts,
                       spatialCovariates = covariates, pointsSpatial = points_spatial_type, 
                       Offset = Offset, Boundary = boundary)
    form_count <- formula(paste('~ . -', paste(cov_pp, collapse = '-')))
    form_po <- formula(paste('~ . -', paste(cov_names, collapse = '-')))
    model$updateFormula(datasetName = "Count", Formula = form_count) 
    model$updateFormula(datasetName = "Presence", Formula = form_po)
    model$changeComponents("bio14_wc30s_pp(main = bio14_wc30s, model = spde_bio14)")
    model$changeComponents("bio14_wc30s(main = bio14_wc30s, model = spde_bio14)")
    
  } else {  
    model <- startISDM(data = train_data, Mesh = mesh, Projection = proj, responseCounts = responseCounts,
                       spatialCovariates = covariates_pc, pointsSpatial = points_spatial_type, 
                       Offset = Offset, Boundary = boundary)
    model$changeComponents("bio14_wc30s(main = bio14_wc30s, model = spde_bio14)")
  }
  
  #--- Step 3: Specify the PC prior for spatial signals and add eventually the bias field -----
  if (points_spatial_type %in% c("shared", "copy", "po", "count")) {
    model$specifySpatial(sharedSpatial = TRUE, PC = TRUE, prior.sigma = prior.sigma, prior.range = prior.range)
    
  } else if (points_spatial_type == "individual") {
    for (ds_name in names(dataset)) {
      model$specifySpatial(sharedSpatial = FALSE, PC = TRUE, datasetName = ds_name, 
                           prior.sigma = prior.sigma, prior.range = prior.range)
    }
  } 
  
  if (grepl("_bias", model_name)) {
    if (!is.null(bias_names)) {
      model$addBias(datasetNames = bias_names, copyModel = FALSE) 
      for (current_bias_name in bias_names) {
        model$specifySpatial(sharedSpatial = FALSE, PC = TRUE, Bias = current_bias_name,
                             prior.sigma = prior.sigma, prior.range = prior.range)
      }
    } else {
      warning(sprintf("Model '%s' specified, but 'bias_names' is NULL. No bias component added.", model_name), call. = FALSE)
    }
  }
  
  #--- Step 4: Fit the model and predict ----
  seed_predict <- seed + 1
  seed_metric  <- seed + 2
  
  spat_fit <- tryCatch({
    fitISDM(model, options = list(control.inla = list(int.strategy = int.strategy, diagonal = diagonal),
                                  control.compute = list(dic = TRUE), 
                                  verbose = verbose)
            )
  }, error = function(e) {
    return(structure(list(error_stage = "fit", message = e$message),
                     class = "model_error"))
  })
  if (inherits(spat_fit, "model_error")) return(spat_fit)
  
  dic_value <- spat_fit$dic$dic
  
  # Predictions
  predictions <- tryCatch({
    predict(spat_fit, data = fm_pix, predictor = TRUE, fun = "linear", 
                           n.samples = 500, seed = seed_predict)
  }, error = function(e) {
    return(structure(list(error_stage = "predict", message = e$message),
                     class = "model_error"))
  })
  if (inherits(predictions, "model_error")) return(predictions)
  
  # Suitability index
  pred_df <- format_predictions(predictions, boundary)
  pred_prob <- suitability_index(pred_df, post_stat = "mean", output_format = "prob", 
                                response_type = response_type, projection = proj, 
                                has_offset = has_offset)
  
  if(response_type == "count") {
    pred_expected <- suitability_index(pred_df, post_stat = "mean", output_format = "response", 
                                        response_type = "count", projection = proj)
  } else if(response_type == "joint.po") {
    predictions_count <- tryCatch({
      predict(spat_fit, data = fm_pix, fun = "linear", spatial = TRUE, intercepts = TRUE, 
              datasets = "Count", covariates = names(covariates_pc), 
              n.samples = 500, seed = seed_predict)
    }, error = function(e) {
      return(structure(list(error_stage = "predict", message = e$message), class = "model_error"))
    })
    if (inherits(predictions_count, "model_error")) return(predictions_count)
    
    pred_count    <- format_predictions(predictions_count, boundary)
    pred_expected <- suitability_index(pred_count, post_stat = "mean", output_format = "response", 
                                       response_type = response_type, projection = proj)
    
  } else {
    pred_expected = NULL
  }
  
  #--- Part 5: Compute the evaluation metrics -----
  metrics_res <- tryCatch({
    switch(response_type,
           joint.po = compute_metrics(test_data, 
                                      prob_raster = pred_prob, 
                                      xy_excluded = xy_excluded, 
                                      expected_response = pred_expected,
                                      metrics = metrics, 
                                      overall_roc_metrics = roc_composite,
                                      best_threshold_policy = "max.f1", 
                                      exposure = Offset, 
                                      is_pred_rate = has_offset,
                                      seed = seed_metric),
           
           po = compute_metrics(list(Presence = test_data$Presence), 
                                prob_raster = pred_prob, 
                                xy_excluded = xy_excluded, 
                                metrics = metrics, 
                                overall_roc_metrics = roc_composite,
                                best_threshold_policy = "max.f1", 
                                seed = seed_metric),
           
           count = compute_metrics(list(Count = test_data$Count), 
                                   prob_raster = pred_prob, 
                                   expected_response =  pred_expected, 
                                   xy_excluded = xy_excluded, 
                                   metrics = metrics, 
                                   overall_roc_metrics = roc_composite,
                                   best_threshold_policy = "max.f1", 
                                   exposure = Offset, 
                                   is_pred_rate = has_offset,
                                   seed = seed_metric)
           
          )
  }, error = function(e) {
    return(structure(list(error_stage = "metrics", message = e$message),
                     class = "model_error"))
  })
  if (inherits(metrics_res, "model_error")) return(metrics_res)
  
  eval_param <- list(Fold = fold, Model = model_name, DIC = dic_value)
  eval_res <- c(eval_param, metrics_res)
  
  return(eval_res)
}


#------------------------------------------------
#--- 2) Logged wrapper for model evaluation 
#------------------------------------------------

#' Logging the model fitting pipeline
#' 
#' @param idx integer. The ID of the current job
#' @param grid data.frame. The combination of the experiment factors across folds
#' @param results_dir character. The directory of the outputs
#' @param ... additional arguments

evaluate_model_logged <- function(idx, grid, results_dir = "results", ...) {
  
  if(!dir.exists(results_dir)) dir.create(results_dir, recursive = TRUE)
  
  fold <- grid$fold[idx]
  model_name <- grid$model[idx]
  model_names <- unique(grid$model)
  response_type <- grid$resp_type[idx]
  
  # Consistent model ID rule
  file_id <- sprintf("%03d", fold + 10 * match(model_name, model_names))
  file_out <- file.path(results_dir, sprintf("model_result_%s.rds", file_id))
  log_file <- file.path(results_dir, "model_run_logging.log")
  
  if (file.exists(file_out)) {
    log_message(sprintf("SKIPPED Fold %s | Model: %s (ID %s)", 
                        fold, model_name, file_id), log_file)
    return(NULL)
  }
  
  start_time <- Sys.time()
  log_message(sprintf("START Fold %s | Model: %s (ID %s)", 
                      fold, model_name, file_id), log_file)
  
  result <- tryCatch({
    eval_res <- evaluate_model(fold, model_name, response_type, ...)
    
    # Detect structured errors 
    if (inherits(eval_res, "model_error")) {
      log_message(sprintf(
        "ERROR Fold %s | Model: %s (ID %s) | Stage: %s | %s",
        fold, model_name, file_id, eval_res$error_stage, eval_res$message
      ), log_file)
      return(eval_res)  
    }
    
    stopifnot(is.list(eval_res))
    saveRDS(eval_res, file_out)
    runtime <- round(difftime(Sys.time(), start_time, units = "mins"), 1)
    log_message(sprintf("DONE Fold %s | Model: %s (ID %s) | Time: %s mins",
                        fold, model_name, file_id, runtime), log_file)
    eval_res
    
  }, error = function(e) {
    log_message(sprintf("ERROR Fold %s | Model: %s (ID %s) | %s",
                        fold, model_name, file_id, e$message), log_file)
    structure(
      list(error_stage = "evaluate_model_logged",
           fold = fold, model = model_name, message = e$message),
      class = "model_error"
    )
  })
  
  return(result)
}

# Logging setup
log_message <- function(msg, log_file = "model_run_logged.log") {
  cat(sprintf("[%s] %s\n", Sys.time(), msg), file = log_file, append = TRUE)
}

#------------------------------------------------------------
#---- 3) Run the full pipeline or resume missing tasks 
#------------------------------------------------------------

#' Resume tasks from the cross-validation pipeline
#' 
#' @param grid A data.frame of the combination of the experiment factors to test 
#' @param results_dir character. The directory of the output
#' @param timeout The time limit to skip an unsuccessful task from the whole pipeline
#' @param workers Number of workers in case of parallel processing (optional)
#' @param ... Additional arguments
resume_tasks <- function(grid, results_dir = "results", timeout = 3600, workers = 5, ...) {
  
  if (!dir.exists(results_dir)) dir.create(results_dir, recursive = TRUE)
  
  # The same ID rule
  model_names <- unique(grid$model)
  all_ids <- sprintf("%03d", grid$fold + 10 * match(grid$model, model_names))
  
  existing_ids <- gsub("model_result_|\\.rds", "", list.files(results_dir, pattern = "\\.rds$"))
  missing_indices <- which(!all_ids %in% existing_ids)
  grid_missing <- grid[missing_indices, ]
  
  message(sprintf("Resuming %d missing tasks", nrow(grid_missing)))
  if (nrow(grid_missing) == 0) return(invisible(list()))
  
  # Run the tasks sequentially (safe version, but can be paralleled)
  results <- lapply(seq_along(missing_indices), function(i) {
    idx <- missing_indices[i]
    fold <- grid$fold[idx]
    model_name <- grid$model[idx]
    
    tryCatch({
      R.utils::withTimeout({
        evaluate_model_logged(idx, grid, results_dir = results_dir, ...)
      }, timeout = timeout, onTimeout = "error")
    }, error = function(e) {
      log_message(sprintf("ERROR [Retry] Fold %s | Model %s: %s", fold, model_name, e$message),
                  file.path(results_dir, "model_run_logging.log"))
      structure(
        list(error_stage = "resume_tasks",
             fold = fold, model = model_name,
             message = e$message),
        class = "model_error"
      )
    })
  })
  
  return(results)
}

#---------------------------------------------------------
#--- 4) Compile results after all tasks are completed
#---------------------------------------------------------

compile_results <- function(results_dir = "results") {
  
  files <- list.files(results_dir, pattern = "^model_result_.*\\.rds$", full.names = TRUE)
  if (length(files) == 0) {
    warning("No result files found in ", results_dir)
    return(NULL)
  }
  res_list <- lapply(files, readRDS)
  res_df <- dplyr::bind_rows(res_list)
  
  return(res_df)
}


