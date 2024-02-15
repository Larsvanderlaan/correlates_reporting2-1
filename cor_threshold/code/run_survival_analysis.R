#-----------------------------------------------
# obligatory to append to the top of each script
renv::activate(project = here::here(".."))

#-----------------------------------------------
# load parameters

source(here::here("code", "params.R"))






run_competing_risk_analysis = function(covariates,
                                       failure_time,
                                       event_type,
                                       marker,
                                       variant_type,
                                       variant_names,
                                       marker_for_thresholds
){


  # hard coded
  ph1 <- "ph1.D29"
  Trt <- "Trt"
  viral_load <- "seq1.log10vl"
  Perprotocol <- "Perprotocol"
  weights_twostage <- "wt"
  TwophasesampIndD29 <- "TwophasesampIndD29variant"
  # get the dataset
  #data <- setDT(fread(paste0("~/repositories/covidanalysis/data/janssen_", subset_region, "_partA_data_processed_with_riskscore_hotdeckv4.csv")))
  # subset ph1 and per protocol
  subset <- which(data[[Perprotocol]] == 1 & data[[ph1]] == 1)
  data <- data[subset]
  data <- data[data$Region == 1]


  data <- data[, c(weights_twostage, marker, event_type, failure_time, covariates, variant_type, Perprotocol, TwophasesampIndD29, Trt , viral_load, marker_for_thresholds), with = FALSE]
  # make competing risk indicators
  variant_strata <- c("Ancestral.Lineage" = 181 ,
                      "Zeta" = 176,
                      "Lambda" = 77,
                      "Mu" = 175,
                      "Gamma" = 181
  )



  #names(variant_strata)
  for(variant in variant_names) {
    print(variant)
    event_type_key <- paste0(event_type, "_", variant_type, "_", variant )
    value <- data[[event_type]]
    value[is.na(value)] <- -1
    value[!is.na(value) & value==1] <- ifelse(data[[variant_type]][!is.na(value) & value==1] == variant, 1, 2)
    # Any remaining NAS are assigned a competing risk
    value[value == -1] <- NA
    value[is.na(value)] <- 2
    data[, (event_type_key) := value]




    # use competing weights
    weights <- weights_twostage

    # Run competing risk analysis
    tf <- tf_by_variant[paste0(variant)]
    event_type_target <- event_type_key
    # for CR only, remove observations without variant information
    #if(event_type != event_type_target) {
    # data <- data[!is.na(data[[variant_type]])]
    #}

    # make datasets of placebo and treated analysis
    subset_treated <- data[[Trt]]==1
    data_placebo <- data[!subset_treated]
    data_treated <- data[subset_treated]
    subset <- which(data_treated[[TwophasesampIndD29]] == 1)
    data_treated <- data_treated[subset]  # assumes TwophasesampIndD29 used only for treatment arm
    # subset to reelvant variables
    marker_values_for_thresh <- data_treated[[marker_for_thresholds]]

    all_thresholds <- sort(as.vector(na.omit(marker_values_for_thresh[data_treated[[event_type]] != 0])))
    data_treated <- data_treated[, c(covariates, failure_time, event_type_target, marker, weights), with = FALSE]
    data_placebo <-  data_placebo[, c(covariates, failure_time, event_type_target), with = FALSE]



    if(TRUE) {
      print(quantile(all_thresholds))
      nbins_threshold <- 20

      # ensure there are at least 5 events above threshold.
      drop_thresh <- min(unique(all_thresholds)[order(unique(all_thresholds), decreasing = TRUE)[1:5]])
      #drop_thresh <- c(drop_thresh, quantile(data_treated[[marker]], 0.95, type =1))
      all_thresholds <- all_thresholds[all_thresholds <= drop_thresh]


      # add minimal threshold and make grid of thresholds
      threshold_list <- min(marker_values_for_thresh, na.rm = TRUE)
      threshold_list <- sort(unique(
        c(
          threshold_list,
          quantile(setdiff(all_thresholds, threshold_list), seq(0, 1, length = nbins_threshold), type = 1)
        )))

      if(any(is.na(marker_values_for_thresh))) {
        stop("NAs in marker for threshold CR")
      }
      marker_data <- marker_values_for_thresh
      n_in_bin <- sapply(threshold_list, function(s) {
        sum(marker_data >= s)
      })

      threshold_list <- threshold_list[n_in_bin >= 50]
    }




    output_treated <- run_survtmle3(tf = tf,
                                    data_surv = data_treated,
                                    covariates = covariates,
                                    failure_time = failure_time,
                                    event_type = event_type_target,
                                    weights = weights,
                                    marker = marker,
                                    nbins_time = 20,
                                    threshold_list = threshold_list,
                                    nbins_threshold = 20)

    # make survival dataset
    #output_treated <- output_treated[-ncol(output_treated)]
    output_treated$estimates <- unlist( output_treated$estimates )
    output_treated$se <- unlist( output_treated$se )
    output_treated$times <- tf
    weights_for_iso <- sqrt(n_in_bin)
    weights_for_iso <- weights_for_iso / sum(weights_for_iso)

    output_treated$estimates_monotone <- -isotone::gpava(output_treated$threshold, -output_treated$estimates, weights = weights_for_iso)$x

    #output_treated$estimates_monotone <- -as.stepfun(isoreg(output_treated$threshold, -output_treated$estimates))(output_treated$threshold)




    # run placebo analysis
    data_placebo <- as.data.table(na.omit(data_placebo))
    #form <- as.formula(paste0("Surv(", failure_time,  ", as.factor(", event_type_target, ")", ") ~ 1"))
    fit <- cmprsk::cuminc(data_placebo[[failure_time]], data_placebo[[event_type_target]])
    fit <- cmprsk::timepoints(fit, as.numeric(tf))

    # get estimate and se for reference time
    est <- fit$est[1,1]
    se <- sqrt(fit$var[1,1])

    # delta method log(est) ~ sd(IF)/est
    output_treated$estimates_placebo <- est
    output_treated$se_placebo <- se
    output_treated$estimates_log_RR <-  log(output_treated$estimates) - log(est)
    output_treated$se_log_RR <- sqrt((output_treated$se/output_treated$estimates)^2 + (se/est)^2)

    saveRDS(output_treated, file = here::here(paste0("output/vaccine_", marker, "_", failure_time, "_", event_type_target, ".RDS")))
    fwrite(data_treated,  here::here(paste0("data_clean/data_treated_", marker, "_", failure_time, "_", event_type_target, ".csv")))
    fwrite(data_placebo,   here::here(paste0("data_clean/data_placebo_", marker, "_", failure_time, "_", event_type_target, ".csv")))
  }
}






run_survtmle3 <- function(tf, data_survival, covariates, failure_time, event_type, weights , marker = NULL, nbins_time = 20, threshold_list = NULL, nbins_threshold = 20) {
  tf <- as.numeric(tf)
  data_survival <- as.data.table(data_survival)
  # effectivelly removed observations with no weights
  data_survival <- (data_survival[, c(covariates, failure_time, event_type, marker, weights), with = FALSE ])


  # discretize
  time_grid <- unique(quantile(data_survival[[failure_time]], seq(0,1, length = nbins_time+1), type = 1))
  # add time of interest to grid
  time_grid <- sort(union(time_grid, tf))
  failure_time_discrete <- findInterval(data_survival[[failure_time]], time_grid, all.inside = TRUE)
  tf_discrete <- findInterval(tf, time_grid, all.inside = FALSE)
  data_survival[[failure_time]] <-failure_time_discrete
  print(tf)
  print(tf_discrete)




  # TODO

  lrnr <- Lrnr_cv$new(Stack$new(Lrnr_glm$new(), Lrnr_mean$new()))
  #lrnr <- make_learner(Pipeline, Lrnr_cv$new(lrnr, full_fit = TRUE), Lrnr_cv_selector$new(loss_squared_error))

  stack.failure <- stack.censoring <- Lrnr_cv$new(Lrnr_gam$new())
  learner.event_type <- Lrnr_cv$new(Lrnr_glmnet$new())
  learner.treatment <-  Lrnr_cv$new(Lrnr_glm$new()) #Stack$new(Lrnr_glm$new(),  Lrnr_gam$new(), Lrnr_mean$new())
  #learner.treatment <-  make_learner(Pipeline, Lrnr_cv$new(learner.treatment, full_fit = TRUE), Lrnr_cv_selector$new(loss_squared_error))
  treatment <- "treatment"

  stack.failure <- stack.censoring <- learner.event_type <- learner.treatment <-  Lrnr_glm$new()

  # TEMP
  # learner.treatment <- stack.failure <- stack.censoring <- learner.event_type <- Lrnr_cv$new(Lrnr_glm$new())

  # get thresholds corresponding to events
  if(is.null(threshold_list)) {
    all_thresholds <- sort(as.vector(na.omit(data_survival[[marker]][data_survival[[event_type]] != 0])))
    # ensure there are at least 5 events above threshold.
    drop_thresh <- min(unique(all_thresholds)[order(unique(all_thresholds), decreasing = TRUE)[1:10]])
    #drop_thresh <- min(drop_thresh, quantile( all_thresholds, 0.95, type = 1))
    all_thresholds <- all_thresholds[all_thresholds <= drop_thresh]


    # add minimal threshold and make grid of thresholds
    threshold_list <- min(data_survival[[marker]], na.rm = TRUE)
    threshold_list <- sort(unique(
      c(
        threshold_list,
        quantile(setdiff(all_thresholds, threshold_list), seq(0, 1, length = nbins_threshold), type = 1)
      )))

    if(any(is.na(data_survival[[marker]]))) {
      stop("NAs in marker for threshold CR")
    }
    marker_data <- data_survival[[marker]]
    n_in_bin <- sapply(threshold_list, function(s) {
      sum(marker_data >= s)
    })
    print(n_in_bin)
    threshold_list <- threshold_list[n_in_bin >= 60]
  }
  #max_cutoff <- quantile(data[[marker]], 0.975, na.rm = TRUE)



  treatment <- "treatment"
  out_list <- list()
  for(threshold in threshold_list  ) {
    print(paste0("THRESHOLD: ", threshold))
    #try({
      data_survival[[treatment]] <- 1*(data_survival[[marker]] >= threshold)

      #print(mean(1*(data_survival[[marker]] >= threshold)))

      survout <- survtmle3_discrete(data_survival[[failure_time]], data_survival[[event_type]],
                                    data_survival[[treatment]], data_survival[, covariates, with = FALSE],
                                    weights = data_survival[[weights]],
                                    learner.treatment =  learner.treatment,
                                    learner.failure_time =  stack.failure,
                                    learner.censoring_time = stack.censoring,
                                    learner.event_type = learner.event_type,
                                    target_failure_time = tf_discrete,
                                    target_treatment = c(1),
                                    target_event_type = 1,
                                    failure_time.stratify_by_time = FALSE,
                                    censoring_time.stratify_by_time = FALSE,
                                    cross_fit = FALSE,
                                    cross_validate = FALSE,
                                    calibrate = FALSE,
                                    verbose = TRUE, max_iter = 100,
                                    tol = 1e-5
      )

      survout$threshold <- threshold
      out_list[[paste0(threshold)]] <- survout
   # })
  }



  output <- rbindlist(out_list)
  return(output)

}



