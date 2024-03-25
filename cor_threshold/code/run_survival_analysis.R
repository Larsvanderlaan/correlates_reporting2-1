#-----------------------------------------------
# obligatory to append to the top of each script
renv::activate(project = here::here(".."))

#-----------------------------------------------
# load parameters

source(here::here("code", "params.R"))






run_competing_risk_analysis = function(data,
                                       covariates,
                                       failure_time,
                                       event_type,
                                       marker,
                                       variant_type,
                                       variant_names,
                                       weights,
                                       Trt,
                                       TwophasesampIndD29,
                                       viral_load,
                                       threshold_list
){





  for(variant in variant_names) {
    print(variant)
    ########################
    #### Make variant-specific time-to-event variables for competing risk analysis
    ########################

    # Depends on the name of event indicator variable, the variable that labels the variants, and the actual variant of interest.
    event_type_key <- paste0(event_type, "_", variant_type, "_", variant )
    value <- data[[event_type]]
    value[is.na(value)] <- -1
    value[!is.na(value) & value==1] <- ifelse(data[[variant_type]][!is.na(value) & value==1] == variant, 1, 2)
    # Any remaining NAS are assigned a competing risk
    value[value == -1] <- NA
    value[is.na(value)] <- 2
    data[, (event_type_key) := value]

    ########################
    #### Run competing risk analysis
    ########################

    # Reference time for CR
    tf <- tf_by_variant[paste0(variant)]
    # Variant/event type for CR
    event_type_target <- event_type_key


    # make datasets of placebo and treated analysis
    subset_treated <- data[[Trt]] == 1
    data_placebo <- data[!subset_treated]
    data_treated <- data[subset_treated]
    subset <- which(data_treated[[TwophasesampIndD29]] == 1);  data_treated <- data_treated[subset]
    # subset to relevant variables
    data_treated <- data_treated[, c(covariates, failure_time, event_type_target, marker, weights), with = FALSE]

    data_placebo <-  data_placebo[, c(covariates, failure_time, event_type_target), with = FALSE]
    tmp_marker <- ifelse(data_treated[[weights]] == 0, 0, data_treated[[marker]])
    set(data_treated, , marker, tmp_marker)

    if(any(is.na(data_treated[[marker]]))) {
      print(marker)
      print(table(is.na(data_treated[[marker]])))
      stop("NAs in marker for threshold CR")
    }


    print(quantile(data_treated[[marker]]))
    print(threshold_list)

    output_treated <- run_survtmle3(tf = tf,
                                    data_surv = data_treated,
                                    covariates = covariates,
                                    failure_time = failure_time,
                                    event_type = event_type_target,
                                    weights = weights,
                                    marker = marker,
                                    nbins_time = 20,
                                    threshold_list = threshold_list)

    # make survival dataset
    #output_treated <- output_treated[-ncol(output_treated)]
    output_treated$estimates <- unlist( output_treated$estimates )
    output_treated$se <- unlist( output_treated$se )
    output_treated$times <- tf



    n_in_bin <- sapply(threshold_list, function(thresh) {
        sum( data_treated[[marker]] >= thresh)
      })

    n_events_in_bin <- sapply(threshold_list, function(thresh) {
      sum(data_treated[data_treated[[marker]] >= thresh , event_type_target, with = FALSE] == 1)
    })

    output_treated$n_in_bin <- n_in_bin
    output_treated$n_events_in_bin <- n_events_in_bin

    n_event_cutoff <- 10


    weights_for_iso <- sqrt(n_in_bin)
    weights_for_iso[n_events_in_bin < n_event_cutoff] <- 0.01
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






run_survtmle3 <- function(tf, data_survival, covariates, failure_time, event_type, weights , marker = NULL, nbins_time = 20, threshold_list) {
  print("HERE")
  print(tf)
   tf <- as.numeric(tf)
   print(data_survival)
  data_survival <- as.data.table(data_survival)
  print("HERE1")
  # effectivelly removed observations with no weights
  data_survival <- (data_survival[, c(covariates, failure_time, event_type, marker, weights), with = FALSE ])


  # discretize
  time_grid <- unique(quantile(data_survival[[failure_time]], seq(0,1, length = nbins_time+1), type = 1))
  # add time of interest to grid
  time_grid <- sort(union(time_grid, tf))
  failure_time_discrete <- findInterval(data_survival[[failure_time]], time_grid, all.inside = TRUE)
  tf_discrete <- findInterval(tf, time_grid, all.inside = FALSE)
  data_survival[[failure_time]] <-failure_time_discrete


  print("HERE2")

  lrnr <- Lrnr_cv$new(Stack$new(Lrnr_glm$new(), Lrnr_mean$new()))
  #lrnr <- make_learner(Pipeline, Lrnr_cv$new(lrnr, full_fit = TRUE), Lrnr_cv_selector$new(loss_squared_error))


  stack.failure <- stack.censoring <-  Lrnr_hal9001$new(max_degree = 1, num_knots = 20, smoothness_orders = 0)
  learner.event_type <- Lrnr_mean$new()
  learner.treatment <-  Lrnr_glm$new() #Stack$new(Lrnr_glm$new(),  Lrnr_gam$new(), Lrnr_mean$new())
  #learner.treatment <-  make_learner(Pipeline, Lrnr_cv$new(learner.treatment, full_fit = TRUE), Lrnr_cv_selector$new(loss_squared_error))
  treatment <- "treatment"

 learner.treatment <-  Lrnr_glm$new()
 learner.event_type <- Lrnr_glmnet$new()
 stack.censoring <-   stack.failure <-  Lrnr_hal9001$new(max_degree = 1, num_knots = 15, smoothness_orders = 0)


 print("HERE3")

  treatment <- "treatment"
  out_list <- list()
  for(threshold in threshold_list ) {
    print(paste0("THRESHOLD: ", threshold))
    try({
      data_survival[[treatment]] <- 1*(data_survival[[marker]] >= threshold)
      if(sum(data_survival[[treatment]]) == 0){
        print(dim(data_survival))
        print(quantile(data_survival[[marker]]))
        print(marker)
        print(event_type)
        stop("There were zero observation above the threshold. ")
      }
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
                                    tol = 1e-3
      )
    })

      survout$threshold <- threshold
      out_list[[paste0(threshold)]] <- survout
   # })
  }



  output <- rbindlist(out_list)
  return(output)

}



