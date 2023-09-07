

process_data = function(
                              covariates = "standardized_risk_score",
                              failure_time = "EventTimePrimaryD29",
                              event_type = "EventIndPrimaryMolecConfirmedD29",
                              marker,
                              variant_type = "seq1.variant.hotdeck1",
                              variants
){
  # hard coded
  Trt <- "Trt"
  viral_load <- "seq1.log10vl"
  Perprotocol <- "Perprotocol"
  weights_twostage <- wt
  TwophasesampIndD29 <- ph2

  # get the dataset
  #data <- setDT(fread(paste0("~/repositories/covidanalysis/data/janssen_", subset_region, "_partA_data_processed_with_riskscore_hotdeckv4.csv")))
  # subset ph1 and per protocol
  subset <- which(data[[Perprotocol]] == 1 & data[[ph1]] == 1)
  data <- data[subset]



  data <- data[, c(weights_twostage, marker, event_type, failure_time, covariates, variant_type, Perprotocol, TwophasesampIndD29, Trt , viral_load, "EventIndPrimaryD29"), with = FALSE]
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
    tf <- variant_strata[paste0(variant)]
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
    data_treated <- data_treated[, c(covariates, failure_time, event_type_target, marker, weights), with = FALSE]
    data_placebo <-  data_placebo[, c(covariates, failure_time, event_type_target), with = FALSE]
    assert_that(all(
      data_treated$Trt == 1
      & data_treated$Perprotocol == 1
      & data_treated$TwophasesampIndD29 == 1
      & !any(is.na(data_treated[[variant_type]]))
    ))
    assert_that(all(
      data_placebo$Trt == 0
      & data_placebo$Perprotocol == 1
      & data_placebo$TwophasesampIndD29 == 1
      & !any(is.na(data_placebo[[variant_type]]))
    ))





    #### Analysis for treated
    ## get reference time for survival analysis
    output_treated <- run_analysis(tf = tf,
                                   data_surv = data_treated,
                                   covariates = covariates,
                                   failure_time = failure_time,
                                   event_type = event_type_target,
                                   weights = weights,
                                   marker = marker,
                                   nbins_time = 20,
                                   nbins_threshold = 20)

    # make survival dataset
    output_treated <- output_treated[-ncol(output_treated)]
    output_treated$estimates <- unlist( output_treated$estimates )
    output_treated$times <- unlist( output_treated$times )
    output_treated$estimates_monotone <- -as.stepfun(isoreg(output_treated$threshold, -output_treated$estimates))(output_treated$threshold)

    fwrite(output_treated, here(".", "output/vaccine_", marker, "_", failure_time, "_", event_type_target, ".csv"))
  }
}

