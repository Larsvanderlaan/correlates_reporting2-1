#-----------------------------------------------
# obligatory to append to the top of each script
renv::activate(project = here::here(".."))

#-----------------------------------------------
# load parameters
print("Starting...")
source(here::here("code", "params.R"))
source(here::here(".", "code/run_survival_analysis.R"))
library(data.table)
print("Files loaded...")
variant_names <- config$variants
print(markers)
print(variant_names)
decks <- 1:10
#decks <- decks[1]
multiply_imputed_markers <- c("Day29pseudoneutid50_Beta",
                              "Day29pseudoneutid50_Gamma")

covariates <- "risk_score"
event_type <- "EventIndPrimary"
failure_time <- "EventTimePrimary"
# TEMPROARY
variant_names <- variant_names[1]
decks <- decks[1:2]
markers <- markers[1]
for(marker in markers){
  for(deck in decks) {
    for(variant in variant_names) {
      variant_type <- paste0("seq1.variant.hotdeck", deck)
      print(marker)
      print(variant_type)
    }


    # Each variant has its own marker.
    run_competing_risk_analysis(covariates = covariates,
                                failure_time = failure_time,
                                event_type = event_type,
                                marker = marker,
                                variant_type = variant_type,
                                variant_names = variant_names)
  }
  # handle multiple imputation replicates
  for(variant in variant_names) {
    event_type_target <- paste0(event_type, "_", variant_type, "_", variant )
    key <- paste0("output/vaccine_", marker, "_", failure_time, "_", event_type_target, ".csv")
    variant_results <- fread(here::here(key))

    results_list <- lapply(decks, function(deck) {
      variant_type <- paste0("seq1.variant.hotdeck", deck)
      event_type_target <- paste0(event_type, "_", variant_type, "_", variant )
      key <- paste0("output/vaccine_", marker, "_", failure_time, "_", event_type_target, ".csv")
      variant_results <- fread(here::here(key))
      variant_results$replicate <- deck
      # variant_results$EIF <- NULL
      return(variant_results)
    })
    results_stacked <- rbindlist(results_list, use.names = TRUE)
    print(head(results_stacked))

    variable_ests <- c("estimates", "estimates_monotone", "estimates_placebo", "estimates_log_RR")
    variable_ses <- c("se", "se", "se_placebo", "se_log_RR")
    names(variable_ses) <- variable_ests




    # add check for thresholds.
    results_combined <- copy(results_stacked[replicate==1])
    new_EIF <- results_stacked[,
                               {
                                 mat_list <- .SD[[1]]
                                 avg_matrix <- Reduce(`+`, mat_list) / length(mat_list)
                                 .(avg_matrix = list(avg_matrix))
                               }
                               , by = threshold, .SDcols = variable_est]
    print(dim(new_EIF))
    print(dim(results_combined))
    set(results_combined, , "EIF", new_EIF)

    stop("hi")


    for(variable_est in variable_ests) {
      variable_se <- variable_ses[[variable_est]]
      set(results_combined,    , variable_est, NA)
      set(results_combined,   , variable_se, NA)
      # get replicate-averaged estimates
      new_estimates <- results_stacked[,  colMeans(.SD) , by = threshold, .SDcols = variable_est]


      indices <- match(results_combined$threshold, new_estimates$threshold)
      set(results_combined, indices  , variable_est, new_estimates[[2]])
      # get new standard errors accounting for replication.
      variable_se <- variable_ses[[variable_est]]
      avg_value <- results_stacked[,  colMeans(.SD) , by = threshold, .SDcols = variable_se][[2]]
      replicate_se <- results_stacked[,  apply(.SD, 2 , sd) , by = threshold, .SDcols = variable_est]
      indices <- match(results_combined$threshold, replicate_se$threshold)
      replicate_se <- replicate_se[[2]]
      replicate_se <- replicate_se[[2]]
      total_se <- sqrt(replicate_se^2 + avg_value^2)
      set(results_combined, indices  , variable_se, total_se)
    }

    key <- paste0("output/vaccine_", marker, "_", failure_time, "_", event_type, "_", variant, ".RData")
    save(results_combined, file = here::here(key))

  }
}



