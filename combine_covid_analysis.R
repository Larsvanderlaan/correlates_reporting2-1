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
                              "bindSpike", "bindSpike_B.1.621", "bindSpike_C.37", "bindSpike_P.1",
                              "bindSpike_B.1.351", "bindSpike_DeltaMDW",
                              "pseudoneutid50", "pseudoneutid50_delta",
                              "pseudoneutid50_beta", "pseudoneutid50_zeta",
                              "pseudoneutid50_mu", "pseudoneutid50_gamma",
                              "pseudoneutid50_lambda")

covariates <- "risk_score"
event_type <- "EventIndPrimary"
failure_time <- "EventTimePrimary"
# TEMPROARY
#variant_names <- c("Ancestral.Lineage", "Gamma") #variant_names
decks <- decks[1:10]
markers <- "Day29pseudoneutid50"
all_outputs_by_deck <- list()
for(marker in markers){
  for(deck in decks) {
    for(variant in variant_names) {

      variant_type <- paste0("seq1.variant.hotdeck", deck)
      if(variant != "Ancestral.Lineage") {
        marker_name <- paste0(marker, "_", variant, deck)
      } else {
        marker_name <- marker
      }

      print(variant_type)
      print(marker_name)


      run_competing_risk_analysis(covariates = covariates,
                                  failure_time = failure_time,
                                  event_type = event_type,
                                  marker = marker_name,
                                  variant_type = variant_type,
                                  variant_names = variant,
                                  marker_for_thresholds = marker)

    }



    estimates_list <- list()
    estimates_placebo_list <- list()
    EIF_list <- list()
    thresholds_all <- c()
    for(variant_ref in variant_names) {
      if(variant_ref != "Ancestral.Lineage") {
        marker_name <- paste0(marker, "_", variant_ref, deck)
      } else {
        marker_name <- "Day29pseudoneutid50"
      }
      event_type_target <- paste0(event_type, "_", variant_type, "_", variant_ref )
      key <- paste0("output/vaccine_", marker_name, "_", failure_time, "_", event_type_target, ".RDS")
      print(key)
      variant_results <- readRDS(here::here(key))
      print(variant_results$estimates)
      estimates_list[[variant_ref]] <- variant_results$estimates
      thresholds_all <- c(thresholds_all, variant_results$threshold)
      estimates_placebo_list[[variant_ref]] <- variant_results$estimates_placebo
      EIF_list[[variant_ref]] <- do.call(cbind, variant_results$EIF)
    }





    # each block of k rows are thresholds estimates, where block corresponds to variant
    estimates <- pmax(unlist(estimates_list), 1e-8)
    thresholds <- variant_results$threshold
    estimates_placebo <- unlist(estimates_placebo_list)




    if(length(estimates) != length(variant_names) * length(thresholds)) {
      stop("Estimates result lengths dont match.")
    }



    EIF <- do.call(cbind, EIF_list)
    estimates_mat <- matrix(estimates, nrow = nrow(EIF), ncol = ncol(EIF), byrow = TRUE)

    boot_sample <- do.call(rbind, lapply(1:10000, function(iter){
      index <- sample(1:nrow(EIF), nrow(EIF), replace = TRUE)
      EIF_tmp <- estimates_mat + EIF
      colMeans(EIF[index , , drop = FALSE])
    }))
    # remove bias
    boot_sample <-  boot_sample+  matrix(estimates - colMeans(boot_sample) , nrow = nrow(boot_sample), ncol = ncol(boot_sample), byrow = TRUE)

    boot_sample <- apply(boot_sample, 2, pmax, 1e-8)

    #print(as.vector(apply(boot_sample, 2, sd))[1:6])
    #print(as.vector(sqrt(diag(var(EIF))) / sqrt(nrow(EIF)))[1:6])


    num_thresh <- length(thresholds)
    all_outputs <- list()
    combos <- combn(variant_names, 2)
    apply(combos, 2, function(combo) {
      variant_ref <- combo[1]
      variant_control <- combo[2]
      index_ref <- match(variant_ref, variant_names)
      index_control <- match(variant_control, variant_names)
      print(index_ref)
      print(index_control)
      boot_ref <- boot_sample[, seq((index_ref-1)*num_thresh + 1, index_ref * num_thresh, 1)]
      boot_control <- boot_sample[, seq((index_control-1)*num_thresh + 1, index_control * num_thresh, 1)]

      # get true placebo estimates
      est_placebo_ref <- estimates_placebo[seq((index_ref-1)*num_thresh + 1, index_ref * num_thresh, 1)]
      est_placebo_control <- estimates_placebo[seq((index_control-1)*num_thresh + 1, index_control * num_thresh, 1)]
      # get true vaccine arm estimates
      est_vaccine_ref <- estimates[seq((index_ref-1)*num_thresh + 1, index_ref * num_thresh, 1)]
      est_vaccine_control <- estimates[seq((index_control-1)*num_thresh + 1, index_control * num_thresh, 1)]
      # get relative VEs
      real_ests <- log(est_vaccine_ref/est_placebo_ref) - log(est_vaccine_control/est_placebo_control)

      # add how to choose variants

      # turn placebo estimates to matrix
      est_placebo_ref <- matrix(est_placebo_ref, nrow = 10000, ncol = num_thresh, byrow = TRUE)
      est_placebo_control <- matrix(est_placebo_control, nrow = 10000, ncol = num_thresh, byrow = TRUE)



      boot_ests <- log(boot_ref/est_placebo_ref) - log(boot_control/est_placebo_control)
      standard_errors <- as.vector(apply(boot_ests, 2, sd))


      output_list <- list(thresholds = thresholds,
                          reference = variant_ref,
                          control = variant_control,
                          bootstrap = boot_ests,
                          estimates = real_ests,
                          standard_errors = standard_errors)

      all_outputs[[paste0(variant_ref, "_", variant_control)]] <<- output_list



    })




    # list of estimates and bootstraps
    #key <- paste0("output/relVE_", marker, "_", failure_time, "_", variant, ".RDS")

    # saveRDS(all_outputs, file = "")

    # stop("hi")
    all_outputs_by_deck[[paste0(deck)]] <- all_outputs

  }
  saveRDS(all_outputs_by_deck, file = here::here(paste0("output/CompRisk_all_", marker, ".RDS")))


}


for(marker in markers){
  results <- readRDS(file = here::here(paste0("output/CompRisk_all_", marker, ".RDS")))
  ndeck <- length(results)
  results_1 <- results[[1]]

  combined_results <- rbindlist(lapply(seq_along(results_1), function(index) {
    estimates_MI <- do.call(cbind, lapply(results, function(results_iter) {
      results_iter[[index]]$estimates
    }))
    se_MI <- do.call(cbind, lapply(results, function(results_iter) {
      results_iter[[index]]$standard_errors
    }))
    thresholds  <- results[[1]][[index]]$thresholds
    estimates_comb <- rowMeans(estimates_MI)
    se_comb <- sqrt(rowMeans(se_MI^2) +
                      (1 + 1/ndeck) *1 / (ndeck - 1) * rowSums((estimates_MI - estimates_comb)^2) )

    return(data.table(comparison = names(results_1)[index], thresholds = thresholds, estimates = estimates_comb, se = se_comb))
  }))
  saveRDS(combined_results, file = here::here(paste0("output/CompRisk_all_", marker, "_combined.RDS")))

  sapply(unique(combined_results$comparison), function(comp) {
    #comp <- unique(results_comb$comparison)[2]
    results <- combined_results[combined_results$comparison==comp,]
    library(ggplot2)
    plt <- ggplot(results, aes(x = thresholds, y = estimates)) + geom_point() + geom_line() #+ geom_hline(yintercept=0)
    ggsave(plot = plt, filename = here::here(paste0("output/CompRisk_all_", marker, "_", comp, "_plot.RDS")))
  })

}

