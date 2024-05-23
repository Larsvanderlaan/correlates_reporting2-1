#-----------------------------------------------
# obligatory to append to the top of each script
renv::activate(project = here::here(".."))

#-----------------------------------------------
# load parameters



source(here::here("code", "params.R"))
source(here::here("code", "params_CR.R"))
library(data.table)



# Iterate
all_outputs_by_deck <- list()
all_outputs_single_by_deck <- list()
for(marker in markers){
  for(deck in decks) {
    variant_type <- paste0("seq1.variant.hotdeck", deck)
    print(variant_type)
    # for(variant in variant_names) {
    #
    #   variant_type <- paste0("seq1.variant.hotdeck", deck)
    #   if(variant != "Ancestral.Lineage") {
    #     marker_name <- paste0(marker, "_", variant, deck)
    #   } else {
    #     marker_name <- marker
    #   }
    #
    #   print(variant_type)
    #   print(marker_name)
    #
    #
    #   run_competing_risk_analysis(covariates = covariates,
    #                               failure_time = failure_time,
    #                               event_type = event_type,
    #                               marker = marker_name,
    #                               variant_type = variant_type,
    #                               variant_names = variant,
    #                               marker_for_thresholds = marker)
    #
    # }



    estimates_list <- list()
    estimates_monotone_list <- list()
    estimates_placebo_list <- list()
    EIF_list <- list()
    thresholds_all <- list()
    thresholds_n_events <- list()
    se_list <- list()
    simult_CI_info_list <- list()
    for(variant_tgt in variant_names) {
      if(variant_tgt != "Ancestral.Lineage") {
        if(marker %in% markers_varying_by_variant) {
          marker_name <- paste0(marker, "_", variant_tgt, "_", deck)
        } else {
          marker_name <- paste0(marker, "_", deck)
        }

      } else {
        marker_name <- paste0(marker, "_", deck)
      }
      event_type_target <- paste0(event_type, "_", variant_type, "_", variant_tgt )
      key <- paste0("output/vaccine_", marker_name, "_", failure_time, "_", event_type_target, ".RDS")
      print(key)
      variant_results <- readRDS(here::here(key))

      estimates_list[[variant_tgt]] <- variant_results$estimates
      estimates_monotone_list[[variant_tgt]] <- variant_results$estimates_monotone

      thresholds_all[[variant_tgt]] <-  variant_results$threshold
      thresholds_n_events[[variant_tgt]] <-  variant_results$n_events_in_bin
      se_list[[variant_tgt]] <-  variant_results$se
      estimates_placebo_list[[variant_tgt]] <- variant_results$estimates_placebo

      #library(ggplot2)
      #plt <- ggplot(data.table(thresholds = variant_results$threshold, estimates = variant_results$estimates), aes(x = thresholds, y = estimates)) + geom_point() + geom_line() #+ geom_hline(yintercept=0)
      #ggsave(plot = plt, filename = here::here(paste0("output/variant_CompRisk_all_", marker, "_", variant_tgt, "_plot.pdf")))


      EIF_list[[variant_tgt]] <- do.call(cbind, variant_results$EIF)
      eif <- do.call(cbind, variant_results$EIF)
      var_mat <- cov(eif)
      rho_mat <- var_mat / sqrt(tcrossprod(diag(var_mat)))
      # Simulunaeous CI for 1-alpha levels 0.8, 0.95
      levels <- c(0.8, 0.95)
      qs <- sapply(levels, function(level) {
        mvtnorm::qmvnorm(level, tail = "both", corr = rho_mat)$quantile
      })
      simult_CI_info <- data.frame(level = levels, quantile = qs)
      simult_CI_info_list[[variant_tgt]] <- simult_CI_info
    }



    output_single_variant <- lapply(variant_names, function(variant) {
      data.table(
        estimates = estimates_list[[variant]],
        estimates_mono = estimates_monotone_list[[variant]],
        estimates_placebo = estimates_placebo_list[[variant]],
        thresholds = thresholds_all[[variant]],
        thresholds_n_events = thresholds_n_events[[variant]],
        se = se_list[[variant]],
        simult_level = simult_CI_info_list[[variant_tgt]]$level,
        simult_quantile = simult_CI_info_list[[variant_tgt]]$quantile
      )
    })
    names(output_single_variant) <- variant_names
    all_outputs_single_by_deck[[paste0(deck)]] <-output_single_variant


    data_treated <- paste0("data_clean/data_treated_", marker_name, "_", failure_time, "_", event_type_target, ".csv")

    # each block of k rows are thresholds estimates, where block corresponds to variant
    estimates <- pmax(unlist(estimates_list), 1e-8)
    estimates_mono <- pmax(unlist(estimates_monotone_list), 1e-8)
    thresholds <- variant_results$threshold
    estimates_placebo <- unlist(estimates_placebo_list)
    thresholds_n_events <- unlist(thresholds_n_events)
    se_list <- unlist(se_list)




    if(length(estimates) != length(variant_names) * length(thresholds)) {
      print(length(estimates))
      print(length(variant_names))
      print(length(thresholds))
      stop("Estimates result lengths dont match.")
    }


    # differing ph2 indicators
    EIF <- do.call(cbind, EIF_list)
    estimates_mat <- matrix(estimates, nrow = nrow(EIF), ncol = ncol(EIF), byrow = TRUE)
    nboot <- 100
    boot_sample <- do.call(rbind, lapply(1:nboot, function(iter){
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
    combos <- t(do.call(rbind, unlist(lapply(variant_names, function(v_tgt) {
      lapply(variant_names, function(v_ref) {
        if(v_tgt == v_ref) {
          return(NULL)
        }
        c(v_tgt, v_ref)
      })
    }), recursive  = FALSE)))


    apply(combos, 2, function(combo) {
      variant_tgt <- combo[1]
      variant_ref <- combo[2]
      index_ref <- match(variant_tgt, variant_names)
      index_reference <- match(variant_ref, variant_names)

      boot_ref <- boot_sample[, seq((index_ref-1)*num_thresh + 1, index_ref * num_thresh, 1)]
      boot_reference <- boot_sample[, seq((index_reference-1)*num_thresh + 1, index_reference * num_thresh, 1)]

      # get true placebo estimates
      est_placebo_ref <- estimates_placebo[seq((index_ref-1)*num_thresh + 1, index_ref * num_thresh, 1)]
      est_placebo_reference <- estimates_placebo[seq((index_reference-1)*num_thresh + 1, index_reference * num_thresh, 1)]
      # get true vaccine arm estimates
      est_vaccine_ref <- estimates[seq((index_ref-1)*num_thresh + 1, index_ref * num_thresh, 1)]
      est_vaccine_reference <- estimates[seq((index_reference-1)*num_thresh + 1, index_reference * num_thresh, 1)]
      # get relative VEs
      real_ests <- log(est_vaccine_ref/est_placebo_ref) - log(est_vaccine_reference/est_placebo_reference)

      est_mono_vaccine_ref <- estimates_mono[seq((index_ref-1)*num_thresh + 1, index_ref * num_thresh, 1)]
      est_mono_vaccine_reference <- estimates_mono[seq((index_reference-1)*num_thresh + 1, index_reference * num_thresh, 1)]
      # get relative VEs
      real_ests_mono <- log(est_mono_vaccine_ref/est_placebo_ref) - log(est_mono_vaccine_reference/est_placebo_reference)


      thresholds_n_events_ref <- thresholds_n_events[seq((index_ref-1)*num_thresh + 1, index_ref * num_thresh, 1)]
      thresholds_n_events_reference <- thresholds_n_events[seq((index_reference-1)*num_thresh + 1, index_reference * num_thresh, 1)]
      thresholds_n_events_min <- pmin(thresholds_n_events_ref, thresholds_n_events_reference)


      # add how to choose variants

      # turn placebo estimates to matrix
      est_placebo_ref <- matrix(est_placebo_ref, nrow = nboot, ncol = num_thresh, byrow = TRUE)
      est_placebo_reference <- matrix(est_placebo_reference, nrow = nboot, ncol = num_thresh, byrow = TRUE)



      boot_ests <- log(boot_ref/est_placebo_ref) - log(boot_reference/est_placebo_reference)
      var_mat <- cov(boot_ests)
      rho_mat <- var_mat / sqrt(tcrossprod(diag(var_mat)))
      # Simulunaeous CI for 1-alpha levels 0.8, 0.95
      levels <- c(0.8, 0.95)
      qs <- sapply(levels, function(level) {
        mvtnorm::qmvnorm(level, tail = "both", corr = rho_mat)$quantile
      })
      simult_CI_info <- data.frame(level = levels, quantile = qs)
      standard_errors <- as.vector(apply(boot_ests, 2, sd))


      output_list <- list(thresholds = thresholds,
                          thresholds_n_events_min=thresholds_n_events_min,
                          target = variant_tgt,
                          reference = variant_ref,
                          bootstrap = boot_ests,
                          estimates = real_ests,
                          estimates_mono = real_ests_mono,
                          standard_errors = standard_errors,
                          simult_CI_info = simult_CI_info)

      all_outputs[[paste0(variant_tgt, "_", variant_ref)]] <<- output_list



    })




    # list of estimates and bootstraps
    #key <- paste0("output/relVE_", marker, "_", failure_time, "_", variant, ".RDS")

    # saveRDS(all_outputs, file = "")

    # stop("hi")
    all_outputs_by_deck[[paste0(deck)]] <- all_outputs

  }
  saveRDS(all_outputs_by_deck, file = here::here(paste0("output/CompRisk_all_", marker, ".RDS")))
  saveRDS(all_outputs_single_by_deck, file = here::here(paste0("output/single_CompRisk_all_", marker, ".RDS")))


}


for(marker in markers){
  results <- readRDS(file = here::here(paste0("output/single_CompRisk_all_", marker, ".RDS")))
  ndeck <- length(results)
  results_1_tmp <- results[[1]]

  combined_results <- rbindlist(lapply(seq_along(results_1_tmp), function(index) {
    estimates_placebo <-  do.call(cbind, lapply(results, function(results_iter) {
      results_iter[[index]]$estimates_placebo
    }))
    estimates_MI <- do.call(cbind, lapply(results, function(results_iter) {
      results_iter[[index]]$estimates
    }))
    estimates_mono_MI <- do.call(cbind, lapply(results, function(results_iter) {
      results_iter[[index]]$estimates_mono
    }))
    se_MI <- do.call(cbind, lapply(results, function(results_iter) {
      results_iter[[index]]$se
    }))
    thresholds_n_events_MI <- do.call(cbind, lapply(results, function(results_iter) {
      results_iter[[index]]$thresholds_n_events
    }))
    simult_CI_info <- as.data.table(rbindlist(lapply(results, function(results_iter) {
      data.frame(level = results_iter[[index]]$simult_level, quantile = results_iter[[index]]$simult_quantile)
    })))



    simult_CI_info <- simult_CI_info[, .(quantile = mean(quantile)), by = level]

    thresholds  <- results[[1]][[index]]$thresholds
    # convert to VE, current version doesn't adjust for randomness in placebo estimate
    estimates_MI <- 1 - estimates_MI / estimates_placebo
    estimates_mono_MI <- 1 - estimates_mono_MI / estimates_placebo
    se_MI <- se_MI / estimates_placebo

    estimates_comb <- rowMeans(estimates_MI)
    estimates_mono_comb <- rowMeans(estimates_mono_MI)
    thresholds_n_events <- rowMeans(thresholds_n_events_MI)
    se_comb <- sqrt(rowMeans(se_MI^2) +
                      (1 + 1/ndeck) *1 / (ndeck - 1) * rowSums((estimates_MI - estimates_comb)^2) )


    return(data.table(comparison = names(results_1_tmp)[index], thresholds = thresholds, thresholds_n_events = thresholds_n_events, estimates = estimates_comb, estimates_mono = estimates_mono_comb,  se = se_comb, simult_CI_info = list(simult_CI_info)))
  }))
  saveRDS(combined_results, file = here::here(paste0("output/single_CompRisk_all_", marker, "_combined.RDS")))

  # sapply(unique(combined_results$comparison), function(comp) {
  #   #comp <- unique(results_comb$comparison)[2]
  #   results <- combined_results[combined_results$comparison==comp,]
  #   library(ggplot2)
  #   plt <- ggplot(results, aes(x = thresholds, y = estimates)) + geom_point() + geom_line() #+ geom_hline(yintercept=0)
  #   ggsave(plot = plt, filename = here::here(paste0("output/CompRisk_all_", marker, "_", comp, "_plot.pdf")))
  # })

}



for(marker in markers){
  results <- readRDS(file = here::here(paste0("output/CompRisk_all_", marker, ".RDS")))
  ndeck <- length(results)
  results_1_tmp <- results[[1]]

  combined_results <- rbindlist(lapply(seq_along(results_1_tmp), function(index) {
    estimates_MI <- do.call(cbind, lapply(results, function(results_iter) {
      results_iter[[index]]$estimates
    }))
    estimates_mono_MI <- do.call(cbind, lapply(results, function(results_iter) {
      results_iter[[index]]$estimates_mono
    }))
    se_MI <- do.call(cbind, lapply(results, function(results_iter) {
      results_iter[[index]]$standard_errors
    }))
    thresholds_n_events_MI <- do.call(cbind, lapply(results, function(results_iter) {
      results_iter[[index]]$thresholds_n_events
    }))
    simult_CI_info <- as.data.table(rbindlist(lapply(results, function(results_iter) {
      results_iter[[index]]$simult_CI_info
    })))
    simult_CI_info <- simult_CI_info[, .(quantile = mean(quantile)), by = level]


    thresholds_n_events <- rowMeans(thresholds_n_events_MI)
    thresholds  <- results[[1]][[index]]$thresholds
    estimates_comb <- rowMeans(estimates_MI)
    estimates_mono_comb <- rowMeans(estimates_mono_MI)
    se_comb <- sqrt(rowMeans(se_MI^2) +
                      (1 + 1/ndeck) *1 / (ndeck - 1) * rowSums((estimates_MI - estimates_comb)^2) )

    return(data.table(comparison = names(results_1_tmp)[index], thresholds = thresholds, estimates = estimates_comb, estimates_mono = estimates_mono_comb, se = se_comb, thresholds_n_events = thresholds_n_events, simult_CI_info = list(simult_CI_info)))
  }))
  saveRDS(combined_results, file = here::here(paste0("output/CompRisk_all_", marker, "_combined.RDS")))

  # sapply(unique(combined_results$comparison), function(comp) {
  #   #comp <- unique(results_comb$comparison)[2]
  #   results <- combined_results[combined_results$comparison==comp,]
  #   library(ggplot2)
  #   plt <- ggplot(results, aes(x = thresholds, y = estimates)) + geom_point() + geom_line() #+ geom_hline(yintercept=0)
  #   ggsave(plot = plt, filename = here::here(paste0("output/CompRisk_all_", marker, "_", comp, "_plot.pdf")))
  # })

}

