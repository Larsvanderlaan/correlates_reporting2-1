#-----------------------------------------------
# obligatory to append to the top of each script
renv::activate(project = here::here(".."))

#-----------------------------------------------
# load parameters
print("Starting...")
source(here::here("code", "params.R"))
source(here::here(".", "code/run_survival_analysis.R"))
source(here::here("code", "params_CR.R"))
print("Files loaded...")


# Iterate over all markers, decks, and variants.
# Subset variants relevant to the specified COR argument.
# Note: different COR arguments correspond to different variants (due to differs ph2 sampling procedures).
if(COR == "D29VLancestral") {
  variant_names <- "Ancestral.Lineage"
} else if(COR == "D29VLvariant") {
  variant_names <- setdiff(variant_names, "Ancestral.Lineage")
}
print(COR)
print(variant_names)
for(marker in markers){
  ({
    print(marker)



    ########################
    ### Determine marker-specific grid of thresholds for the threshold risk analysis.
    ### Note:
    ### - The threshold grid differs by markers.
    ### - For a given marker, the threshold grid is identical for each variant and imputation hot deck iteration.
    ########################



    # # make competing risk indicators
    # variant_strata <- c("Ancestral.Lineage" = 181 ,
    #                     "Zeta" = 176,
    #                     "Lambda" = 77,
    #                     "Mu" = 175,
    #                     "Gamma" = 181
    # )


    # Get threshold grid based on marker values for ancestral strain (hot deck 1).

    subset <- which((data[[Trt]] == 1 & data[["ph2.D29"]] == 1))
    data_anc_ref <- data[subset]
    marker_values_for_thresh <- na.omit(data_anc_ref[[paste0(marker, "")]])
    min_thresholds <- sapply(union(variant_names, "Ancestral.Lineage"), function(vname) {
      min(na.omit(data_anc_ref[[paste0(marker, "_", vname, "_1")]]))
    })


    # Only keep marker values for observations where an event occurs.
    all_thresholds <- sort(as.vector(na.omit(marker_values_for_thresh[data_anc_ref[[event_type]] != 0])))
    all_thresholds <- union(min_thresholds, all_thresholds)
    # CHECK FOR NAs in marker values





    #########
    ###### Make threshold grid
    ###### Subset to thresholds with at least 5 total events and at least 50 observations above threshold
    #######

    # Subset to thresholds with at least 20 total events
    drop_thresh <- min(unique(all_thresholds)[order(unique(all_thresholds), decreasing = TRUE)[1:10]])
    all_thresholds <- all_thresholds[all_thresholds <= drop_thresh]

    # add minimal threshold and make grid of thresholds
    thresholds_to_include <- min(marker_values_for_thresh, na.rm = TRUE)
    threshold_list <- sort(unique(
      c(
        thresholds_to_include,
        quantile(setdiff(all_thresholds, thresholds_to_include), seq(0, 1, length = nbins_threshold), type = 1)
      )))
    marker_data <- marker_values_for_thresh

    # subset to thresholds with at least 50 observations above threshold
    n_in_bin <- sapply(threshold_list, function(s) {
      sum(marker_data >= s)
    })
    threshold_list <- threshold_list[n_in_bin >= 60]




    # iterate over imputation decks and variant_names
    for(deck in decks) {
      for(variant in variant_names) {
        # Get
        variant_type_deck <- paste0("seq1.variant.hotdeck", deck)
        if(variant != "Ancestral.Lineage") {
          if(marker %in% markers_varying_by_variant) {
            marker_name <- paste0(marker, "_", variant, "_", deck)
          } else {
            marker_name <- paste0(marker, "_", deck)
            stop("TEMP")
          }

        } else {
          marker_name <-  paste0(marker, "_", deck) #marker
        }


        print(variant_type_deck)
        print(marker_name)


        run_competing_risk_analysis(data = data, covariates = covariates,
                                    failure_time = failure_time,
                                    event_type = event_type,
                                    marker = marker_name,
                                    variant_type = variant_type_deck,
                                    variant_names = variant,
                                    weights = weights_twostage,
                                    Trt = Trt,
                                    TwophasesampIndD29 = TwophasesampIndD29,
                                    viral_load = viral_load,
                                   threshold_list = threshold_list)

      }

    }
  })
}

