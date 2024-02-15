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

if(COR == "D29VLancestral") {
  variant_names <- "Ancestral.Lineage"
} else if(COR == "D29VLvariant") {

}
print(COR)
print(variant_names)
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

  }
}

