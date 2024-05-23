
#########################
#######  ARGUMENTS
#########################

# Index numbers for the multiple imputations. For the final analysis, it should be decks <- 1:10
decks <- 1:10
# Numeric-coded region to subset by in "Region" variable of data.
region_number <- 1
# Variants/strains of interest for competing risk analysis.
variant_names <- c("Ancestral.Lineage",  "Lambda",  "Mu", "Gamma")
#variant_names <- c("Lambda")
#variant_names <- "Gamma" # "Gamma"
# Immune-response markers of interest for analysis.
marker_bindSpike <-  paste0("Day29", c(
  "bindSpike", # ancestral
  "bindSpike_B.1.621", # Mu
  "bindSpike_C.37", # Lambda
  "bindSpike_P.1", #Gamma
  #"bindSpike_B.1.351", # Beta
  "bindSpike_DeltaMDW")) # Delta))
names(marker_bindSpike) <- paste0("Day29", c(
  "bindSpike_Ancestral.Lineage", # ancestral
  "bindSpike_Mu", # Mu
  "bindSpike_Lambda", # Lambda
  "bindSpike_Gamma", #Gamma
  # "bindSpike_Beta", # Beta
  "bindSpike_Delta")) # Delta ))
for(name in names(marker_bindSpike)) {
  for(deck in decks) {
    new_name <- paste0(name, "_", deck)
    old_name <- paste0(marker_bindSpike[[name]], "_", deck)
    data[[new_name]] <- data[[old_name]]
  }
}
markers <- c("Day29bindSpike")
# Markers that differ for each competing variant of interest.
# Example: For variant "variant" and imputation deck "deck", we use
# the marker corresponding to the variable ``marker_name <- paste0(marker, "_", variant, "_", deck)".
markers_varying_by_variant <-  c("Day29pseudoneutid50", "Day29bindSpike")

# Covariates to adjust for in the analysis.
covariates <- "risk_score"
# The name of the event type indicator variable (1 is an observed event and 0 is censoring).
event_type <- "EventIndPrimary"
# The name of the time to event variable.
failure_time <- "EventTimePrimary"






# Hard coded, relevant variables.
ph1 <- "ph1" # phase-1 sampling indicator
Trt <- "Trt" # Vaccine arm indicator
viral_load <- "seq1.log10vl" # viral load
Perprotocol <- "Perprotocol" # perperocol indicator
weights_twostage <- "wt" # weighting for analysis
TwophasesampIndD29 <- "ph2" # two phase sampling indicator.

##################
### Subset to perprotol, phase 1, and region of interest
##################
subset <- which(data[[Perprotocol]] == 1 & data[[ph1]] == 1 & data$Region == region_number)
data <- data[subset]
# Handle the fact that there are two ph2 groups.
# To make downstream analysis more simple, we want to do the analysis using the same observations,
# So we handle this internally by setting weights to zero for
data$ph2.any <- 1*(data[["ph2.D29"]] == 1 | data[["ph2.D29variant"]] == 1)
set(data, , weights_twostage, data[[weights_twostage]] * data[[TwophasesampIndD29]])
TwophasesampIndD29 <- "ph2.any"
##################
### Arguments for threshold analysis
##################
# Max number of thresholdds to consider in grid.
nbins_threshold <- 20

