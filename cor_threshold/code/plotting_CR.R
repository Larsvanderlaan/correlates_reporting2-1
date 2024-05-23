#-----------------------------------------------
# obligatory to append to the top of each script
renv::activate(project = here::here(".."))

#-----------------------------------------------
library(cowplot)
library(scales)
library(knitr)
library(dplyr)
library(magrittr)
library(ggplot2)
ident <- function(x) x
event_type <- "EventIndPrimary"
failure_time <- "EventTimePrimary"
source(here::here("code", "params.R"))
source(here::here("code", "params_CR.R"))
source(here::here("code", "plotting_helpers_CR.R"))

for(marker in markers) {
  for(variant in variant_names) {
    get_plot_CR(marker, variant = variant, monotone = FALSE)
    get_plot_CR(marker, variant = variant, monotone = TRUE)
    get_plot_CR(marker, variant = variant, monotone = TRUE, simultaneous_CI = TRUE, level = 0.8)
    get_plot_CR(marker, variant = variant, monotone = TRUE, simultaneous_CI = TRUE, level = 0.95)
  }
}



for(marker in markers) {
  for(variant_tgt in variant_names) {
    for(variant_ref in variant_names) {
      if(variant_tgt!=variant_ref) {
        get_plot_CR_quotient(marker, variant_tgt = variant_tgt, variant_ref = variant_ref)
        get_plot_CR_quotient(marker, variant_tgt = variant_tgt, variant_ref = variant_ref, monotone = TRUE)
        get_plot_CR_quotient(marker, variant_tgt = variant_tgt, variant_ref = variant_ref, monotone = TRUE, simultaneous_CI = TRUE, level = 0.8)
        get_plot_CR_quotient(marker, variant_tgt = variant_tgt, variant_ref = variant_ref, monotone = TRUE, simultaneous_CI = TRUE, level = 0.95)

      }
    }
  }
}
