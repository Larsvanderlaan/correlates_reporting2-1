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
source(here::here("code", "plotting_helpers.R"))
variant_names <- config$variants
for (marker in markers) {
  for(variant in variant_names) {
    get_plot(marker, failure_time = failure_time, event_type = event_type, variant = variant, simultaneous_CI = F, monotone = F, above = TRUE)
    #ggsave(plot, file = here::here(key))
    }
  }
