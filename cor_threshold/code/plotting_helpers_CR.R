
source(here::here("code", "params_CR.R"))

get_plot_CR <- function(marker, variant = "", simultaneous_CI = F, monotone = F, level = 0.95) {

  ###### Some variables.
  variant_type <- paste0("seq1.variant.hotdeck", 1) # Used to extract marker data, any hot deck suffices.
  failure_time = "EventTimePrimary" # failure time variable name
  event_type = "EventIndPrimary" # failure type variable name

  library(survival)
  library(data.table)
  # Variant-specific RESULTS for threshold competing risk analysis.
  key <- paste0("output/single_CompRisk_all_", marker, "_combined.RDS")
  print(key)
  results <- as.data.table(readRDS(here::here(key)))


  # Get data containing marker information for RCDF marker plots.
  event_type_target <- paste0(event_type, "_", variant_type, "_", variant )
  data_treated <- fread(here::here(paste0("data_clean/data_treated_", marker, "_", failure_time, "_", event_type_target, ".csv")))
  print(names(data_treated))
  set(data_treated, , event_type ,data_treated[, grep("EventInd", names(data_treated), value = TRUE), with = FALSE])
  set(data_treated, , failure_time ,data_treated[, grep("EventTime", names(data_treated), value = TRUE), with = FALSE])
  data_treated <- na.omit(data_treated)
  data_treated <- as.data.frame(data_treated)
  event_type_target <- as.data.frame(event_type_target)
  print(names(data_treated))


  key <- marker
  laby <- paste0("Vaccine Efficacy")
  labx <- plotting_assay_label_generator(marker, above = TRUE)
  main <- plotting_assay_title_generator(marker)
  ident <- function(x) x
  col <- c(col2rgb("olivedrab3")) # orange, darkgoldenrod2
  col <- rgb(col[1], col[2], col[3], alpha = 255 * 0.4, maxColorValue = 255)



  a <- marker_to_assay[[marker]]
  llod <- lloxs[a]



  results <- results[results$comparison == variant,]

  if(monotone) {
    estimates <- results$estimates_mono
  } else {
    estimates <- results$estimates
  }

  standard_errors <- results$se

  thresholds <- results$thresholds
  if(simultaneous_CI) {
    qval <- results$simult_CI_info[[1]]
    index <- which(qval$level == level)
    qval <- as.vector(unlist(qval[index,]$quantile))
  } else {
    qval <- qnorm(1 - (1-level)/2)
  }

  CI_left <- estimates -qval * standard_errors
  CI_right <- estimates +qval* standard_errors
  #print(cbind(CI_left,  estimates, CI_right))
  #print(standard_errors)
  #stop("hi")
  #CI_left <- pmax(CI_left, 0)
  #CI_right <- pmin(CI_right, 0.1)
  #estimates <- pmin(estimates, 0.1)
  plot_data <- data.table(thresholds = thresholds, est= estimates, se= standard_errors,
                          lower = CI_left, upper = CI_right)
  keep <- which(results$thresholds_n_events >= 5)
  plot_data <- plot_data[keep,]

  library(ggplot2)

  # Initial PLOT
  v <- ggplot(plot_data, aes(x = thresholds, y = est)) + geom_point() + geom_line() +
    geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.3, color = NA) + scale_y_continuous()

  # Add RCDF
  max_thresh <- max(data_treated[[marker]][data_treated[[event_type]] == 1])
  # out <- hist(na.omit(data_treated[[marker]]), col=col,axes=F,labels=F,main="",xlab="",ylab="",breaks=10,border=0,freq=F)
  scale_coef <- max(v$data$upper, na.rm = T) * 1

  # coef <- max(out$density)/scale_coef
  RCDF <- function(a) {
    sum(data_treated$wt * (data_treated[[marker]] >= a)) / sum(data_treated$wt) * scale_coef
  }
  RCDF <- Vectorize(RCDF)


  min_range <- min(thresholds)
  max_range <- max(thresholds)
  plot <- v + ggtitle(main) +
    stat_function(fun = RCDF, color = col, geom = "area", fill = col, alpha = 0.2) +
    scale_y_continuous(
      name = laby,
      sec.axis = sec_axis(~ . / scale_coef, name = "Reverse CDF"), n.breaks = 10
    )  +
    theme(plot.title = element_text(size = 25), axis.text.x = element_text(angle = 0, hjust = 1, size = 18), axis.text.y = element_text(angle = 0, hjust = 1, size = 18)) +
    #geom_text(alpha = 0.75, aes(quantile(v$data$thresholds, 0.1),min(max(v$data$upper),risk_plac),label = paste0("placebo overall risk: ", risk_plac)), vjust = 0, size = 5)+
    scale_x_continuous(
      #breaks = xx,
      #labels = do.call(expression,labels),
      name =labx,
      limits = c(min_range, max_range))#,
  #limits = xlimits    )



  append_end <- ""
  append_start <- "PLOT"
  append <- ""
  folder <- ""
  if (monotone) {
    append_start <- paste0(append_start, "_monotone_")
  } else {
    append_start <- paste0(append_start, "_")
  }
  if (simultaneous_CI) {
    append_end <- paste0(append_end, "_", "simultCI",append)
    folder <- "simultaneous_CI"
  } else {
    append_end <- paste0(append_end, "_", "pointwiseCI",append)
    folder <- "pointwise_CI"
  }



  path <- here::here(
    "figs", folder,
    paste0(append_start,  COR, "_", marker, "_", variant, append_end, "_level_", level, ".pdf")
  )
  print(path)
  ggsave(
    filename = path,
    plot = plot, height = 7, width = 9
  )

  return(plot)
}





get_plot_CR_quotient <- function(marker, variant_tgt, variant_ref, simultaneous_CI = F, monotone = F, level = 0.95) {
  print("get_plot_CR_quotient")
  ###### Some variables.
  variant_type <- paste0("seq1.variant.hotdeck", 1) # Used to extract marker data, any hot deck suffices.
  failure_time = "EventTimePrimary" # failure time variable name
  event_type = "EventIndPrimary" # failure type variable name

  library(survival)
  library(data.table)
  # Variant-specific RESULTS for threshold competing risk analysis.
  key <- paste0("output/CompRisk_all_", marker, "_combined.RDS")

  results <- as.data.table(readRDS(here::here(key)))


  # Get data containing marker information for RCDF marker plots.
  event_type_target <- paste0(event_type, "_", variant_type, "_", variant_tgt )
  data_treated <- fread(here::here(paste0("data_clean/data_treated_", marker, "_", failure_time, "_", event_type_target, ".csv")))

  set(data_treated, , event_type ,data_treated[, grep("EventInd", names(data_treated), value = TRUE), with = FALSE])
  set(data_treated, , failure_time ,data_treated[, grep("EventTime", names(data_treated), value = TRUE), with = FALSE])
  data_treated <- na.omit(data_treated)
  data_treated <- as.data.frame(data_treated)
  event_type_target <- as.data.frame(event_type_target)


  key <- marker
  laby <- paste0("Quotients of relative probability:", variant_tgt, " over ", variant_ref)
  labx <- plotting_assay_label_generator(marker, above = TRUE)
  main <- plotting_assay_title_generator(marker)
  ident <- function(x) x
  col <- c(col2rgb("olivedrab3")) # orange, darkgoldenrod2
  col <- rgb(col[1], col[2], col[3], alpha = 255 * 0.4, maxColorValue = 255)



  a <- marker_to_assay[[marker]]
  llod <- lloxs[a]

  print("hey")
  results <- results[results$comparison == paste0(variant_tgt, "_", variant_ref)]
  if(monotone) {
    estimates <- results$estimates_mono
  } else {
    estimates <- results$estimates
  }
  standard_errors <- results$se
  if(simultaneous_CI) {
    qval <- results$simult_CI_info[[1]]
    index <- which(qval$level == level)
    qval <- as.vector(unlist(qval[index,]$quantile))
  } else {
    qval <- qnorm(1 - (1-level)/2)
  }
  thresholds <- results$thresholds
  CI_left <- estimates - qval * standard_errors
  CI_right <- estimates + qval * standard_errors
  #print(cbind(CI_left,  estimates, CI_right))
  #print(standard_errors)
  #stop("hi")
  CI_left <- pmax(CI_left, -2)
  CI_right <- pmin(CI_right, 2)
  #estimates <- pmin(estimates, 0.1)
  plot_data <- data.table(thresholds = thresholds, est= estimates, se= standard_errors,
                          lower = CI_left, upper = CI_right)
  keep <- which(results$thresholds_n_events >= 5)
  plot_data <- plot_data[keep,]

  print(head(results))




  library(ggplot2)

  # Initial PLOT
  v <- ggplot(plot_data, aes(x = thresholds, y = est)) + geom_point() + geom_line() +
    geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.3, color = NA) + scale_y_continuous()

  # Add RCDF
  max_thresh <- max(data_treated[[marker]][data_treated[[event_type]] == 1])
  # out <- hist(na.omit(data_treated[[marker]]), col=col,axes=F,labels=F,main="",xlab="",ylab="",breaks=10,border=0,freq=F)
  scale_coef <- max(v$data$upper, na.rm = T) * 1

  # coef <- max(out$density)/scale_coef
  RCDF <- function(a) {
    sum(data_treated$wt * (data_treated[[marker]] >= a)) / sum(data_treated$wt) * scale_coef
  }
  RCDF <- Vectorize(RCDF)


  min_range <- min(thresholds)
  max_range <- max(thresholds)
  plot <- v + ggtitle(main) +
    stat_function(fun = RCDF, color = col, geom = "area", fill = col, alpha = 0.2) +
    scale_y_continuous(
      name = laby,
      sec.axis = sec_axis(~ . / scale_coef, name = "Reverse CDF"), n.breaks = 10
    )  +
    theme(plot.title = element_text(size = 25), axis.text.x = element_text(angle = 0, hjust = 1, size = 18), axis.text.y = element_text(angle = 0, hjust = 1, size = 18)) +
    #geom_text(alpha = 0.75, aes(quantile(v$data$thresholds, 0.1),min(max(v$data$upper),risk_plac),label = paste0("placebo overall risk: ", risk_plac)), vjust = 0, size = 5)+
    scale_x_continuous(
      #breaks = xx,
      #labels = do.call(expression,labels),
      name =labx,
      limits = c(min_range, max_range))#,
  #limits = xlimits    )



  append_end <- ""
  append_start <- "PLOT_relative"
  append <- ""
  folder <- ""
  if (monotone) {
    append_start <- paste0(append_start, "_monotone_")
  } else {
    append_start <- paste0(append_start, "_")
  }
  if (simultaneous_CI) {
    append_end <- paste0(append_end, "_", "simultCI",append)
    folder <- "simultaneous_CI"
  } else {
    append_end <- paste0(append_end, "_", "pointwiseCI",append)
    folder <- "pointwise_CI"
  }



  path <- here::here(
    "figs", folder,
    paste0(append_start,  COR, "_", marker, "_", variant_tgt, "_", variant_ref, append_end, "_level_", level, ".pdf")
  )
  print(path)
  ggsave(
    filename = path,
    plot = plot, height = 7, width = 9
  )

  return(plot)
}

