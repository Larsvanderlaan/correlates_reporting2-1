
#' Generates plot of threshold-response function.
#' @param marker The marker variable to generate plots for.
#' @param simultaneous_CI True if simultaneous CI should be plotted. Otherwise if False pointwise CI are plotted.
#' @param monotone True if monotone correction should be done on estimates. Otherwise no monotone correction. This should be done if one expects monotonicity.
#' Assume monotone nonincreasing


plot_threshold_response <- function(results, simultaneous_CI = FALSE, monotone = FALSE) {
   if(monotone) {
     estimates <- results$estimates_monotone
   } else {
     estimates <- results$estimates
   }
  standard_errors <- results$se

  thresholds <- results$threshold
  CI_left <- estimates - qnorm(0.975) * standard_errors
  CI_right <- estimates + qnorm(0.975) * standard_errors
  CI_left <- pmax(CI_left, 0)
  CI_right <- pmin(CI_right, 0.02)
  estimates <- pmin(estimates, 0.02)
  plot_data <- data.table(thresholds = thresholds, est= estimates, se= standard_errors,
                          lower = CI_left, upper = CI_right)

  library(ggplot2)
  plot <- ggplot(plot_data, aes(x = thresholds, y = est)) + geom_point() + geom_line() +
    geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.3, color = NA) + scale_y_continuous()
  return(plot)

}



get_plot <- function(marker, failure_time = "EventTimePrimary", event_type = "EventIndPrimary", variant = "", simultaneous_CI = F, monotone = F, above = TRUE) {
  variant_type <- paste0("seq1.variant.hotdeck", 1)
  library(survival)
  library(data.table)


  key <- paste0("output/vaccine_", marker, "_", failure_time, "_", event_type, "_", variant, ".csv")
  results <- fread(here::here(key))
  event_type_target <- paste0(event_type, "_", variant_type, "_", variant )
  data_treated <- fread(here::here(paste0("data_clean/data_treated_", marker, "_", failure_time, "_", event_type_target, ".csv")))
  print(names(data_treated))
  set(data_treated, , event_type ,data_treated[, grep("EventInd", names(data_treated), value = TRUE), with = FALSE])
  set(data_treated, , failure_time ,data_treated[, grep("EventTime", names(data_treated), value = TRUE), with = FALSE])
  data_treated <- na.omit(data_treated[, c(marker, failure_time, event_type, "wt"), with = F])
  risk_plac <- mean(results$estimates_placebo) # this is a constant
  data_treated <- as.data.frame(data_treated)
  event_type_target <- as.data.frame(event_type_target)
  print(names(data_treated))

  max_t <- unique(results$times)[1]
  key <- marker

  time <- tpeak
  day <- ""
  if(TRIAL == "hvtn705"){
    laby <- paste0("Probability of HIV by Day ",max_t)
  } else {
    laby <- paste0("Probability of COVID by Day ",max_t)
  }

  labx <- plotting_assay_label_generator(marker, above)
  # subtitle_main <- "Nonparametric estimate of Threshold-response function"
  main <- plotting_assay_title_generator(marker)
  if(length(grep("start", key)) > 0) {
    main <- paste0(main , " (1-day-post)")
  }
  ident <- function(x) x
  col <- c(col2rgb("olivedrab3")) # orange, darkgoldenrod2
  col <- rgb(col[1], col[2], col[3], alpha = 255 * 0.4, maxColorValue = 255)
  # Get initial threshold-response plot with simultaneous CI
  v <- plot_threshold_response(results, simultaneous_CI = simultaneous_CI, monotone = F)
  max_thresh <- max(data_treated[[marker]][data_treated[[event_type]] == 1])
  # out <- hist(na.omit(data_treated[[marker]]), col=col,axes=F,labels=F,main="",xlab="",ylab="",breaks=10,border=0,freq=F)
  scale_coef <- max(v$data$upper, na.rm = T) * 1

  # coef <- max(out$density)/scale_coef
  RCDF <- function(a) {
    sum(data_treated$wt * (data_treated[[marker]] >= a)) / sum(data_treated$wt) * scale_coef
  }
  RCDF <- Vectorize(RCDF)
  if (!simultaneous_CI) {
    #main <- paste0(main, " with point-wise confidence intervals")
  } else {
    #main <- paste0(main, " with simultaneous confidence bands")
  }
  a <- marker_to_assay[[marker]]
  #xlimits <- get.range.cor(data, a,   tpeak)
  llod <- lloxs[a] #llox_labels
  #xx <- report.assay.values(data_treated[[marker]], assay)
  #labels_info <- draw.x.axis.cor(xlimits, lloxs[a], llox_labels[a], for.ggplot=T)
  #xx <- labels_info$ticks
  #labels <- as.list(labels_info$labels)


  plot <- v + ggtitle(main) +
    stat_function(fun = RCDF, color = col, geom = "area", fill = col, alpha = 0.2) +
    scale_y_continuous(
      name = laby,
      sec.axis = sec_axis(~ . / scale_coef, name = "Reverse CDF"), n.breaks = 10
    )  +
    theme(plot.title = element_text(size = 25), axis.text.x = element_text(angle = 0, hjust = 1, size = 18), axis.text.y = element_text(angle = 0, hjust = 1, size = 18)) +
    geom_text(alpha = 0.75, aes(quantile(v$data$thresholds, 0.1),min(max(v$data$upper),risk_plac),label = paste0("placebo overall risk: ", risk_plac)), vjust = 0, size = 5)+ scale_x_continuous(
      #breaks = xx,
      #labels = do.call(expression,labels),
      name =labx)#,
      #limits = xlimits    )

  try({
  if(above  && max_thresh < log10(uloqs[a]) - 0.05) {
    plot <- plot + geom_vline(xintercept = max_thresh, colour = "red", linetype = "longdash")
  } else if(!above && risk_plac<= max(v$data$upper, na.rm = T)) {
    plot <- plot + geom_hline(aes(yintercept=risk_plac), alpha = 0.4, linetype = "longdash", colour = "red")
  }
  })
  plot <- plot + ggtitle(paste0(COR, "_", marker, "_", variant))

  #data_interp <- as.data.frame(approx(plot$data$thresholds, plot$data$est, xout = data_treated[data_treated[[event_type]]==1,marker], rule = 2  ))


  #plot <- plot  + geom_point(data=data_interp,aes(x=x, y=y), colour = "blue")

  #+  geom_text(aes(x=max_thresh *(1.01), label="No observed events", y=0.002), colour="black", angle=90, text=element_text(size=11))
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


  print(here::here(
    "figs", folder,
    paste0(append_start,  COR, "_", marker, "_", variant, append_end, ".pdf")
  ))

  path <- here::here(
    "figs", folder,
    paste0(append_start,  COR, "_", marker, "_", variant, append_end, ".pdf")
  )
  print(path)
  print(plot)
  ggsave(
    filename = path,
    plot = plot, height = 7, width = 9
  )

  return(plot)
}
