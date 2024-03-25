#' Delta Method Custom Function
#'
#' This function provides inference for functions of the survival estimates outputted by `survinf` using the delta method.
#'
#' @param output The output object from a survival analysis using `survinf`.
#' @param transform A function that transforms the treatment and control estimates. Default is the difference between treatment and control estimates.
#' @param transform_grad A function that computes the gradient of the transformation function. Default computes the gradient for difference in estimates.
#' @param control_level The level of the control group in the treatment vector. Default is the minimum value of the treatment vector.
#' @param treatment_levels The levels of the treatment groups to perform inference for. Default is all unique treatment levels except the control level.
#' @param name An optional name for the transformation being applied. Default is "custom".
#'
#' @return A modified output object with transformed estimates, standard errors, confidence intervals, and EIF (Efficient Influence Function) for the specified transformation.
#' @export
delta_method_contrast <- function(output,
                                transform = function(est_treatment, est_control) {est_treatment - est_control},
                                transform_grad = function(est_treatment, est_control) {c(1, -1)}, control_level = min(output$treatment),
                                treatment_levels = setdiff(output$treatment, control_level), name = "custom") {
  treatment_levels <- sort(unique(treatment_levels))
  control_level <- control_level[1]

  event_type_levels <-  unique(output$event_type)
  target_ftime <- sort(unique(unlist(output$times)))
  new_output <- as.data.table(expand.grid(treatment_levels, event_type_levels))
  names(new_output) <- c("treatment", "event_type")
  new_output$name <- name
  new_output$control <- control_level
  new_output$times <- list(sort(unique(unlist(output$times))))
  new_output$estimates = list()
  new_output$se = list()
  new_output$CI = list()
  new_output$EIF = list()
  tmp <- lapply(treatment_levels, function(level){
    lapply(event_type_levels, function(type){

      # get treatment info
      out_treatment <- output[output$treatment==level & output$event_type==type,]
      est_treatment <- unlist(out_treatment$estimates)
      EIF_treatment <- out_treatment$EIF[[1]]
      # get control info
      out_control <- output[output$treatment==control_level & output$event_type==type,]
      est_control <- unlist(out_control$estimates)
      EIF_control <- out_control$EIF[[1]]
      # compute transformed estimates and EIF
      est_transformed <- as.vector(t(mapply(FUN = transform, est_treatment, est_control)))
      names(est_transformed) <- names(output$estimates)
      EIF_transformed <- as.matrix(do.call(cbind, lapply(1:ncol(EIF_treatment), function(index){
        cbind(EIF_treatment[,index], EIF_control[,index]) %*% mapply(FUN = transform_grad, est_treatment[index], est_control[index])
      })))
      row_index <- new_output$event_type==type & new_output$treatment==level
      new_output[row_index, "estimates"] <<- list(est_transformed)
      # new standard errors

      ses <- apply(EIF_transformed, 2, sd)/sqrt(nrow(EIF_transformed))
      # new 95% confidence intervals
      CI <- cbind(est_transformed - qnorm(0.975) * ses, est_transformed + qnorm(0.975) * ses)
      colnames(CI) <- c("CI_left", "CI_right")
      rownames(CI) <- paste0("time_", target_ftime)
      new_output[row_index, "se"] <<- list(ses)
      new_output[row_index, "CI"] <<- list(CI)
      new_output[row_index, "EIF"] <<- list(EIF_transformed)

      return(NULL)
    } )
  })
  return(new_output)

}

#' Delta Method for Average Treatment Effect (ATE)
#'
#' This function provides inference for the Average Treatment Effect (ATE) using the delta method.
#'
#' @param output The output object from a survival analysis using `survinf`.
#' @param control_level The level of the control group in the treatment vector.
#' @param treatment_levels The levels of the treatment groups to perform inference for. Default is all unique treatment levels except the control level.
#'
#' @return A modified output object with transformed estimates, standard errors, confidence intervals, and EIF for the ATE.
#' @export
delta_method_ATE <- function(output, control_level = min(output$treatment), treatment_levels = setdiff(output$treatment, control_level)) {
  transform <- function(est_treatment, est_control) {est_treatment - est_control}
  transform_grad <- function(est_treatment, est_control) {c(1, -1)}
  return(delta_method_contrast(output,
                             transform = transform,
                             transform_grad = transform_grad,
                             control_level = control_level,
                             treatment_levels = treatment_levels,
                             name = "ATE"))

}


#' Delta Method for Log Relative Risk (LRR)
#'
#' This function provides inference for the Log Relative Risk (LRR) using the delta method.
#'
#' @param output The output object from a survival analysis using `survinf`.
#' @param control_level The level of the control group in the treatment vector.
#' @param treatment_levels The levels of the treatment groups to perform inference for. Default is all unique treatment levels except the control level.
#' @param exp_transform Logical, indicating whether to exponentiate the transformed estimates and CI. Default is TRUE.
#'
#' @return A modified output object with transformed estimates, standard errors, confidence intervals, and EIF for the LRR.
#' @export
delta_method_LRR <- function(output, control_level = min(output$treatment), treatment_levels = setdiff(output$treatment, control_level), exp_transform = TRUE) {
  transform <- function(est_treatment, est_control) {log(est_treatment) - log(est_control)}
  transform_grad <- function(est_treatment, est_control) {c(1/est_treatment,  -1/est_control)}
  output_transformed <- delta_method_contrast(output,
                                            transform = transform,
                                            transform_grad = transform_grad,
                                            control_level = control_level,
                                            treatment_levels = treatment_levels,
                                            name = "LRR")
  if(exp_transform) {
    output_transformed$estimates <- list(exp(unlist(output_transformed$estimates)))
    output_transformed$CI <- lapply(output_transformed$CI, function(item){
      list(exp(item))
    })
  }
  return(output_transformed)

}

