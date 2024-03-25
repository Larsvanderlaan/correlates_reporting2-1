
#' internal use
.target_event_type_distr <- function(target_failure_time, target_treatment, target_event_type, data_long, treatment_levels, event_type_levels, treatment.hat, censor.hazard.hats, total.hazard.hats, event_type.distr.hats, weights) {
  n <- sum(weights)
  time_grid <- 1:max(data_long$failure_time)

  # treatment
  treatment_ind <- 1*(data_long$treatment == target_treatment)
  treatment.a0.hat <- as.vector(treatment.hat[, match(target_treatment, treatment_levels)])


  # censoring
  censor_ind <- 1*(data_long$failure_time >= data_long$t & data_long$event_type != 0)
  S_censor <- long_hazard_to_survival_mats(censor.hazard.hats, time_grid = time_grid, left_cont = TRUE)
  names(S_censor) <- as.character(treatment_levels)
  S_censor_a0 <- as.vector(S_censor[[as.character(target_treatment)]])
  S_censor_a0 <- pmax(S_censor_a0, 25/sqrt(n)/log(n))

  # event_type
  event_type.distr.hat_a0_j0 <- as.vector(event_type.distr.hats[[as.character(target_treatment)]][, match(target_event_type, event_type_levels)])


  # clever covariate
  treatment.a0.hat <- pmax(treatment.a0.hat, 1e-8)
  H1.hat <- (treatment_ind/treatment.a0.hat) * (censor_ind/S_censor_a0)
  dN <- data_long$dN
  ind_j0 <- 1*(data_long$event_type == target_event_type)
  # tmle update


  event_type.distr.hat_a0_j0 <- pmin(pmax(event_type.distr.hat_a0_j0, 1e-8), 1 - 1e-8)
  eps.hat <- coef(glm.fit(rep(1,length(ind_j0)), ind_j0 ,family = binomial(), weights = weights*dN*H1.hat, offset = qlogis(event_type.distr.hat_a0_j0)))
  event_type.distr.hat_a0_j0 <- plogis(qlogis(event_type.distr.hat_a0_j0) + eps.hat)

  score <- mean(weights * dN*H1.hat* (ind_j0 - event_type.distr.hat_a0_j0))
  if(abs(score) > 1e-5) {
    print("TMLE score for event_typeer distribution is solved at low level. There may be convergence issues")
    print(mean(weights * dN*H1.hat* (ind_j0 - event_type.distr.hat_a0_j0)))
  }

  EIF <- weights *dN*H1.hat* (ind_j0 - event_type.distr.hat_a0_j0)
  event_type.distr.hats[[as.character(target_treatment)]][, match(target_event_type, event_type_levels)] <- event_type.distr.hat_a0_j0
  return(list(event_type.distr.hats = event_type.distr.hats, score = score, EIF = EIF))

}
