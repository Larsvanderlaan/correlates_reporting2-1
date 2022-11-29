#'
#' Targets the J component of the likelihood for the parameter of interest
#' @param likelihoods The list of likelihood
#' @param fits The list of learner fits the likelihood
#' @param data The data in long format
#' @param target_times The times at which to estimte the cumulative incidence
#' @param node_list A list named list of nodes/variable dictionary
#' @import stats
target_J <- function(likelihoods, fits, data, target_times, node_list) {
  dNt <- as.numeric(data$Nt==1 & data$Event==1)
  Ct_left <- (data$Ct==1) - (data$Ct==1 & data$Event==1)
  nt <- max(data$t)

  surv_C_left <- hazard_to_survival(likelihoods$C, nt, left = T)
  H_list <- list()
  zero_by <-  dNt *(1-Ct_left)

  for(i in 1:ncol(likelihoods$outcomes_A)) {
    gG <- likelihoods$A[[i]]*surv_C_left
    H_i <-likelihoods$outcomes_A[[i]]/ bound(gG, 0.0025)
    H_list[[i]] <-  H_i
  }
  H <- do.call(cbind, H_list)

  H_list <- list()
  for(t_tgt in target_times) {
    ind <- data$t <= t_tgt
    H_list[[as.character(t_tgt)]] <- H*ind
  }
  H <- do.call(cbind, H_list)

  epsilons_J <- list()
  EIC_J <- list()
  if(!is.null(node_list$weights)) {
    weights <- data[[node_list$weights]]
  } else {
    weights <- rep(1, nrow(H))
  }


  for(i in 1:ncol(likelihoods$outcomes_J)) {
    fit_J <- stats::glm(Y~X-1, data = list(Y = likelihoods$outcomes_J[[i]], X = as.matrix(H*zero_by)), family = binomial(), weights = weights*zero_by, offset = stats::qlogis(bound(likelihoods$J[[i]],0.00001)))

    eps <- stats::coef(fit_J)
    eps[is.na(eps)] <- 0
    epsilons_J[[i]] <- eps
    likelihoods$J[[i]] <- as.vector(stats::plogis(stats::qlogis(bound(likelihoods$J[[i]],0.00001)) + as.matrix(H) %*% eps))


    EIC_J[[i]] <- as.matrix(H*zero_by)*(likelihoods$outcomes_J[[i]] - likelihoods$J[[i]])
  }

  fits$epsilons_J <- epsilons_J
  EIC_J <- do.call(cbind, EIC_J)
  EIC_J <- apply(EIC_J, 2, function(v) {
    rowSums(matrix(v, ncol = nt))
  })
  fits$EIC_J <- EIC_J * weights[data$t==1]
  return(list(fits = fits, likelihoods = likelihoods))
}
#'
#' Updates the J component of the likelihood for the parameter of interest.
#' Requires that the likelihood has already been targeted by a call to \code{target_J}
#' @param likelihoods The list of likelihoods
#' @param fits The list of learner fits the likelihood
#' @param data The data in long format
#' @param target_times The times at which to estimte the cumulative incidence
#' @import stats
update_J <- function(likelihoods, fits, data, target_times) {
  epsilons_J <- fits$epsilons_J
  if(is.null(epsilons_J)) {
    stop("J is not targeted yet.")
  }
  nt <- max(data$t)
  Ct_left <- (data$Ct==1) - (data$Ct==1 & data$Event==1)
  surv_C_left <- hazard_to_survival(likelihoods$C, nt, left = T)
  H_list <- list()
  for(i in 1:ncol(likelihoods$outcomes_A)) {
    H_i <-likelihoods$outcomes_A[[i]]/likelihoods$A[[i]]/surv_C_left
    H_list[[i]] <-  H_i
  }
  H <- do.call(cbind, H_list)
  H_list <- list()
  for(t_tgt in target_times) {
    ind <- data$t <= t_tgt
    H_list[[as.character(t_tgt)]] <- H*ind
  }
  H <- do.call(cbind, H_list)
  for(i in 1:ncol(likelihoods$J)) {
    eps <- epsilons_J[[i]]
    likelihoods$J[[i]] <- stats::plogis(stats::qlogis(bound(likelihoods$J[[i]],0.00001)) + as.matrix(H) %*% eps)
  }
  return(likelihoods)
}



#' Performs one iteration of targeting for the N component of the likelihood for the parameter of interest.
#' @param likelihoods The list of likelihoods
#' @param fits The list of learner fits the likelihood
#' @param data The data in long format
#' @param target_times The times at which to estimte the cumulative incidence
#' @param n_full_sample The total sample size (including those with missing A values)
#' @param max_eps Max epsilon value to consider when performing MLE
#' @param force_converge Forces convergence and computes EIF components. This is necessary if for whatever reason there is no convergence before the max number of iteratons.
#' @import stats
target_N <- function(data, likelihoods, fits, node_list, target_times, n_full_sample = NULL, max_eps = NULL, force_converge = FALSE) {
  if(!is.null(fits$max_eps)) {
    max_eps <- fits$max_eps
  }
  nt <- max(data$t)
  n <- (nrow(data)/nt)
  if(is.null(n_full_sample)) {
    n_full_sample <- n
  }
  #print(apply(matrix(as.vector(likelihoods$N), ncol = nt), 2, function(v) {print(quantile(v))}))
  
  #print(quantile(likelihoods$N))
  survs <- compute_survival_functions(likelihoods, nt)
  #print(head(as.matrix(matrix(likelihoods$C, ncol =nt))))
   
  surv_C <- pmax(0.005,hazard_to_survival(likelihoods$C, nt = nt, left = T))
  print("SurvCCCC")
  print(apply(matrix(surv_C,ncol = nt),2,quantile))
  
  
  
  Ft <- survs$Ft
  St <-  survs$St 
  print(colMeans(matrix(St,ncol = nt)))
  
   #stop("hi")
  Ft <- lapply(1:ncol(Ft), function(i) {
    v <- Ft[,i]
    v <- matrix(v, ncol = nt)
    Hv <- list()
    for(k in seq_along(target_times) ){
      t <- target_times[k]
      Hv[[k]] <- (data$t <= t) * (likelihoods$J[[i]] - as.vector(v[,t] - v)/St)

    }
    return(do.call(cbind, Hv))
  })
  Ft <- do.call(cbind, Ft)
  
  H <-  Ft /surv_C
  
  H_list <- list()
  for(i in 1:ncol(likelihoods$outcomes_A)) {
     
    g <- pmax(likelihoods$A[[i]],0.00025)
     
    H_list[[i]] <- H*(likelihoods$outcomes_A[[i]]/g)
  }
  H <- as.matrix(do.call(cbind,H_list))
  H <- bound(H, c(-100,100))
  rm(H_list)
  dNt <- data$Event*data$Nt
  if(!is.null(node_list$weights)) {
    weights <- as.vector(data[[node_list$weights]])
  } else {
    weights <- rep(1, nrow(data))
  }

  D_N <- data$at_risk * H * (dNt - likelihoods$N) * weights
  # Empirical mean of EIF
  direction <- apply(D_N, 2, sum) / n_full_sample

  D_weights <- fits$D_weights
  if(is.null(D_weights)) {

    D_weights <- 1/apply(D_N,2,function(v) {
      res <- sd(c(rep(0, n_full_sample - n), rowSums(matrix(v, ncol = nt))))
      res[res < 1e-6] <- 1e-6
      res
    })
    D_weights <- pmin(D_weights, 1.2*sqrt(n_full_sample)/log(n_full_sample))
    D_weights[is.infinite(D_weights)] <- 0
    D_weights[is.na(D_weights)] <- 0
    fits$D_weights <- D_weights

  }
  #print(D_weights)

  #print(direction*D_weights)
  norm <- sqrt(sum((direction*D_weights)^2)/length(direction))
  print("Current scores: N")
  print(apply(D_N, 2, mean))
  if(norm == 0){
      norm <- 1
  }
  # print(fits$epsilons_N)
  
  #if(all(abs(direction) <= pmax(1/(weights*sqrt(n)*log(n)), 1/n))) {
  if(!is.null(fits$epsilons_N) && length(fits$epsilons_N) > 1 && (all(abs(direction*D_weights) <= 0.05/sqrt(n_full_sample)/log(n_full_sample)) || (norm <= 0.1/(sqrt(n_full_sample)*log(n_full_sample)))) || force_converge) {
    D_N <- apply(D_N,2,function(v) {
      rowSums(matrix(v, ncol = nt))
    })

    fits$EIC_N <- D_N
    fits$Ft <- survs$Ft
    return(list(fits = fits, likelihoods = likelihoods, converged = TRUE))
  }
  #which_converged <- abs(direction*D_weights) <= 0.1/sqrt(n_full_sample)/log(n_full_sample)
  #direction[which_converged] <- 0

  direction <- direction*D_weights
  direction <- direction / sqrt(mean(direction^2))#sqrt(mean(D_weights^2))
  direction[is.na(direction)|is.infinite(direction)] <- 0
  print("HHHHHHHHHHh")
  print(quantile(H))
  H <- bound(H, c(-250,250))
  H_orig <- H

  H <- H %*% direction
  keep <- data$at_risk==1
  lik_N <- bound(likelihoods$N[keep], 0.0000001)
   
  dNt_train <- dNt[keep]
  H_train <- H[keep,]
  if(!is.null(node_list$weights)) {
    weights <- as.vector(data[[node_list$weights]][keep])
  } else {
    weights <- rep(1, length(dNt_train))
  }

  risk <- function(epsilon) {


      update <- bound(plogis(qlogis(lik_N) + as.vector(H_train) * epsilon ),0.000000001)

    loss <- -1 * ifelse(dNt_train == 1, log(update), log(1 - update))
    
    return(stats::weighted.mean(loss, weights))

  }

  optim_fit <- optim(
    par = list(epsilon = max_eps), fn = risk,
    lower = -max_eps, upper = max_eps,
    method = "Brent"
  )
  epsilon <- optim_fit$par
  
  print("EPSILON")
  print(epsilon)
  beta <- coef(glm(Y~ H_train - 1, data = list(Y = dNt_train, H_train = H_train), offset = qlogis(lik_N), family = binomial(), weights = weights))
  print(beta)
  if(is.null(epsilon) || is.na(epsilon) || is.nan(epsilon) || !is.numeric(epsilon)) {
      epsilon <- 0
  }
  #if(abs(epsilon) <= max_eps*0.85) {
   # fits$max_eps <- abs(max_eps)*0.85
  #}
  # if(abs(epsilon) > max_eps*0.975) {
  #   fits$max_eps <- abs(max_eps)*1.3
  # }

  likelihoods$N <- plogis(qlogis(bound(likelihoods$N, 0.0000001)) + as.vector(H) * epsilon )
    
  #print(quantile(H))
  full_epsilon <- epsilon * direction
  fits$epsilons_N <-c(fits$epsilons_N, list(full_epsilon))
  return(list(fits = fits, likelihoods = likelihoods, converged = F))
}

#' Performs one iteration of targeting for the N component of the likelihood for the parameter of interest.
#' Requires that the likelihood has already been targeted by a call to \code{target_N}
#' @param likelihoods The list of likelihoods
#' @param fits The list of learner fits the likelihood
#' @param data The data in long format
#' @param target_times The times at which to estimte the cumulative incidence
#' @param node_list A list named list of nodes/variable dictionary
#' @param step Step number
#' @import stats
update_N <- function(data, likelihoods, fits, node_list, target_times, step) {
  nt <- max(data$t)
  survs <- compute_survival_functions(likelihoods, nt)
  surv_C <- hazard_to_survival(likelihoods$C, nt = nt, left = T)
  Ft <- survs$Ft
  St <- survs$St
  Ft <- lapply(1:ncol(Ft), function(i) {
    v <- Ft[,i]
    v <- matrix(v, ncol = nt)
    Hv <- list()
    for(k in seq_along(target_times) ){
      t <- target_times[k]
      Hv[[k]] <- (data$t <= t) * (likelihoods$J[[i]] - as.vector(v[,t] - v)/St)

    }
    return(do.call(cbind, Hv))
  })
  Ft <- do.call(cbind, Ft)
  H <-  as.matrix(Ft /surv_C)

  H_list <- list()
  for(i in 1:ncol(likelihoods$outcomes_A)) {
    H_list[[i]] <- H*(likelihoods$outcomes_A[[i]]/likelihoods$A[[i]])
  }
  H <- as.matrix(do.call(cbind,H_list))
  H <- bound(H, c(-50,50))
  direction <- as.vector(fits$epsilons_N[[step]])
  likelihoods$N <- plogis(qlogis(likelihoods$N) + H %*% direction)
  likelihoods$Ft <- survs$Ft
  return(likelihoods)
}



