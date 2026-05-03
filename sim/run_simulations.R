library(survival)
library(parallel)
library(flexsurv)
library(furrr)
library(aftgee)

set.seed(1)
NSIM <- 3000
N_CORES <- max(1, availableCores() - 1)
plan(multisession, workers = N_CORES)

scenarios <- expand.grid(
  dgp = c("weibull", "loglogistic"),
  n   = c(50, 200),
  cens_rate = c(0.2, 0.5, 0.8),
  rho = c(1.0, 1.25)
)

ci_multiplicative <- function(data,
                              log_lower = -3,
                              log_upper = 3,
                              alpha = 0.05,
                              ngrid = 1001) {
  grid_rho <- exp(seq(log_lower, log_upper, length.out = ngrid))
  pval <- vapply(grid_rho, function(rho) {
    d <- data
    d$Y[d$Z == 0] <- d$Y[d$Z == 0] * rho
    survdiff(Surv(Y, event) ~ Z, data = d, rho = 0)$pvalue
  }, numeric(1))
  data.frame(
    est   = grid_rho[which.max(pval)],
    lower = min(grid_rho[pval > alpha]),
    upper = max(grid_rho[pval > alpha])
  )
}

ci_semi_aft <- function(data, alpha = 0.05) {
  fit <- 
    tryCatch(
      aftgee(
        Surv(Y, event) ~ Z,
        data = data
      ),
      error = function(e)
        NULL
    )
  if (is.null(fit))
    return(data.frame(
      est = NA_real_,
      lower = NA_real_,
      upper = NA_real_
    ))
  b  <- coef(fit)["Z"]
  se <- sqrt(vcov(fit)[2, 2])
  z  <- qnorm(1 - alpha / 2)
  data.frame(
    est   = exp(b),
    lower = exp(b - z * se),
    upper = exp(b + z * se)
  )
}
ci_weibull_aft <- function(data, alpha = 0.05) {
  fit <-
    tryCatch(
      survreg(
        Surv(Y, event) ~ Z,
        data = data,
        dist = "weibull",
        control = survreg.control(maxiter = 100)
      ),
      error = function(e)
        NULL
    )
  if (is.null(fit))
    return(data.frame(
      est = NA_real_,
      lower = NA_real_,
      upper = NA_real_
    ))
  b  <- coef(fit)["Z"]
  se <- sqrt(vcov(fit)["Z", "Z"])
  z  <- qnorm(1 - alpha / 2)
  data.frame(
    est   = exp(b),
    lower = exp(b - z * se),
    upper = exp(b + z * se)
  )
}


find_theta <- function(shape, scale, rho, cens_rate) {
  p_event <- function(nu, theta) {
    x <- (theta / nu) ^ shape
    (nu / (shape * theta)) * pgamma(x, 1 / shape) * gamma(1 / shape)
  }
  uniroot(function(theta) {
    0.5 * p_event(scale, theta) +
      0.5 * p_event(scale * rho, theta) - cens_rate
  }, interval = c(1e-6, 1e6 * scale))$root
}

find_theta_loglogistic <-
  function(shape, scale, rho, cens_rate, n = 1e6) {
    T0 <- rllogis(n, shape = shape, scale = scale)
    Z  <- rbinom(n, 1, 0.5)
    Tt <- ifelse(Z == 1, T0 * rho, T0)
    U  <- runif(n)
    uniroot(function(theta)
      mean(Tt > U * theta) - cens_rate,
      interval = c(1e-6, 1e6 * scale))$root
  }

scenarios$theta <- mapply(function(dgp, rho, cens_rate) {
  if (dgp == "weibull") {
    find_theta(
      shape = 1.5,
      scale = 100,
      rho = rho,
      cens_rate = cens_rate
    )
  } else {
    find_theta_loglogistic(
      shape = 4,
      scale = 100,
      rho = rho,
      cens_rate = cens_rate
    )
  }
}, scenarios$dgp, scenarios$rho, scenarios$cens_rate)

gen_weibull <-
  function(n,
           shape = 1.5,
           scale = 100,
           rho = 1.25,
           theta) {
    Z  <- rbinom(n, 1, 0.5)
    T0 <- rweibull(n, shape = shape, scale = scale)
    Tt <- ifelse(Z == 1, T0 * rho, T0)
    C  <- runif(n, 0, theta)
    data.frame(Y = pmin(Tt, C),
               Z = Z,
               event = as.integer(Tt <= C))
  }

gen_loglogistic <-
  function(n,
           shape = 4,
           scale = 100,
           rho = 1.25,
           theta) {
    Z  <- rbinom(n, 1, 0.5)
    T0 <- rllogis(n, shape = shape, scale = scale)
    Tt <- ifelse(Z == 1, T0 * rho, T0)
    C  <- runif(n, 0, theta)
    data.frame(Y = pmin(Tt, C),
               Z = Z,
               event = as.integer(Tt <= C))
  }

covers <- function(lo, hi, val) {
  ! is.na(lo) && lo <= val && val <= hi
}

run_one <- function(dat, true_rho) {
  if (sum(dat$event) == 0) return(NULL)
  
  np   <- ci_multiplicative(dat)
  waft <- ci_weibull_aft(dat)
  semi <- ci_semi_aft(dat)
  
  data.frame(
    cens_rate    = 1 - mean(dat$event),
    np_covered   = covers(np$lower, np$upper, true_rho),
    np_width     = np$upper - np$lower,
    np_est       = np$est,
    waft_covered = covers(waft$lower, waft$upper, true_rho),
    waft_width   = waft$upper - waft$lower,
    waft_est     = waft$est,
    semi_covered = covers(semi$lower, semi$upper, true_rho),
    semi_width   = semi$upper - semi$lower,
    semi_est     = semi$est
  )
}

dgp_list <- list(
  weibull = gen_weibull, 
  loglogistic = gen_loglogistic)

results <-
  do.call(rbind, lapply(seq_len(nrow(scenarios)), function(s) {
    dgp_name <- scenarios$dgp[s]
    n <- scenarios$n[s]
    rho <- scenarios$rho[s]
    theta <- scenarios$theta[s]
    target_cens <- scenarios$cens_rate[s]
    dgp <- dgp_list[[dgp_name]]
    
    out <- future_map(
      seq_len(NSIM),
      function(i) run_one(dgp(n = n, rho = rho, theta = theta), rho),
      .options = furrr_options(seed = TRUE)   
    ) 
    
    out <- Filter(Negate(is.null), out)
    sims <- do.call(rbind, out)
    cat("Finished iteration", s, "\n")
    
    data.frame(
      dgp          = dgp_name,
      n            = n,
      rho          = rho,
      target_cens  = target_cens,
      n_sims       = nrow(sims),
      cens_rate    = round(mean(sims$cens_rate,    na.rm = TRUE), 3),
      np_coverage  = round(mean(sims$np_covered,   na.rm = TRUE), 3),
      np_width     = round(median(sims$np_width,   na.rm = TRUE), 3),
      np_est       = round(median(sims$np_est,     na.rm = TRUE), 3),
      aft_coverage = round(mean(sims$waft_covered, na.rm = TRUE), 3),
      aft_width    = round(median(sims$waft_width, na.rm = TRUE), 3),
      aft_est      = round(median(sims$waft_est,   na.rm = TRUE), 3),
      semi_coverage = round(mean(sims$semi_covered, na.rm = TRUE), 3),
      semi_width    = round(median(sims$semi_width, na.rm = TRUE), 3),
      semi_est      = round(median(sims$semi_est,   na.rm = TRUE), 3)
    )
  }))

plan(sequential) 
save(results, file = "sim/results/results.rda")
