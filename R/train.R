#' Train methods
#'
#' Train the mrfMCUSUM, gruMCUSUM method. Also works for the resMCUSUM method
#' from Bodnar et al. (2017), which uses a VARMA(1, 1). 
#' 
#'
#' @param data A multivariate time series in dataframe or matrix form. 
#' @param lags The number of lags of each variable to be included in the design matrix.
#' @param k A tuning parameter for the MCUSUM, large k results in shorter memory.
#' @param r A tuning parameter for MEWMA, large r results in shorter memory.
#' @param center_scale A logical, whether or not to center and scale data before modeling.
#' @return A named list including the plotting statistic, trained model, residuals, and constants.
#' @name train
#' @rdname train
#' @export
train <- function(data, method = "gruMCUSUM", lags = 1, k = 1.1, r = .3,
                  center_scale = TRUE) {
  # Constants
  method <- tolower(method)
  l <- ifelse(grepl("varma", method), 1, lags)
  p <- ncol(data)
  mean_sd <- list(mean = colMeans(data),
                  sd = purrr::map_dbl(data, \(x) sd(x)))
  
  # Center and Scale
  if(center_scale) {
    data <- 
      purrr::pmap_dfc(list(data, mean_sd$mean, mean_sd$sd), \(x, y, z) (x - y)/z) |> 
      as.matrix()
  }
  
  # Constant for MCUSUM or MEWMA
  if(grepl("mcusum", method)) {
    constants <- c(k, l, p)
    names(constants) <- c("k", "lags", "p")
  } else if(grepl("mewma", method)) {
    constants <- c(r, l, p)
    names(constants) <- c("r", "lags", "p")
  } else if(grepl("htsquare", method)) {
    constants <- c(0, l, p)
    names(constants) <- c("NA", "lags", "p")
  }
  
  # Design and Prediction Matrices
  X <- create_X(data, lags = l)
  Y <- create_Y(data, lags = l)
  
  # Methods: GRU
  if(grepl("gru", method)) {
    fit <- train_gru(X, Y, l) # python
    
    preds <- pred_gru(fit, X) # python
    colnames(preds) <- colnames(data)
    
  # Methods: MRF
  } else if(grepl("mrf", method)) {
    fit <-
      MultivariateRandomForest::build_single_tree(X, Y,
                                                  m_feature = floor(sqrt(l*p)),
                                                  min_leaf = 10,
                                                  Inv_Cov_Y = solve(cov(Y)),
                                                  Command = 2)
    preds <- MultivariateRandomForest::single_tree_prediction(fit, X, p)
    colnames(preds) <- colnames(data)
  
  # Methods: VARMA
  } else if(grepl("varma", method)) {
    fit <-
      MTS::VAR(data, p = 1, include.mean = TRUE)
    
    preds <- Y - fit$residuals
    colnames(preds) <- colnames(data)
    
  # Methods: Hotelling's T Square
  } else if(grepl("htsquare", method)) {
    fit <- NA
    
    preds <- matrix(0, nrow = nrow(Y), ncol = ncol(Y))
  }
  
  
  # Get Tau
  tau <- calc_tau(Y - preds)
  mu_tau <- colMeans(tau)
  sigma_tau_inv <- solve(cov(tau))
  
  # Calc Pstat
  if(grepl("mcusum", method)) {
    # MCUSUM
    S <- calc_S(tau, constants[1], mu_tau, sigma_tau_inv)
    pstat <- calc_PStat(S, sigma_tau_inv)
  } else if(grepl("mewma", method)) {
    # MEWMA
    D <- calc_D(tau, mu_tau, sigma_tau_inv)
    pstat <- calc_PStat_MEWMA(r, D, p)
  } else if(grepl("htsquare", method)) {
    # Hotelling's T^2
    pstat <- calc_D(tau, mu_tau, sigma_tau_inv)
  }
  
  # Return
  list(pstat = pstat,
       model = fit,
       method = method, 
       residuals = Y - preds,
       center_scale = center_scale,
       mean_sd = mean_sd,
       mu_tau = mu_tau,
       sigma_tau_inv = sigma_tau_inv,
       constants = constants)
}

#' @rdname train
#' @export
train_gruMCUSUM <- function(data, lags = 1, k = .9) {
  l <- lags
  p <- ncol(data)
  
  constants <- c(k, l, p)
  names(constants) <- c("k", "lags", "p")
  
  X <- create_X(data, lags = l)
  Y <- create_Y(data, lags = l)
  
  
  fit_gru <- train_gru(X, Y, l) # python
  
  gru_preds <- pred_gru(fit_gru, X) # python
  colnames(gru_preds) <- colnames(data)

  # Get Tau
  tau <- calc_tau(Y - gru_preds)
  
  mu_tau <- colMeans(tau)
  sigma_tau_inv <- solve(cov(tau))

  # Get S
  S <- calc_S(tau, k, mu_tau, sigma_tau_inv)

  # Plotting Statistic
  pstat <- calc_PStat(S, sigma_tau_inv)
  
  # Return
  list(pstat = pstat,
      model = fit_gru,
      residuals = Y - gru_preds,
      mu_tau = mu_tau,
      sigma_tau_inv = sigma_tau_inv,
      constants = constants)
}

#' @rdname train
#' @export
train_mrfMCUSUM <- function(data, lags = 1, k = .9) {
  l <- lags
  p <- ncol(data)
  
  constants <- c(k, l, p)
  names(constants) <- c("k", "lags", "p")
  
  X <- create_X(data, lags = l)
  Y <- create_Y(data, lags = l)
  
  
  fit_mrf <-
    MultivariateRandomForest::build_single_tree(X, Y,
                                                m_feature = floor(sqrt(l*p)),
                                                min_leaf = 10,
                                                Inv_Cov_Y = solve(cov(Y)),
                                                Command = 2)
  
  mrf_preds <- MultivariateRandomForest::single_tree_prediction(fit_mrf, X, p)
  colnames(mrf_preds) <- colnames(data)
  
  # Get Tau
  tau <- calc_tau(Y - mrf_preds)
  
  mu_tau <- colMeans(tau)
  sigma_tau_inv <- solve(cov(tau))
  
  # Get S
  S <- calc_S(tau, k, mu_tau, sigma_tau_inv)
  
  # Plotting Statistic
  pstat <- calc_PStat(S, sigma_tau_inv)
  
  # Return
  list(pstat = pstat,
       model = fit_mrf,
       residuals = Y - mrf_preds,
       mu_tau = mu_tau,
       sigma_tau_inv = sigma_tau_inv,
       constants = constants)
}


#' @rdname train
#' @export
train_varmaMCUSUM <- function(data, k = 0.9) {
  l <- 1
  p <- ncol(data)
  
  constants <- c(k, l, p)
  names(constants) <- c("k", "lags", "p")
  
  # Fit Model
  fit_varma <-
    MTS::VARMACpp(data, p = 1, q = 1, include.mean = TRUE)
  
  residuals <- fit_varma$residuals
  colnames(residuals) <- colnames(data)
  
  # Get Tau
  tau <- calc_tau(residuals)
  
  mu_tau <- colMeans(tau)
  sigma_tau_inv <- solve(cov(tau))
  
  # Get S
  S <- calc_S(tau, k, mu_tau, sigma_tau_inv)
  
  # Plotting Statistic
  pstat <- calc_PStat(S, sigma_tau_inv)
  
  # Return
  list(pstat = pstat,
       model = fit_varma,
       residuals = residuals,
       mu_tau = mu_tau,
       sigma_tau_inv = sigma_tau_inv,
       constants = constants)
}


#' @rdname train
#' @export
train_varmaMEWMA <- function(data, r = 0.3) {
  l <- 1
  p <- ncol(data)
  
  constants <- c(r, l, p)
  names(constants) <- c("r", "lags", "p")
  
  fit_varma <-
    MTS::VARMACpp(data, p = 1, q = 1, include.mean = TRUE)
  
  residuals <- fit_varma$residuals
  colnames(residuals) <- colnames(data)
  
  # Get Tau
  tau <- calc_tau(residuals)
  
  mu_tau <- colMeans(tau)
  sigma_tau_inv <- solve(cov(tau))
  
  # Get D
  D <- calc_D(tau, mu_tau, sigma_tau_inv)
  
  # Plotting Statistic
  pstat <- calc_PStat_MEWMA(r, D, p)
  
  # Return
  list(pstat = pstat,
       D = D,
       model = fit_varma,
       residuals = residuals,
       mu_tau = mu_tau,
       sigma_tau_inv = sigma_tau_inv,
       constants = constants)
}
