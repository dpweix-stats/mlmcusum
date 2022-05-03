#' Create tau vector
#'
#' A helper function to create the tau vector given a
#' vector of residuals from one of our methods.
#' 
#'
#' @param residuals The errors from one of our prediction methods on a Y matrix.
#' @return The tau vector to be used for calculating the S vector.
#' @export
calc_tau <- function(residuals) {
  p <- ncol(residuals)
  
  1:nrow(residuals) |> 
    purrr::map(
      \(x) {
        ks::vech(residuals[x, ] %*% t(residuals[x, ]))
      }) |> 
    unlist() |> 
    matrix(ncol = p*(p+1)/2, byrow = TRUE)
}


#' Calculate S matrix
#'
#' A helper function to calculate the S vector given tau, k,
#' mu_tau, and sigma_tau_inv. While k is a tuning parameter, the rest
#' are calculated from the data.
#' 
#'
#' @param tau A vector representing covariance.
#' @param k A tuning parameter.
#' @param mu_tau The mean of tau. Estimated from the training data.
#' @param sigma_tau_inv The inverse covariance matrix of tau. Estimated from the training data.
#' @return A vector S which is used to calculate the plotting statistic.
#' @export
calc_S <- function(tau, k, mu_tau, sigma_tau_inv) {
  N <- nrow(tau)
  p <- ncol(tau)
  
  S <- matrix(0, nrow = N, ncol = ncol(tau))
  C <- vector("numeric", N)
  C[1] <- sqrt(t(tau[1, ] - mu_tau) %*% sigma_tau_inv %*% (tau[1, ] - mu_tau))
  
  if(!any(is.na(tau))) {
    2:N |> 
      purrr::walk(
        \(x) {
          c <- S[x-1, ] + tau[x, ] - mu_tau
          C[x] <<- sqrt(t(c) %*% sigma_tau_inv %*% c)
          
          if(C[x] > k) {
            s <- (S[x-1, ] + tau[x, ] - mu_tau)*(1-k/C[x])
            1:p |> 
              purrr::walk(
                \(y) {
                  S[x, y] <<- s[y]
                })
          }
        })
  }
  
  S
}


#' Calculate D matrix
#'
#' A helper function to calculate the D vector given tau,
#' mu_tau, and sigma_tau_inv. This is the mahalanobis distance for tau.
#' 
#'
#' @param tau A vector representing covariance.
#' @param mu_tau The mean of tau. Estimated from the training data.
#' @param sigma_tau_inv The inverse covariance matrix of tau. Estimated from the training data.
#' @return A vector D which is used to calculate the plotting statistic.
#' @export
calc_D <- function(tau, mu_tau, sigma_tau_inv) {
  N <- nrow(tau)
  p <- ncol(tau)
  
  D <- vector("numeric", length = N)
  
  1:N |> 
    purrr::walk(
      \(i) {
        D[i] <<- t(tau[i, ] - mu_tau) %*% sigma_tau_inv %*% (tau[i, ] - mu_tau)
      })
  
  D
}


#' Calculate plotting statistic
#'
#' A helper function to calculate the plotting statistic given
#' S and sigma_tau_inv.
#' 
#'
#' @param S A matrix which tracks the MCUSUM through observations 
#' @param sigma_tau_inv The inverse covariance matrix of tau. Estimated from the training data.
#' @return The plotting statistic which determines if the process is considered in-control or out-of-control at a given observation.
#' @export
calc_PStat <- function(S, sigma_tau_inv) {
  N <- nrow(S)
  
  1:N |> 
    purrr::map_dbl(
      \(x) {
        sqrt(t(S[x, ]) %*% sigma_tau_inv %*% S[x, ])
      })
}


#' Calculate plotting statistic
#'
#' A helper function to calculate the plotting statistic given
#' D and r.
#' 
#'
#' @param r The weight of the current observation. The previous observation has weight 1-r.
#' @param D The Mahalanobis distance of each tau vector.
#' @param p The number of variables in your multivariate time series.
#' @return The plotting statistic which determines if the process is considered in-control or out-of-control at a given observation.
#' @export
calc_PStat_MEWMA <- function(r, D, p) {
  N <- length(D)
  
  # First Entry
  pstat <- vector(mode = "numeric", length = N)
  pstat[1] <- r*D[1] + (1-r)*(p*(p+1)/2)
  
  # Other Entries
  2:N |> 
    purrr::walk(
      \(i) {
        pstat[i] <<- r*D[i]+(1-r)*pstat[i-1]
      })
  
  pstat
}


#' Create Design Matrix
#'
#' A helper function to create a design matrix to be the 
#' input for the machine learning MCUSUM methods. Requires
#' choosing the number of lags in addition to each variable.
#' Potentially later the number of lags could be a vector
#' of length p. That way the number of lags could be different
#' for each variable.
#' 
#'
#' @param data A dataframe or matrix of observations. 
#' @param lags The number of lags of each variable to be included in the design matrix.
#' @return A design matrix for training/predictions with the machine learning MCUSUM methods.
#' @export
create_X <- function(data, lags) {
  N <- nrow(data)
  p <- ncol(data)
  l <- lags
  
  
  1:(N-l) |> 
    purrr::map(
      \(x) {
        data[x:(x+l-1), 1:p] |> 
          unlist() |> 
          as.numeric()
      }) |> 
    unlist() |> 
    matrix(ncol = l*p, byrow = TRUE)
}


#' Create Prediction Matrix
#'
#' A helper function to create a prediction matrix for 
#' use when training the machine learning MCUSUM methods.
#' Requires choosing the number of lags so that the prediction
#' and design matrices line up.
#' 
#'
#' @param data A dataframe or matrix of observations. 
#' @param lags The number of lags of each variable, the first row of this matrix will be row lag + 1 of the data.
#' @return A prediction matrix which can be used for training/testing the machine learning MCUSUM methods.
#' @export
create_Y <- function(data, lags) {
  N <- nrow(data)
  p <- ncol(data)
  l <- lags
  
  data[(l+1):N, ] |> 
    as.matrix()
}


#' Scale a TS
#'
#' Scale a univariate time series between lower bound a and upper bound b.
#' Used for generating t for a simulation study.
#' 
#'
#' @param ts A univariate time series.
#' @param a Lower bound of time series, default is 0.01.
#' @param b Upper bound of time series, default is 2.
#' @return Returns the scaled time series.
#' @export

scale_t <- function(ts, a = .01, b = 2) {
  as.numeric((b-a)*(ts - min(ts))/(max(ts) - min(ts)) + a)
}


#' Calculate Run Length
#'
#' Calculates the run length for a single instance of a method given a 
#' plotting statistic and an upper limit by comparing the total number
#' of flagged observations to the total number of observations.
#' 
#'
#' @param pstat The plotting statistic from a fault detection method.
#' @param h The upper in-control limit for this plotting statistic.
#' @return Returns the run length.
#' @export
get_run_length <- function(pstat, h) {
  denom <- sum(pstat > h)
  numer <- length(pstat)
  
  ifelse(denom == 0, numer, numer/denom)
}

#' Calculate Run Length
#'
#' Calculates the run length for a single instance of a method given a 
#' plotting statistic and an upper limit by returning the index of the
#' first flagged observation.
#' 
#'
#' @param pstat The plotting statistic from a fault detection method.
#' @param h The upper in-control limit for this plotting statistic.
#' @return Returns the run length.
#' @export
get_first_fault <- function(pstat, h) {
  ret <- which(pstat > h)[1]
  
  if(length(ret) == 0, length(pstat), ret)
}