#' Apply Method
#'
#' Predicts the values of new data using the output from the train_fd
#' function. Returns the calculated residuals and the plotting statistic.
#'
#'
#' @param model Output from the train_fd function.
#' @param new_data A multivariate time series in dataframe or matrix form.
#' @return A named list including the plotting statistic and residuals.
#' @name predict_fd
#' @rdname predict_fd
#' @export
predict_fd <- function(model, new_data) {

  # Constants
  l <- model$constants[2]

  # Center and Scale
  if(model$center_scale) {
    new_data <-
      purrr::pmap_dfc(list(new_data, model$mean_sd$mean, model$mean_sd$sd), \(x, y, z) (x-y)/z) |>
      as.matrix()
  }

  # Design and Prediction Matrices
  X <- create_X(new_data, lags = l)
  Y <- create_Y(new_data, lags = l)

  # Predictions: GRU
  if(grepl("gru", model$method)) {

    preds <- pred_gru(model$model, X) # python
    colnames(preds) <- colnames(new_data)

    # Predictions: MRF
  } else if(grepl("mrf", model$method)) {
    preds <- randomForestSRC::get.mv.predicted(predict(model$model, as.data.frame(cbind(Y, X))))
    colnames(preds) <- colnames(new_data)

    # Predictions: VARMA
  } else if(grepl("var", model$method)) {
    # Extract VAR(1) Model Info
    Ph0 <- model$model$Ph0
    Phi <- model$model$Phi
    Theta <- matrix(0, nrow = dim(Phi)[1], ncol = dim(Phi)[2])

    # Convert New Data to Matrix
    new_data <- as.matrix(new_data)

    # Prep Residual Matrix
    residuals <- matrix(nrow = nrow(new_data), ncol = ncol(new_data))

    # Calc First Row
    pred <- as.numeric(Ph0+Phi %*% as.numeric(tail(model$model$data, 1))-Theta %*% as.numeric(tail(model$model$residuals, 1)))

    1:ncol(new_data) |>
      purrr::walk(\(j) {
        residuals[1, j] <<- new_data[1, j] - pred[j]
      })

    # Calculate Residuals
    2:nrow(new_data) |>
      purrr::walk(\(i) {

        pred <- as.numeric(Ph0+Phi %*% new_data[i-1, ]-Theta %*% residuals[i-1, ])

        1:ncol(new_data) |>
          purrr::walk(\(j){
            residuals[i,j] <<- new_data[i, j] - pred[j]
          })
      })

    # Predictions
    preds <- Y - residuals[-1, ]

  } else {
    preds <- matrix(0, nrow = nrow(Y), ncol = ncol(Y))
  }


  # Get Residuals
  residuals = Y - preds

  # Get Tau
  tau <- calc_tau(residuals)

  if(grepl("mcusum", model$method)) {
    # Get S
    S <- calc_S(tau, model$constants[1], model$mu_tau, model$sigma_tau_inv)

    # Plotting Statistic
    pstat <- calc_PStat(S, model$sigma_tau_inv)
  } else if(grepl("mewma", model$method)) {
    # Get D
    D <- calc_D(tau, model$mu_tau, model$sigma_tau_inv)

    # Plotting Statistic
    pstat <- calc_PStat_MEWMA(model$constants[1], D, model$constants[3])
  } else if(grepl("htsquare", model$method)) {

    pstat <- calc_D(tau, model$mu_tau, model$sigma_tau_inv)
  }



  # Return
  list(pstat = pstat,
       residuals = residuals)
}


#' Predict Multivariate Random Forest
#'
#' Similar to the build_forest_predict function from MultivariateRandomForest. However,
#' this function takes in output from build_mrf and makes a new prediction.
#'
#'
#' @param model Output from build_mrf.
#' @param testX The design matrix for predictions. Can be made with create_X.
#' @return A prediction matrix
#' @name predict_mrf
#' @rdname predict_mrf
#' @export
predict_mrf <- function (model, testX) {
  # Prediction setup
  Variable_number = model$constants[["p"]]
  n_tree = model$constants[["n_tree"]]

  # Set up prediction matrices
  Y_HAT = matrix(0 * (1:Variable_number * nrow(testX)), ncol = Variable_number, nrow = nrow(testX))
  Y_pred = NULL

  # Make predictions
  1:n_tree |>
   purrr::walk(\(i) {
     # Prediction of single tree
     Y_pred <<- MultivariateRandomForest::single_tree_prediction(model$trees[[i]],
                                                                 testX,
                                                                 Variable_number)

     # Add to other tree predictions
     1:Variable_number |>
       purrr::walk(\(j) {
         Y_HAT[, j] <<- Y_HAT[, j] + Y_pred[, j]
       })
   })

  # Average predictions
  Y_HAT = Y_HAT/n_tree
  return(Y_HAT)
}
