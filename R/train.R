#' Train methods
#'
#' Train the mrfMCUSUM, gruMCUSUM method. Also works for the resMCUSUM method
#' from Bodnar et al. (2017), which uses a VARMA(1, 1).
#'
#'
#' @param data A multivariate time series in dataframe or matrix form.
#' @param method An indicator of which model and fault detection method to use. Options are gruMCUSUM, mrfMCUSUM, varMCUSUM, varMEWMA, or htsquare.
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
  l <- ifelse(grepl("var", method), 1, lags)
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

    print(c("GRU:", center_scale))

  # Methods: MRF
  } else if(grepl("mrf", method)) {
    fit <- randomForestSRC::rfsrc(formula = as.formula(paste0("Multivar(", paste0(colnames(Y), collapse = ","), ") ~ .")),
                                  data = as.data.frame(cbind(Y, X)),
                                  ntree = 500,
                                  splitrule = "mahalanobis")


    preds <- randomForestSRC::get.mv.predicted(fit)
    colnames(preds) <- colnames(data)

    print(c("MRF:", center_scale))

  # Methods: VARMA
  } else if(grepl("var", method)) {
    fit <-
      MTS::VAR(data, p = 1, include.mean = TRUE)

    preds <- Y - fit$residuals
    colnames(preds) <- colnames(data)

    print(c("VAR:", center_scale))

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

#' Train Multivariate Random Forest
#'
#' Similar to the build_forest_predict function from MultivariateRandomForest. However,
#' this function will save the model for future use.
#'
#'
#' @param trainX The design matrix for predictions. Can be made with create_X.
#' @param trainY The value of the response variables.
#' @param n_tree Number of trees in the forest, which must be positive integer.
#' @param m_feature Number of randomly selected features considered for a split in each regression tree node, which must be positive integer and less than N (number of input features)
#' @param min_leaf Minimum number of samples in the leaf node. If a node has less than or equal to min_leaf samples, then there will be no splitting in that node and this node will be considered as a leaf node. Valid input is positive integer, which is less than or equal to M (number of training samples)
#' @return A named list including the trained model and it's predictions.
#' @name build_mrf
#' @rdname build_mrf
#' @export
build_mrf <- function (trainX, trainY, n_tree, m_feature, min_leaf) {
  # Stopping conditions
  if (class(n_tree) == "character" || n_tree%%1 != 0 ||
      n_tree < 1)
    stop("Number of trees in the forest can not be fractional or negative integer or string")
  if (class(m_feature) == "character" || m_feature%%1 !=
      0 || m_feature < 1)
    stop("Number of randomly selected features considered for a split can not be fractional or negative integer or string")
  if (class(min_leaf) == "character" || min_leaf%%1 !=
      0 || min_leaf < 1 || min_leaf > nrow(trainX))
    stop("Minimum leaf number can not be fractional or negative integer or string or greater than number of samples")

  # Bootstrap setup
  theta <- function(trainX) {
    trainX
  }
  results <- bootstrap::bootstrap(1:nrow(trainX), n_tree, theta)
  b = results$thetastar
  Variable_number = ncol(trainY)

  # MRF or plain RF
  if (Variable_number > 1) {
    Command = 2
  }
  else if (Variable_number == 1) {
    Command = 1
  }

  # Set up prediction matrices
  Y_HAT = matrix(0 * (1:Variable_number * nrow(trainX)), ncol = Variable_number,
                 nrow = nrow(trainX))
  Y_pred = NULL

  # Build single trees
  single_trees <- 1:n_tree |>
    purrr::map(\(i) {
      Single_Model = NULL
      X = trainX[b[, i], ]
      Y = matrix(trainY[b[, i], ], ncol = Variable_number)
      Inv_Cov_Y = solve(cov(Y))

      if (Command == 1) {
        Inv_Cov_Y = matrix(rep(0, 4), ncol = 2)
      }

      Single_Model = MultivariateRandomForest::build_single_tree(X, Y, m_feature, min_leaf, Inv_Cov_Y, Command)
      Y_pred <<- MultivariateRandomForest::single_tree_prediction(Single_Model, trainX, Variable_number)

      1:Variable_number |>
        purrr::walk(\(j) {
          Y_HAT[, j] <<- Y_HAT[, j] + Y_pred[, j]
        })

      print(paste0("Fit tree ", i))

      Single_Model
    })

  # Divide by predictions to get average
  Y_HAT = Y_HAT/n_tree

  # Return value (for use in predict_mrf)
  list(preds = Y_HAT,
       trees = single_trees,
       constants = c(n_tree = n_tree,
                     m_feature = m_feature,
                     min_leaf = min_leaf,
                     p = Variable_number))

}
