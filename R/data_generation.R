#' Generate Sample Data
#'
#' Generates a three dimensional multivariate time series for use in
#' simulation studies. Multivariate time series data are generated with
#' either linear or non-linear relationships between variables. Further,
#' the option to produce nonstationary data by adding a sine wave is included.
#' \n Each time a multivariate time series is generated four data frames are
#' produced. There are three different faults that are be introduced, and one
#' data frame in which no fault is introduced.
#'
#'
#' @param n_ic The number of observations which are in-control.
#' @param n_oc The number of observations which are out-of-control, i.e. faulty.
#' @param phi The strength of the autocorrelation in the base time series
#' @param sin_scale The scale of the sine wave for non-stationarity.
#' @param a_1 Increase in variance of x2 for fault 1.
#' @param a_2 Increase in variance of x1, x2, and x3 for fault 3.
#' @return A list of four data frames where the first n_ic observations are
#'  identical, but the last n_oc observations have experienced fault 1, 2, 3,
#'  or no fault.
#' @export
gen_dat_lin <- function(n_ic = 1500, n_oc = 500, phi = .8,
                        a_1 = 3, a_2 = 2) {
  # Create t
  t <- arima.sim(model = list(ar = phi), n = (n_ic + n_oc)) |> scale_t()

  index_ic <- 1:n_ic
  index_oc <- (n_ic+1):(n_ic+n_oc)

  # Generate Training Observations
  dat <-
    tibble::tibble(x1 = t[index_ic] + rnorm(n_ic, 0, .1),
                   x2 = 2*t[index_ic] + rnorm(n_ic, 0, .1),
                   x3 = -.5*t[index_ic] + rnorm(n_ic, 0, .1))

  # Generate Testing Observations
  dat_nf <-
    dplyr::bind_rows(dat,
                     tibble::tibble(x1 = t[index_oc] + rnorm(n_oc, 0, .1),
                                    x2 = 2*t[index_oc] + rnorm(n_oc, 0, .1),
                                    x3 = -.5*t[index_oc] + rnorm(n_oc, 0, .1))
    )

  dat_f1 <-
    dplyr::bind_rows(dat,
                     tibble::tibble(x1 = t[index_oc] + rnorm(n_oc, 0, .1),
                                    x2 = 2*t[index_oc] + rnorm(n_oc, 0, a_1*.1), # increased variance
                                    x3 = -.5*t[index_oc] + rnorm(n_oc, 0, .1))
    )
  dat_f2 <-
    dplyr::bind_rows(dat,
                     tibble::tibble(x1 = t[index_oc] + rnorm(n_oc, 0, .1),
                                    x2 = 2*t[index_oc] + rnorm(n_oc, 0, .1),
                                    x3 = -1*t[index_oc] + 0.5 + rnorm(n_oc, 0, .1)) # alter x3
    )

  dat_f3 <-
    dplyr::bind_rows(dat,
                     tibble::tibble(x1 = t[index_oc] + rnorm(n_oc, 0, a_2*.1), # increased variance
                                    x2 = 2*t[index_oc] + rnorm(n_oc, 0, a_2*.1), # for all
                                    x3 = -.5*t[index_oc] + rnorm(n_oc, 0, a_2*.1)) # variables
    )

  # Return
  list(none = dat_nf,
       f1 = dat_f1,
       f2 = dat_f2,
       f3 = dat_f3)
}

#' @rdname gen_dat_lin
#' @export
gen_dat_ltl <- function(n_ic = 1500, n_oc = 500, phi = .8, sin_scale = 1.5,
                        a_1 = 3, a_2 = 2) {

  # Create t
  t <- (arima.sim(model = list(ar = phi), n = (n_ic + n_oc)) +
                   rep(sin_scale*sin(seq(0, 2*pi, length.out = 51))[-1], (n_ic+n_oc)/50)) |> scale_t()

  index_ic <- 1:n_ic
  index_oc <- (n_ic+1):(n_ic+n_oc)

  # Generate Training Observations
  dat <-
    tibble::tibble(x1 = t[index_ic] + rnorm(n_ic, 0, .1),
                   x2 = 2*t[index_ic] + rnorm(n_ic, 0, .1),
                   x3 = -.5*t[index_ic] + rnorm(n_ic, 0, .1))

  # Generate Testing Observations
  dat_nf <-
    dplyr::bind_rows(dat,
                     tibble::tibble(x1 = t[index_oc] + rnorm(n_oc, 0, .1),
                                    x2 = 2*t[index_oc] + rnorm(n_oc, 0, .1),
                                    x3 = -.5*t[index_oc] + rnorm(n_oc, 0, .1))
    )

  dat_f1 <-
    dplyr::bind_rows(dat,
                     tibble::tibble(x1 = t[index_oc] + rnorm(n_oc, 0, .1),
                                    x2 = 2*t[index_oc] + rnorm(n_oc, 0, a_1*.1), # increased variance
                                    x3 = -.5*t[index_oc] + rnorm(n_oc, 0, .1))
    )
  dat_f2 <-
    dplyr::bind_rows(dat,
                     tibble::tibble(x1 = t[index_oc] + rnorm(n_oc, 0, .1),
                                    x2 = 2*t[index_oc] + rnorm(n_oc, 0, .1),
                                    x3 = -1*t[index_oc] + 0.5 + rnorm(n_oc, 0, .1)) # alter x3
    )

  dat_f3 <-
    dplyr::bind_rows(dat,
                     tibble::tibble(x1 = t[index_oc] + rnorm(n_oc, 0, a_2*.1), # increased variance
                                    x2 = 2*t[index_oc] + rnorm(n_oc, 0, a_2*.1), # for all
                                    x3 = -.5*t[index_oc] + rnorm(n_oc, 0, a_2*.1)) # variables
    )

  # Return
  list(none = dat_nf,
       f1 = dat_f1,
       f2 = dat_f2,
       f3 = dat_f3)
}

#' @rdname gen_dat_lin
#' @export
gen_dat_nlr <- function(n_ic = 1500, n_oc = 500, phi = .8,
                        a_1 = 3, a_2 = 2) {
  # Create t
  t <- arima.sim(model = list(ar = phi), n = (n_ic + n_oc)) |> scale_t()

  index_ic <- 1:n_ic
  index_oc <- (n_ic+1):(n_ic+n_oc)

  # Generate Training Observations
  dat <-
    tibble::tibble(x1 = exp(t[index_ic]) + rnorm(n_ic, 0, .1),
                   x2 = t[index_ic]^2 - 3*t[index_ic] + rnorm(n_ic, 0, .1),
                   x3 = -t[index_ic]^3 + 3*t[index_ic]^2 + rnorm(n_ic, 0, .1))

  # Generate Testing Observations
  dat_nf <-
    dplyr::bind_rows(dat,
                     tibble::tibble(x1 = exp(t[index_oc]) + rnorm(n_oc, 0, .1),
                                    x2 = t[index_oc]^2 - 3*t[index_oc] + rnorm(n_oc, 0, .1),
                                    x3 = -t[index_oc]^3 + 3*t[index_oc]^2 + rnorm(n_oc, 0, .1))
    )

  dat_f1 <-
    dplyr::bind_rows(dat,
                     tibble::tibble(
                       x1 = exp(t[index_oc]) + rnorm(n_oc, 0, .1),
                       x2 = t[index_oc]^2 - 3*t[index_oc] + rnorm(n_oc, 0, a_1*.1),  # increase variance
                       x3 = -t[index_oc]^3 + 3*t[index_oc]^2 + rnorm(n_oc, 0, .1))
    )

  dat_f2 <-
    dplyr::bind_rows(dat,
                     tibble::tibble(
                       x1 = exp(t[index_oc]) + rnorm(n_oc, 0, .1),
                       x2 = t[index_oc]^2 - 3*t[index_oc] + rnorm(n_oc, 0, .1),
                       x3 = -2*t[index_oc]^3 + 4.5*t[index_oc]^2 + rnorm(n_oc, 0, .1)) # alter x3
    )

  dat_f3 <-
    dplyr::bind_rows(dat,
                     tibble::tibble(
                       x1 = exp(t[index_oc]) + rnorm(n_oc, 0, a_2*.1), # increased variance
                       x2 = t[index_oc]^2 - 3*t[index_oc] + rnorm(n_oc, 0, a_2*.1), # for all
                       x3 = -t[index_oc]^3 + 3*t[index_oc]^2 + rnorm(n_oc, 0, a_2*.1)) # variables
    )

  # Return
  list(none = dat_nf,
       f1 = dat_f1,
       f2 = dat_f2,
       f3 = dat_f3)
}

#' @rdname gen_dat_lin
#' @export
gen_dat_ltm <- function(n_ic = 1500, n_oc = 500, phi = .8, sin_scale = 1.5,
                        a_1 = 3, a_2 = 2) {
  # Create t (heavy linear weight on 22nd lag)
  # t <- arima.sim(model = list(ar = c(.05, rep(0, 20), phi)), n = (n_ic+n_oc)) |> scale_t()

  # Create t (sin wave within bounds of training)
  t <- (arima.sim(model = list(ar = phi), n = (n_ic + n_oc)) +
          rep(sin_scale*sin(seq(0, 2*pi, length.out = 51))[-1], (n_ic+n_oc)/50)) |> scale_t()

  index_ic <- 1:n_ic
  index_oc <- (n_ic+1):(n_ic+n_oc)

  # Generate Training Observations
  dat <-
    tibble::tibble(x1 = exp(t[index_ic]) + rnorm(n_ic, 0, .1),
                   x2 = t[index_ic]^2 - 3*t[index_ic] + rnorm(n_ic, 0, .1),
                   x3 = -t[index_ic]^3 + 3*t[index_ic]^2 + rnorm(n_ic, 0, .1))

  # Generate Testing Observations
  dat_nf <-
    dplyr::bind_rows(dat,
                     tibble::tibble(x1 = exp(t[index_oc]) + rnorm(n_oc, 0, .1),
                                    x2 = t[index_oc]^2 - 3*t[index_oc] + rnorm(n_oc, 0, .1),
                                    x3 = -t[index_oc]^3 + 3*t[index_oc]^2 + rnorm(n_oc, 0, .1))
    )

  dat_f1 <-
    dplyr::bind_rows(dat,
                     tibble::tibble(
                       x1 = exp(t[index_oc]) + rnorm(n_oc, 0, .1),
                       x2 = t[index_oc]^2 - 3*t[index_oc] + rnorm(n_oc, 0, a_1*.1),  # increase variance
                       x3 = -t[index_oc]^3 + 3*t[index_oc]^2 + rnorm(n_oc, 0, .1))
    )
  dat_f2 <-
    dplyr::bind_rows(dat,
                     tibble::tibble(
                       x1 = exp(t[index_oc]) + rnorm(n_oc, 0, .1),
                       x2 = t[index_oc]^2 - 3*t[index_oc] + rnorm(n_oc, 0, .1),
                       x3 = -2*t[index_oc]^3 + 4.5*t[index_oc]^2 + rnorm(n_oc, 0, .1)) # alter x3
    )

  dat_f3 <-
    dplyr::bind_rows(dat,
                     tibble::tibble(
                       x1 = exp(t[index_oc]) + rnorm(n_oc, 0, a_2*.1), # increased variance
                       x2 = t[index_oc]^2 - 3*t[index_oc] + rnorm(n_oc, 0, a_2*.1), # for all
                       x3 = -t[index_oc]^3 + 3*t[index_oc]^2 + rnorm(n_oc, 0, a_2*.1)) # variables
    )

  # Return
  list(none = dat_nf,
       f1 = dat_f1,
       f2 = dat_f2,
       f3 = dat_f3)

}
