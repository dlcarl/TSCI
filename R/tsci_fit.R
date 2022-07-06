#' Two Stage Curvature Identification
#'
#' @description \code{tsci_fit} performs two stage curvature identification for pre-selected treatment model hyperparameters. The function randomly splits the data in two sets \code{A1} and \code{A2}. \code{A2} is used to fit the treatment model given the hyperparameters. \code{A1} is used to fit the outcome model, identify the violation space and to estimate the treatment effect.
#'
#' @param df_treatment a data frame containing the subset of observations of the treatment model.
#' @param Y outcome vector with dimension n by 1.
#' @param D treatment vector with dimension n by 1.
#' @param Z instrument variable with dimension n by 1.
#' @param X baseline covariates with dimension n by p.
#' @param vio_space vio_space a matrix or a list. If a matrix, then each column corresponds to a violation form of Z; If a list, then each element corresponds to a violation form of Z and must be a matrix of n rows, e.g. (Z^3,Z^2); If NULL, then default by the n by 3 matrix (Z^3, Z^2, Z). Violation form selection will be performed according to provided violation forms, for example, null violation space vs Z vs (Z^2, Z) vs (Z^3, Z^2, Z) in the default case.
#' @param intercept logical, including the intercept or not in the outcome model, default by TRUE.
#' @param split_prop split_prop numeric, proportion of observations used to fit the outcome model.
#' @param str_thol minimal value of the threshold of IV strength test, default by 10.
#' @param alpha the significance level, default by 0.05.
#' @param params a list containing the hyperparameters of the treatment model fitting method.
#' @param function_hatmatrix a function to get the hat matrix of the treatment model.
#'
#' @return The output of \code{tsci_selection}
#'
#' @noRd
#'
tsci_fit <- function(df_treatment,
                     Y,
                     D,
                     Z,
                     X,
                     list_vio_space,
                     intercept,
                     str_thol,
                     split_prop,
                     alpha,
                     params,
                     function_hatmatrix) {

  # split data
  n <- NROW(df_treatment)
  n_A1 <- round(split_prop * n)
  A1_ind <- sample(seq_len(n), n_A1)
  df_treatment_A1 <- df_treatment[A1_ind, ]
  df_treatment_A2 <- df_treatment[-A1_ind, ]
  n <- length(Y)
  n_A1 <- length(A1_ind)

  # fit treatment model and get hat matrix.
  model_treatment <- function_hatmatrix(
    df_treatment_A1 = df_treatment_A1,
    df_treatment_A2 = df_treatment_A2,
    params = params
  )

  Y_A1 <- Y[A1_ind, ]
  D_A1 <- D[A1_ind, ]
  if (is.null(X)){
    X_A1 <- NULL
  } else {
    X_A1 <- X[A1_ind, ]
  }

  vio_space <- list_vio_space$vio_space[A1_ind, ]
  rm_ind <- list_vio_space$rm_ind
  Q <- list_vio_space$Q


  # estimate treatment effect on outcome using TSCI.
  outputs <- tsci_selection(
    Y = Y,
    D = D,
    X = X,
    Y_A1 = Y_A1,
    D_A1 = D_A1,
    X_A1 = X_A1,
    vio_space = vio_space,
    rm_ind = rm_ind,
    Q = Q,
    weight = model_treatment$weight,
    intercept = intercept,
    str_thol = str_thol,
    alpha = alpha
  )

  return(outputs)
}
