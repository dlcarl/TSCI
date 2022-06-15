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
#' @return
#'     \item{\code{Coef_all}}{a series of point estimators of treatment effect corresponding to different violation spaces and the OLS}
#'     \item{\code{sd_all}}{standard errors of Coef_all}
#'     \item{\code{CI_all}}{confidence intervals for the treatment effect corresponding to different violation spaces and the OLS}
#'     \item{\code{Coef_robust}}{the point estimators corresponding to the violation space selected by the robust comparison}
#'     \item{\code{sd_robust}}{the standard errors of Coef_robust}
#'     \item{\code{CI_robust}}{confidence intervals for the treatment effect with the violation space selected by the robust comparison}
#'     \item{\code{iv_str}}{IV strength corresponding to different violation spaces}
#'     \item{\code{iv_thol}}{the threshold of IV strength test corresponding to different violation spaces}
#'     \item{\code{Qmax}}{the index of largest violation space selected by IV strength test. If -1, the IV strength test fails for null violation space and run OLS. If 0, the IV Strength test fails for the first violation space and run TSRF only for null violation space. In other cases, violation space selection is performed}
#'     \item{\code{q_hat}}{the index of estimated violation space corresponding to Qmax}
#'     \item{\code{invalidity}}{invalidity of TSLS. If TRUE, the IV is invalid; Otherwise, the IV is valid}
#' @noRd
#'
tsci_fit <- function(df_treatment,
                     Y,
                     D,
                     Z,
                     X,
                     vio_space,
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

  # fit treatment model and get hat matrix.
  model_treatment <- function_hatmatrix(
    df_treatment_A1 = df_treatment_A1,
    df_treatment_A2 = df_treatment_A2,
    params = params
  )


  # estimate treatment effect on outcome using TSCI.
  outputs <- tsci_secondstage(
    Y = Y,
    D = D,
    Z = Z,
    X = X,
    vio_space = vio_space,
    A1_ind = A1_ind,
    weight = model_treatment$weight,
    intercept = intercept,
    str_thol = str_thol,
    alpha = alpha
  )

  return(outputs)
}
