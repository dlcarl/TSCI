#' TSCI with single data splitting.
#'
#' @param df_treatment a data frame containing the subset of observations of the treatment model.
#' @param Y outcome vector with dimension n by 1.
#' @param D treatment vector with dimension n by 1.
#' @param Z instrument variable with dimension n by 1.
#' @param X baseline covariates with dimension n by p.
#' @param vio_space vio_space a matrix or a list. If a matrix, then each column corresponds to a violation form of Z; If a list, then each element corresponds to a violation form of Z and must be a matrix of n rows, e.g. (Z^3,Z^2); If NULL, then default by the n by 3 matrix (Z^3, Z^2, Z). Violation form selection will be performed according to provided violation forms, for example, null violation space vs Z vs (Z^2, Z) vs (Z^3, Z^2, Z) in the default case.
#' @param intercept logical, including the intercept or not in the outcome model, default by TRUE.
#' @param A1_ind the indices of samples in A1, used for constructing the point estimator and the confidence interval.
#' @param str_thol minimal value of the threshold of IV strength test, default by 10.
#' @param alpha the significance level, default by 0.05.
#' @param params a list containing the hyperparameters of the treatment model fitting method.
#' @param function_hatmatrix a function to get the hat matrix of the treatment model.
#' @param save_model logcial, if \code{TRUE} the fitted boosting model is returned.
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
#'     \item{\code{treatment_model}}{the fitted boosting model. Will be returned only if save_model is \code{TRUE}}
#' @noRd
#'
single_split <- function(df_treatment,
                         Y,
                         D,
                         Z,
                         X,
                         vio_space,
                         A1_ind,
                         intercept,
                         str_thol,
                         alpha,
                         params,
                         function_hatmatrix,
                         save_model) {
  # split sample
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
  if (save_model == TRUE) {
    outputs <- append(outputs, list("treatment_model" = model_treatment$model))
  }
  return(outputs)
}
