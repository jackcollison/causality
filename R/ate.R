#' @description Help function to filter out NA values
#' 
#' @param data a dataframe object containing the variables and values.
#' 
#' @return dataframe object that is filtered to drop NAs
#' 
.filter_nas <- function(data){
  #Find NA values
  complete_cases <- complete.cases(data)
  
  #Warning
  if (sum(complete_cases) != nrow(data)) {
    warning(paste0("There are ", nrow(data) - sum(complete_cases), 
                   " rows with missing values. These have been removed"))
  }
  
  #Return
  data[complete_cases, ]
}

#' @description Helper function to return ATE with and without confidence interval
#'
#' @param ATE the average treatment effect as calculated with another function.
#' @param se_hat an estimate of the standard error.
#' @param cf logical; if TRUE, then include confidence interval on ATE.
#' 
#' @return a list of ATE with upper- and lower-bounds or just ATE, depending on user input of \code{cf}.
#' 
.return_helper <- function(ATE, se_hat, cf){
  if (cf == T) {
    #Confidence intervals
    lb <- ATE - 1.96 * se_hat
    ub <- ATE + 1.96 * se_hat
    return(c(ATE = ATE, Lowerbound = lb, Upperbound = ub))
  } else {
    return(ATE)
  }
}

#' Estimate naive average treatment effect (ATE).
#'
#' @param data a dataframe object containing the variables and values.
#' @param y a character vector specifying the response variable.
#' @param w a character vector specifying the treatment status.
#' @param cf logical; if TRUE then includes confidence interval on ATE.
#' 
#' @details Computes a simple difference in means between the treatment group and the control group, 
#' \eqn{\tau = E[Y_i | W = 1] - E[Y_i | W = 0]}. 
#' 
#' @return a list of ATE, 95 percent confidence interval upperbound and lowerbound or just ATE, depending on user input of \code{cf}
#' 
#' @export
#' 
#' @examples
#' data("lalonde")
#' ate <- naive_ate(data = lalonde, y = "re78", w = "treat")
#' 
naive_ate <- function(data, y, w, cf = TRUE){
  #Filter for NAs
  data <- .filter_nas(data)
  
  #Split into control and treatment
  control <- data[data[w] == 0, y]
  treatment <- data[data[w] == 1, y]
  
  #Sample sizes
  c_size <- length(control)
  t_size <- length(treatment)
  
  #Average treatment effect
  ATE <- mean(treatment) - mean(control)
  se_hat <- sqrt(var(control)/(c_size-1) + var(treatment)/(t_size-1))
  
  .return_helper(ATE, se_hat, cf)
}


#' Etimate average treatment effect (ATE) with
#' inverse propensity score weighting.
#'
#' @param data a dataframe object containing the variables and values.
#' @param y a character vector specifying the response variable.
#' @param w a character vector specifying the treatment status.
#' @param p a vector containing propensity score values.
#' @param cf logical; if TRUE then includes confidence interval on ATE.
#' 
#' @details Computes an estimate of the ATE \eqn{\tau} using propensity score weighting: 
#' \deqn{\tau = E \left[ \frac{Y_i W_i}{e(X_i)} - \frac{Y_i (1 - W_i)}{1 - e(X_i)} \right]}
#' where \eqn{e(X_i) = P(W_i = 1 | X_i = x)} is the propensity score. 
#' 
#' @return a list of ATE, 95 percent confidence interval upperbound and lowerbound or just ATE, depending on user input of \code{cf}
#' 
#' @references Robins, James M., Andrea Rotnitzky and Lue Ping Zhao. 1994. "Estimation of Regression Coefficients 
#' When Some Regressors Are Not Always Observed." \emph{Journal of the American Statistical Association}. 
#' Vol. 89(427), Sep., pgs. 846-866. \url{https://doi.org/10.1080/01621459.1994.10476818}

#' Lunceford JK, Davidian M. 2004. ``Stratification and weighting via the propensity score in estimation of causal 
#' treatment effects: a comparative study." \emph{Statistics in Medicine}. Vol. 23(19), Aug., pgs. 2937–2960. 
#' \url{https://doi.org/10.1002/sim.1903}
#' 
#' @export
#' 
#' @examples
#' data("lalonde")
#' 
#' # calculate propensity scores
#' p <- propensity_score(data = lalonde, y = "re78", w = "treat", model = "logistic")
#' ate <- ipw_ate(data = lalonde, y = "re78", w = "treat", p = p)
#' 
ipw_ate <- function(data, y, w, p, cf = T) {
  #Filter for NAs
  data <- .filter_nas(data)
  
  #Create variables for brevity
  W <- data[[w]]
  Y <- data[[y]]
  
  if (length(W) != length(p)) {
    stop("Data and propensity scores do not have the same length")
  }
  
  #Create weighting
  G <- ((W - p) * Y) / (p * (1 - p))
  
  #Estimate average treatment effect and standard error
  ATE <- mean(G)
  se_hat <- sqrt(var(G) / (length(G) - 1))
  
  .return_helper(ATE, se_hat, cf)
}


#' Estimate average treatment effect (ATE) with
#' inverse propensity score weighting using OLS.
#'
#' @param data a dataframe object containing the variables and values.
#' @param y a character vector specifying the response variable.
#' @param w a character vector specifying the treatment status.
#' @param p a vector containing propensity score values.
#' @param cf logical; if TRUE then includes confidence interval on ATE.
#' 
#' @details Computes an estimate of the ATE \eqn{\tau} via weighted OLS regression using propensity score weighting
#' The weights are given by:
#' \deqn{w_i = \frac{W_i}{e(X_i)} - \frac{1 - W_i}{1 - e(X_i)}}
#' where \eqn{e(X_i) = P(W_i = 1 | X_i = x)} is the propensity score. 
#' 
#' @return a list of ATE, 95 percent confidence interval upperbound and lowerbound or just ATE, depending on user input of \code{cf}
#' 
#' @references  Aronow, Peter M.; Samii, Cyrus. ``Estimating average causal effects under general interference, with application to a social network experiment." 
#' \emph{Annals of Applied Statistics} 11 (2017), no. 4, 1912--1947
#' \url{https://arxiv.org/pdf/1305.6156.pdf}
#' 
#' @export
#' 
#' @examples 
#' data("lalonde")
#' 
#' # calculate propensity scores
#' p <- propensity_score(lalonde, y = "re78", w = "treat", model = "logistic")
#' ate <- prop_weighted_ols_ate(data = lalonde, y = "re78", w = "treat", p = p)
#' 
prop_weighted_ols_ate <- function(data, y, w, p, cf = TRUE) {
  #Filter for NAs
  data <- .filter_nas(data)
  
  #Create variables for brevity
  W <- data[[w]]
  Y <- data[[y]]
  
  if (length(W) != length(p)) {
    stop("Data and propensity scores do not have the same length")
  }
  
  weights <- (W / p) + ((1 - W) / (1 - p))
  
  #Initial weighted linear fit
  linear_fit <- lm(Y ~ W, weights = weights)
  
  #Estimate average treatment effect and standard error
  ATE <- as.numeric(coef(linear_fit)["W"])
  se_hat <- as.numeric(sqrt(sandwich::vcovHC(linear_fit)["W", "W"]))
  
  .return_helper(ATE, se_hat, cf)
}

#' Estimate average treatment effect (ATE) with
#' augmented inverse propensity score weighting.
#'
#' @param data a dataframe object containing the variables and values.
#' @param x a list of character vectors specifying variables to be included in 
#' the model (columns in the data). If unspecified, then it is assumed to be all 
#' columns in the data besides y and w.
#' @param y a character vector specifying the response variable.
#' @param w a character vector specifying the treatment status.
#' @param p a vector containing propensity score values.
#' @param cf logical; if TRUE then includes confidence interval on ATE.
#' 
#' @details Estimates the ATE \eqn{\tau} using the doubly robust method described in Scharfstein, Robins and Rotznitzky (1998) that combines
#' both regression and propensity score weighting. 
#' \deqn{\tau = E \left[  W_i \frac{Y_i-\tau(1,X_i)}{e(X_i)} + (1-W_i) \frac{Y_i-\tau(0,X_i)}{1-e(X_i)} + \tau(1,X_i) - \tau(0,X_i)\right]}
#' where \eqn{e(X_i) = P(W_i = 1 | X_i = x)} is the propensity score and \eqn{\tau(1, X_i)} and \eqn{\tau{0, X_i}} are estimated in
#' the first stage via OLS regression. 
#' 
#' @return a list of ATE, 95 percent confidence interval upperbound and lowerbound or just ATE, depending on user input of \code{cf}.
#' 
#' @references Rotnitzky, Andrea, James M. Robins, and Daniel O. Scharfstein. 1998. ``Semiparametric 
#' Regression for Repeated Outcomes with Nonignorable Nonresponse." \emph{Journal of the American Statistical Association}. 
#' Vol. 93, No. 444, Dec. pgs. 1321-1339.
#' 
#' Cao, Weihua, Anastasios A. Tsiatis, and Marie Davidian. 2009. ``Improving efficiency and robustness of the doubly 
#' robust estimator for a population mean with incomplete data". \emph{Biometrika}. Vol. 96(3). pgs. 723–734.
#' 
#' @export
#' 
#' @examples 
#' data("lalonde")
#' 
#' # calculate propensity scores
#' p <- propensity_score(lalonde, y = "re78", w = "treat", model = "logistic")
#' ate <- aipw_ate(data = lalonde, y = "re78", w = "treat", p = p)
#' 
aipw_ate <- function(data, x, y, w, p, cf = TRUE){
  #Filter for NAs
  data <- .filter_nas(data)
  
  if (nrow(data) != length(p)) {
    stop("Data and propensity scores do not have the same length")
  }
  
  #Check if missing
  if (missing(x)) x <- setdiff(names(data), c(w, y))
  
  data <- data[c(x, y, w)]
  
  formula <- as.formula(paste0(y, " ~ ", w, " * ."))
  ols.fit <-  lm(formula, data = data)
  
  data.treatall <- data
  data.treatall[, w] <-  rep(1, nrow(data))
  treated_pred <-  predict(ols.fit, data.treatall)
  
  data.controlall <-  data
  data.controlall[, w] <- rep(0, nrow(data))
  control_pred <- predict(ols.fit, data.controlall)
  
  actual_pred <-  predict(ols.fit, data)
  
  G <- treated_pred - control_pred +
    ((data[, w] - p) * (data[, y] - actual_pred)) / (p * (1 - p))
  
  ATE <- mean(G)
  se_hat <- sqrt(var(G) / (length(G) - 1))
  
  .return_helper(ATE, se_hat, cf)
}


#' Estimate average treatment effect (ATE) with
#' propensity score stratification.
#'
#' @param data a dataframe object containing the variables and values.
#' @param y a character vector specifying the response variable.
#' @param w a character vector specifying the treatment status.
#' @param p a vector containing propensity score values.
#' @param n_strata the number of strata to use when splitting propensity scores. Default: 4
#' @param min_obs Minimum number of observations per strata. 
#' If this is specified as a valid nontrivial number (greater than 1, less than non NA data size), 
#' then n_strata will be set to the minimum of the specified n_strata, and the n_strata such that min_obs is satisfied. 
#' @param model a model to estimate ATE (one of "naive", "ipw", "weighted ols", or "aipw")
#' 
#' @details Estimates average treatment effects across different groups of propensity scores specified
#' by the user (either with a minimum number of observations in a strata or number of strata). The user
#' can also specify the method of estimation.
#' 
#' The naive estimator is given by:
#' 
#' \eqn{\tau = E[Y_i | W = 1] - E[Y_i | W = 0]}.
#' 
#' The inverse propensity score weighting estimator is given by:
#' 
#' \deqn{\tau = E \left[ \frac{Y_i W_i}{e(X_i)} - \frac{Y_i (1 - W_i)}{1 - e(X_i)} \right]}
#' where \eqn{e(X_i) = P(W_i = 1 | X_i = x)} is the propensity score.
#' 
#' The weighted OLS model is a regression with weights given by:
#' 
#' \deqn{w_i = \frac{W_i}{e(X_i)} - \frac{1 - W_i}{1 - e(X_i)}}
#' where \eqn{e(X_i) = P(W_i = 1 | X_i = x)} is the propensity score. 
#' 
#' Finally, the augmented inverse propensity score weighted method is given by:
#' 
#' \deqn{\tau = E \left[  W_i \frac{Y_i-\tau(1,X_i)}{e(X_i)} + (1-W_i) \frac{Y_i-\tau(0,X_i)}{1-e(X_i)} + \tau(1,X_i) - \tau(0,X_i)\right]}
#' where \eqn{e(X_i) = P(W_i = 1 | X_i = x)} is the propensity score and \eqn{\tau(1, X_i)} and \eqn{\tau{0, X_i}} are estimated in
#' the first stage via OLS regression.
#' 
#' More details on each of these functions can be found in their individual estimators. The ATE is given by
#' the average of the estimates over the strata:
#' 
#' \eqn{\frac{1}{n}\sum_{i=1}^{n}\tau_i}
#' where \eqn{\tau_i} is the average treatment effect estimated in strata \eqn{i} and \eqn{n} is the number
#' of strata specified by the user (or calculated with the minimum number of observations per strata specified).
#' 
#' @return numeric estimate of the average treatment effect.
#' 
#' @references Rosenbaum, Paul, and Donald Rubin. 1984. ``Reducing Bias in Observational Studies Using Subclassification 
#' on the Propensity Score." \emph{Journal of the American Statistical Association}. Vol. 79(387), Feb. pgs 516-524.
#' 
#' Lunceford JK, Davidian M. 2004. ``Stratification and weighting via the propensity score in estimation of causal 
#' treatment effects: a comparative study." \emph{Statistics in Medicine}. Vol. 23(19), Aug., pgs. 2937–2960. 
#' \url{https://doi.org/10.1002/sim.1903}
#' 
#' @export
#' 
#' @examples 
#' data("lalonde")
#' 
#' # calculate propensity scores
#' p <- propensity_score(lalonde, y = "re78", w = "treat", model = "logistic")
#' ate <- propensity_strat_ate(data = lalonde, 
#'                             y = "re78", 
#'                             w = "treat",
#'                             p = p, 
#'                             n_strata = 4, 
#'                             min_obs = 1)
#' 
propensity_strat_ate <- function(data, y, w, p, n_strata = 4, min_obs = 1, model = "naive"){
  
  #Filter for NAs
  data <- .filter_nas(data)
  
  if (nrow(data) != length(p)) {
    stop("Data and propensity scores do not have the same length")
  }
  
  if (n_strata %% 1 != 0 | min_obs %% 1 != 0 | (min_obs < 0) | (n_strata <= 0)) {
    stop("n_strata and min_obs must be valid integers")
  }
  
  if (nrow(data) < n_strata) {
    stop("Invalid n_strata, n_strata greater than number of rows")
  }
  
  if (min_obs > nrow(data)) {
    stop("Minimum number of observations exceeds dataset size")
  }
  
  #min_obs is specified to a positive integer, then manually set n_strata if the specified n strata violates min_obs
  if (min_obs > 1) {
    n_strata_limit = floor(nrow(data) / min_obs)
    n_strata = min(n_strata, n_strata_limit)
  }
  #Instantiate aggregate ATE
  ate_aggregate <- 0
  
  #Split into quantiles
  quantile <- as.numeric(cut(p, breaks = quantile(p, probs = seq(0, 1, length = n_strata), na.rm=TRUE),
                             include.lowest=TRUE))
  
  #Loop over quantiles and estimate ATE for each stratum
  for (i in seq(1:n_strata)) {
    indicies <- which(quantile == i)
    stratum <- data[indicies, ]
    
    #Check if there exists data or if all of one class
    if (nrow(stratum) == 0 | length(unique(stratum[[w]])) < 2) {
      next
    }
    
    #Calculate ATE and aggregate
    if (model == "naive") {
      ate_current <- naive_ate(stratum, y, w, cf = F)
    } else if (model == "ipw") {
      ate_current <- ipw_ate(data, y, w, p, cf = F)
    } else if (model == "weighted ols") {
      ate_current <- prop_weighted_ols_ate(data, y, w, p, cf = F)
    } else if (model == "aipw") {
      ate_current <- aipw_ate(data = data, y = y, w = w, p = p, cf = F)
    } else {
      errorCondition("User did not input valid ATE estimator.")
    }
    ate_aggregate <- ate_aggregate + ate_current
  }
  
  return(ate_aggregate / n_strata)
}

#' Estimate average treatment effect (ATE) with
#' double selection methodology.
#' 
#' @param data a dataframe object containing the variables and values.
#' @param x a list of character vectors specifying variables to be included in 
#' the model (columns in the data). If unspecified, then it is assumed to be all 
#' columns in the data besides y and w.
#' @param y a character vector specifying the response variable.
#' @param w a character vector specifying the treatment status.
#' @param nfolds integer describing number of folds to use in k-fold cross-validation; preset to 5
#' 
#' @details Estimates the average treatment effect \eqn{\tau} using the methodology developed in
#' Belloni, Chernozhukov, and Hansen (2014), which they term the ``post-double-selection" method. 
#' The general procedure performed here is as follows: 
#' \enumerate{
#'   \item Predict the treatment \eqn{W_i} using the covariates \eqn{X_i} using lasso regression (where \eqn{\lambda} is tuned 
#'   using cross-validation). 
#'   Select the covariates that have non-zero coefficients in the lasso model. 
#'   \item Predict the outcome \eqn{Y_i} using the covariates \eqn{X_i} using lasso regression (where \eqn{\lambda} is tuned 
#'   using cross-validation). 
#'   Select the covariates that have non-zero coefficients in the lasso model. 
#'   \item Estimate the treatment effect \eqn{\tau} by the linear regression of \eqn{Y_i} on the treatment \eqn{W_i}
#'   and the union of the set of variables selected in the two covariate selection steps.
#' }
#' 
#' @return numeric estimate of the average treatment effect.
#' 
#' @references Belloni, Alexandre, Victor Chernozhukov, and Christian Hansen. 2014. ``High-Dimensional Methods 
#' and Inference on Structural and Treatment Effects." \emph{Journal of Economic Perspective}.
#' Vol. 28, Num. 2, Spring 2014. pgs. 29–50.
#' \url{https://www.aeaweb.org/articles?id=10.1257/jep.28.2.29}
#' 
#' @export
#' 
#' @examples 
#' data("lalonde")
#' ate <- double_selection_ate(data = lalonde, y = "re78", w = "treat")
#' 
double_selection_ate <- function(data, x, y, w, nfolds = 5){
  # 1. predict w with x values using lasso, lambda cv selection (or specify lambdas)
  # 2. predict y with x values using lasso on treated and control, lambda cv (or specify lambda)
  # 3. run ols on y with union of x variables and w
  if (missing(x)) x <- setdiff(names(data), c(w, y))
  
  #Filter for NAs
  data <- .filter_nas(data)
  
  #First step: predict w with x using lasso (CV)
  formula <- as.formula(paste0(w, " ~ ", paste(x, collapse = ' + ')))
  X <- scale(model.matrix(formula, data = data))[, -1]
  step1 <- glmnet::cv.glmnet(X, data[[w]], alpha = 1, nfolds = nfolds)
  coefs1 <- coef(step1, s = "lambda.min")
  selected_covariates1 <- coefs1@Dimnames[[1]][(coefs1@i + 1)[-1]]
  
  #Second step: predict y with x using lasso on treated and control
  formula <- as.formula(paste0(y, " ~ ", paste(x, collapse = ' + ')))
  X <- scale(model.matrix(formula, data = data))[, -1]
  step2 <- glmnet::cv.glmnet(X, data[[y]], alpha = 1, nfolds = nfolds)
  coefs2 <- coef(step2, s = "lambda.min")
  selected_covariates2 <- coefs2@Dimnames[[1]][(coefs2@i + 1)[-1]]
  
  #Third step: run OLS on y with union of x and w variables
  covariates <- unique(c(selected_covariates1, selected_covariates2))
  formula <- as.formula(paste0(y, " ~ ", paste(w), "+", paste(covariates, collapse = ' + ')))
  fit <- lm(formula, data)
  
  #Return ATE
  ATE <- coef(fit)[w]
  ATE
}



