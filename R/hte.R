#' @description Helper function to filter out NA values
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

#' @description Helper function to fit regression forests, OLS
#' 
#' @param model a string that describes the model being fit.
#' @param X a matrix of independent variables
#' @param Y a matrix of dependent variables
#' @param ... additional arguments to regression_forest
#' 
#' @return an lm or random forest object
#' 
.fit_model <- function(model, X, Y, ...){
  if (model == "regression forest") {
    fit <- grf::regression_forest(X = X, Y = Y, ...)
  } else if (model == "OLS") {
    data <- cbind(X, y = Y)
    fit <- lm(y ~ ., data = data, ...)
  } else {
    stop("must specify 'regression forest' or 'OLS' for `model`.")
  }
  fit
}

#' @description Helper function to predict with regression forests, OLS
#' 
#' @param model a string that describes the model being fit.
#' @param ... additional arguments to predict
#' 
#' @return vector of predictions
#' 
.predict <- function(model, ...){
  if (model == "regression forest") {
    preds <- predict(...)$predictions
  } else if (model == "OLS") {
    preds <- predict(...)
  }
  preds
}

#' @description Helper function to coerce nonnumeric columns
#' 
#' @param data a dataframe object containing the variables and values.
#' 
#' @return a dataframe with coerced nonnumeric columns
#' 
.coerce_nonnumeric_cols <- function(data){
  nonnumeric_to_numeric <- function(data, x){
    if(!is.numeric(data[[x]])) {
      data[[x]] <- as.numeric(as.factor(data[[x]]))  
      warning(paste0("`regression_forest` does not currently support non-numeric features. ",
      "Coerced column `", x, "` to numeric. If one-hot encoding is preferred, do so prior to calling this function."))
    }
    data[[x]]
  }
  data <- sapply(colnames(data), function(x) nonnumeric_to_numeric(data, x) )
  data
}

#' Estimate heterogeneous treatment effects (HTEs)
#' using the S-Learner strategy.
#'
#' @param data a dataframe object containing the variables and values.
#' @param x a list of character vectors specifying variables to be included in 
#' the model (columns in the data). If unspecified, then it is assumed to be all 
#' columns in the data besides y and w.
#' @param y a character vector specifying the response variable.
#' @param w a character vector specifying the treatment status.
#' @param base_learner a character vector specifying the base learner. One of "regression forest"
#' or "OLS". Default is "regression forest".
#' @param plot logical; if TRUE, then plots a histogram of treatment effects.
#' @param ... additional arguments to the base learner.
#' 
#' @details Implements the S-learner algorithm described in Künzel et al. (2019) for estimating
#' conditional average treatment effects (CATE). In the S-learner algorithm, the treatment \eqn{W}
#' is included as a feature similar to all of the other covariates without the indicator being given any special role.
#' The combined response function 
#' \deqn{\mu(x, w) = E \left[ Y^{obs} | X = x, W = w \right]}
#' can then be estimated using any base learner (supervised machine learning or regression algorithm) on the entire dataset. 
#' Here we implement the S-learner with the option forA linear regression or a regression forest (see Athey, Tibshirani, and Wager (2016)) 
#' as the base learner.
#' 
#' The CATE estimator is then given by 
#' \deqn{\hat{\tau}(x) = \hat{\mu}(x, 1) - \hat{\mu}(x, 0).}
#' 
#' 
#' @return a list of two. The first element is a vector of conditional average treatment effect for each observation. 
#' The second element is the estimated average treatment effect. 
#' 
#' @references Künzel, Sören R., Jasjeet S. Sekhon, Peter J. Bickel, and Bin Yu. 2019. 
#' ``Metalearners for estimating heterogeneous treatment effects using machine learning." 
#' \emph{Proceedings of the National Academy of Sciences of the United States of America.}
#' Mar. 116(10): 4156–4165. \url{https://doi.org/10.1073/pnas.1804597116}
#' 
#' Athey, Susan, Julie Tibshirani, and Stefan Wager. 2016. ``Generalized Random Forests."
#' Working paper; Forthcoming in the Annals of Statistics. \url{https://arxiv.org/abs/1610.01271}
#' 
#' @seealso \code{\link{t_learner}}, \code{\link{x_learner}}
#' 
#' @export
#' 
#' @examples 
#' data("lalonde")
#' hte <- s_learner(data = lalonde, y = "re78", w = "treat", num.trees = 100, mtry = 3)
#' 
s_learner <- function(data, x, y, w, base_learner = "regression forest", plot = TRUE, ...) {
  #Filter for NAs
  data <- .filter_nas(data)
  
  if (base_learner != "regression forest" && base_learner != "OLS"){
    stop("Base learner invalid")
  }
  
  #Check if missing
  if (missing(x)) x <- setdiff(names(data), c(w, y))
  
  if (base_learner == "regression forest") {
    data[names(data) != w] <- .coerce_nonnumeric_cols(data[names(data) != w])
  }
  
  if(base_learner == "regression forest") 
    data[names(data) != w] <- .coerce_nonnumeric_cols(data[names(data) != w])
  
  #Separate data for brevity
  X <- data[c(x, w)]
  Y <- data[[y]]
  
  #Fit regression forest
  s_fit <- .fit_model(model = base_learner, X = X, Y = Y, ...)
  
  #Create treatment and control
  all_treat <- X
  all_treat[w] <- 1
  all_control <- X
  all_control[w] <- 0
  
  #Estimate HTE using difference in predict functions
  y_treat <- .predict(base_learner, s_fit, newdata = all_treat)
  y_control <- .predict(base_learner, s_fit, newdata = all_control)
  
  tauhat <- y_treat - y_control
  
  #Plot histogram of estimate
  if (plot == T) {
    hist(tauhat, main = "CATE Estimated by S-Learner", xlab = "CATE", ylab = "Frequency")
  }
  
  list(CATE = tauhat, ATE = mean(tauhat))
}

#' Estimate heterogeneous treatment effects (HTEs)
#' using the T-Learner strategy.
#'
#' @param data a dataframe object containing the variables and values.
#' @param x a list of character vectors specifying variables to be included in 
#' the model (columns in the data). If unspecified, then it is assumed to be all 
#' columns in the data besides y and w.
#' @param y a character vector specifying the response variable.
#' @param w a character vector specifying the treatment status.
#' @param base_learner a character vector specifying the base learner. One of "regression forest"
#' or "OLS". Default is "regression forest".
#' @param plot logical; if TRUE, then plots a histogram of treatment effects.
#' @param ... additional arguments passed to the base learner.
#' 
#' @details Implements the T-learner algorithm described in Künzel et al. (2019) for estimating
#' conditional average treatment effects (CATE). In the T-learner algorithm, the control response
#' function is estimated using all units in the control group as
#' \deqn{\mu_0 = E [ Y(0) | X = x],}
#' and the treatment response function is estimates using all units in the treatment group as
#' \deqn{\mu_1 = E [ Y(1) | X = x].}
#' Both \eqn{\mu_0 } and \eqn{\mu_1 } are estimated using any base learner (supervised machine learning or 
#' regression algorithm). Here we implement the T-learner with the option for linear regression or regression forest 
#' (see Athey, Tibshirani, and Wager (2016)) as the base learner.
#' 
#' The CATE is then estimated in the second stage as 
#' \deqn{\hat{\tau}(x) =  \hat{\mu}(x, 1) - \hat{\mu}(x, 0).}
#' 
#' 
#' @return a list of two. The first element is a vector of conditional average treatment effect for each observation. 
#' The second element is the estimated average treatment effect. 
#' 
#' @references Künzel, Sören R., Jasjeet S. Sekhon, Peter J. Bickel, and Bin Yu. 2019. 
#' ``Metalearners for estimating heterogeneous treatment effects using machine learning." 
#' \emph{Proceedings of the National Academy of Sciences of the United States of America.}
#' Mar. 116(10): 4156–4165. \url{https://doi.org/10.1073/pnas.1804597116}
#' 
#' Athey, Susan, Julie Tibshirani, and Stefan Wager. 2016. ``Generalized Random Forests."
#' Working paper; Forthcoming in the Annals of Statistics. \url{https://arxiv.org/abs/1610.01271}
#' 
#' @seealso \code{\link{s_learner}}, \code{\link{x_learner}}
#' 
#' @export
#' 
#' @examples 
#' data("lalonde")
#' hte <- t_learner(data = lalonde, y = "re78", w = "treat", num.trees = 100, mtry = 3)
#' 
t_learner <- function(data, x, y, w, base_learner = "regression forest", plot = TRUE, ...) { 
  #Filter for NAs
  data <- .filter_nas(data)
  
  #Check if missing
  if (missing(x)) x <- setdiff(names(data), c(w, y))
  
  if (base_learner != "regression forest" && base_learner != "OLS"){
    stop("Base learner invalid")
  }
  
  if (base_learner == "regression forest") {
    data[names(data) != w] <- .coerce_nonnumeric_cols(data[names(data) != w])
  }
  
  if(base_learner == "regression forest") 
    data[names(data) != w] <- .coerce_nonnumeric_cols(data[names(data) != w])
  
  #Separate treatment for brevity
  treatment <- data[data[,w]==1,] 
  treatment[w] <- NULL
  treat_X <- treatment[x]
  treat_Y <- treatment[[y]]
  
  #Separate control for brevity
  control <- data[data[,w]==0,] 
  control[w] <- NULL
  control_X <- control[x]
  control_Y <- control[[y]]
  
  #Fit regression forest on treatment data
  treatment_fit <- .fit_model(model = base_learner, X = treat_X, Y = treat_Y, ...) 
  
  #Fit regression forest on control data
  control_fit <- .fit_model(model = base_learner, X = control_X, Y = control_Y, ...)
  
  #Predict with results from regression forests
  data[w] <- NULL
  data[y] <- NULL
  y_treat <- .predict(base_learner, treatment_fit, data[x])
  y_control <- .predict(base_learner, control_fit, data[x])
  
  #Estimate HTE using difference in predict functions
  tauhat <- y_treat - y_control
  
  #Plot histogram
  if (plot == T) {
    hist(tauhat, main = "CATE Estimated by T-Learner", xlab = "CATE", ylab = "Frequency")
  }
  
  list(CATE = tauhat, ATE = mean(tauhat))
}

#' Estimate heterogeneous treatment effects (HTEs)
#' using the X-Learner strategy.
#'
#' @param data a dataframe object containing the variables and values.
#' @param x a list of character vectors specifying variables to be included in 
#' the model (columns in the data). If unspecified, then it is assumed to be all 
#' columns in the data besides y and w.
#' @param y a character vector specifying the response variable.
#' @param w a character vector specifying the treatment status.
#' @param base_learner a character vector specifying the base learner. One of "regression forest"
#' or "OLS". Default is "regression forest".
#' @param cate_model a character vector specifying the model used to estimate the CATE. One of "regression forest"
#' or "OLS". Default is "regression forest".
#' @param propensity_model  character string naming the model used to estimate propensity scores (one of "logistic", 
#' "lasso", or "causal forest").
#' @param plot logical; if TRUE, then plots a histogram of treatment effects.
#' @param ... additional arguments passed to the base learner.
#' 
#' @details Implements the X-learner algorithm proposed in Künzel et al. (2019) for estimating
#' conditional average treatment effects (CATE). 
#' 
#' In the first stage of the X-learner, the control response
#' function is estimated using all units in the control group as
#' \deqn{\mu_0 = E [ Y(0) | X = x],}
#' and the treatment response function is estimates using all units in the treatment group as
#' \deqn{\mu_1 = E [ Y(1) | X = x].}
#' Both \eqn{\mu_0 } and \eqn{\mu_1 } are estimated using any base learner (supervised machine learning or 
#' regression algorithm). Here we implement the X-learner with linear regression or regression forest 
#' (see Athey, Tibshirani, and Wager (2016)) as the base learner. 
#' 
#' In the second stage, the treatment effect for each observation is then imputed by estimating the counterfactual
#' outcome for each observation using the first-stage base learner models:
#' \deqn{\tilde{D}^1_i := Y^1_i - \hat{\mu}_0(X^1_i)}
#' and 
#' \deqn{\tilde{D}^0_i := Y^0_i - \hat{\mu}_0(X^0_i)}
#' where \eqn{\tilde{D}^1_i} and \eqn{\tilde{D}^1_i} are the imputed treatment effects (two for each observation).
#' The CATE is then estimated in two ways: 
#' \deqn{\hat{\tau}_1 = E[\tilde{D}^1 | X = x]}
#' and 
#' \deqn{\hat{\tau}_0 = E[\tilde{D}^0 | X = x].}
#' Currently, we include the option to estimate \eqn{\hat{\tau}_1} and \eqn{\hat{\tau}_0} with linear regression 
#' or regression forests.
#' 
#' In the third stage, estimate the CATE by a weighted average of the two estimates from the second stage:
#' \deqn{\hat{\tau} = g(x) \hat{\tau}_0(x) + (1 - g(x)) \hat{\tau}_1(x).}
#' 
#' Here, we choose propensity scores to be the weighting function \eqn{g{x}}.
#' 
#' 
#' @return a list of two. The first element is a vector of conditional average treatment effect for each observation. 
#' The second element is the estimated average treatment effect. 
#' 
#' @references Künzel, Sören R., Jasjeet S. Sekhon, Peter J. Bickel, and Bin Yu. 2019. 
#' ``Metalearners for estimating heterogeneous treatment effects using machine learning." 
#' \emph{Proceedings of the National Academy of Sciences of the United States of America.}
#' Mar. 116(10): 4156–4165. \url{https://doi.org/10.1073/pnas.1804597116}
#' 
#' Athey, Susan, Julie Tibshirani, and Stefan Wager. 2016. ``Generalized Random Forests."
#' Working paper; Forthcoming in the Annals of Statistics. \url{https://arxiv.org/abs/1610.01271}
#' 
#' @seealso \code{\link{s_learner}}, \code{\link{t_learner}}
#' 
#' @export
#' 
#' @examples 
#' data("lalonde")
#' hte <- x_learner(data = lalonde, y = "re78", w = "treat", num.trees = 100, mtry = 3)
#' 
x_learner <- function(data, 
                      x, 
                      y,
                      w, 
                      base_learner = "regression forest", 
                      cate_model = "regression forest", 
                      propensity_model = "logistic", 
                      plot = TRUE, 
                      ...){

  #Filter for NAs
  data <- .filter_nas(data)
  
  if (base_learner != "regression forest" && base_learner != "OLS"){
    stop("Base learner invalid")
  }
  
  if (cate_model != "regression forest" && cate_model != "OLS"){
    stop("Base learner invalid")
  }
  
  if (propensity_model != "causal forest" && propensity_model != "lasso" && propensity_model != "logistic"){
    stop("Propensity Model invalid")
  }
  
  if (base_learner == "regression forest") {
    data[names(data) != w] <- .coerce_nonnumeric_cols(data[names(data) != w])
  }
  
  #Separate data for brevity
  X <- data[x]
  Y <- data[[y]]
  W <- data[[w]]
  
  #Separate into treatment and control
  X_control <- X[W == 0, ]
  Y_control <- Y[W == 0]
  X_treat <- X[W == 1, ]
  Y_treat <- Y[W == 1]
  
  #Fit regression forest on control and predict treatment
  control_fit <- .fit_model(model = base_learner, X = X_control, Y = Y_control, ...)
  yhat_control <- .predict(base_learner, control_fit, X_treat)
  delta_treat <- Y_treat - yhat_control
  cate_fit1 <-  .fit_model(model = cate_model, X = X_treat, Y = delta_treat, ...)
  tauhat_treat <- .predict(cate_model, cate_fit1, X)
  
  #Fit regression forest on treatment and predict control
  treat_fit <- .fit_model(model = base_learner, X = X_treat, Y = Y_treat)
  yhat_treat <- .predict(base_learner, treat_fit, X_control)
  delta_control <- yhat_treat - Y_control
  cate_fit0 <- .fit_model(model = cate_model, X = X_control, Y = delta_control)
  tauhat_control <- .predict(cate_model, cate_fit0, X)
  
  #Estimate propensity scores
  p <- propensity_score(data = data, x = x, y = y, w = w, model = propensity_model, plot = F)
  tauhat <- (1 - p) * tauhat_treat + p * tauhat_control
  
  #Plot histogram
  if (plot == T) {
    hist(tauhat, main = "CATE Estimated by X-Learner", xlab = "CATE", 
         ylab = "Frequency")
  }
  
  list(CATE = tauhat, ATE = mean(tauhat))
}
