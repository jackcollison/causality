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

#' @description Helper function to return ATT with and without confidence interval
#'
#' @param ATT the average treatment effect as calculated with another function.
#' @param se_hat estimate of the standard error.
#' @param cf logical; if TRUE, then include confidence interval on ATT.
#' 
#' @return a list of ATT, 95 percent confidence interval upperbound and lowerbound or just ATT, depending on user input of \code{cf}.
#' 
.att_return_helper <- function(ATT, se_hat, cf){
  if (cf == T) {
    #Confidence intervals
    lb <- ATT - 1.96 * se_hat
    ub <- ATT + 1.96 * se_hat
    return(c(ATT = ATT, Lowerbound = lb, Upperbound = ub))
  } else {
    return(ATT)
  }
}

#' Uses random forests to naively estimate the average 
#' treatment effect on the treated (ATT) without weighting
#' 
#' @param data a dataframe object containing the variables and values.
#' @param x a list of character vectors specifying variables to be included in 
#' the model (columns in the data). If unspecified, then it is assumed to be all 
#' columns in the data besides y and w.
#' @param y a character vector specifying the response variable.
#' @param w a character vector specifying the treatment status.
#' @param cf logical; if TRUE then includes confidence interval on ATT.
#' @param ... additional arguments to causal_forest.
#' 
#' @details Computes an estimate of the ATT \eqn{\tau_T} with a naive estimate on just the treated group
#' (see naive_ate for more details).
#' 
#' @return a list of ATT, 95 percent confidence interval upperbound and lowerbound or just ATT, depending on user input of \code{cf}.
#' 
#' @references Athey, Susan, Imbens, Guido, Pham, Thai, and Wager, Stefan. 2017. ``Estimating Average Treatment Effects: Supplementary Analyses and Remaining Challenges." 
#' \emph{The American Economic Review}, Vol. 107(5). pgs. 278–281.
#' url{www.jstor.org/stable/44250405}
#' 
#' @export
#' 
#' @examples
#' data("lalonde")
#' att <- rf_att(data = lalonde, y = "re78", w = "treat", num.trees = 100, mtry = 3)
#' 
rf_att <- function(data, x, y, w, cf = TRUE, ...) {
  #Filter for NAs
  data <- .filter_nas(data)
  
  #Separate data into treatment, response, and covariates
  W <- data[[w]]
  Y <- data[[y]]
  
  if (missing(x)) x <- setdiff(names(data), c(w, y))
  X <- data[, x]
  
  #Fit causal forest
  cf_fit <- grf::causal_forest(X, Y, W, ...) 
  
  #Calculate ATE, estimate standard errors, and return
  att <- grf::average_treatment_effect(cf_fit, target.sample = "treated") 
  ATT <- unname(att["estimate"])
  se_hat <- unname(att["std.err"])
  
  .att_return_helper(ATT, se_hat, cf)
}

#' Estimate the average treatment effect on the treated (ATT) via augmented inverse propensity 
#' weighting.
#' 
#' @param data a dataframe object containing the variables and values.
#' @param x a list of character vectors specifying variables to be included in 
#' the model (columns in the data). If unspecified, then it is assumed to be all 
#' columns in the data besides y and w.
#' @param y a character vector specifying the response variable.
#' @param w a character vector specifying the treatment status.
#' @param p a vector containing propensity score values.
#' @param cf logical; if TRUE then includes confidence interval on ATT.
#' @param ... additional arguments to causal_forest.
#' 
#' @details Computes an estimate of the ATT \eqn{\tau_T} with a augmented inverse propensity score weighted 
#' estimate on just the treated group (see aipw_ate for more details).
#' 
#' @return a list of ATT, 95 percent confidence interval upperbound and lowerbound or just ATT, depending on user input of \code{cf}.
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
#' att <- aipw_rf_att(data = lalonde, y = "re78", w = "treat", p = p, num.trees = 100, mtry = 3)
#' 
aipw_rf_att <- function(data, x, y, w, p, cf = TRUE, ...) {
  #Filter for NAs
  data <- .filter_nas(data)
  
  if (nrow(data) != length(p)) {
    stop("Data and propensity scores do not have the same length")
  }
  
  #Separate data into treatment, response, and covariates
  #Create variables for brevity
  W <- data[[w]]
  Y <- data[[y]]
  
  if (missing(x)) x <- setdiff(names(data), c(w, y))
  X <- data[, x]
  
  #Fit causal forest
  cf_fit <- grf::causal_forest((p / (1 - p))*X, Y, W, ...)
  
  #Calculate ATE, estimate standard errors, and return
  att <- grf::average_treatment_effect(cf_fit, target.sample = "treated") 
  ATT <- unname(att["estimate"])
  se_hat <- unname(att["std.err"])
  
  .att_return_helper(ATT, se_hat, cf)
}

