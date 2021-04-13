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
  data[complete_cases,]
}

#' @description Helper function that calculates propensity scores with logistic 
#' regression
#'
#' @param data A dataframe object containing the variables and values
#' @param x A list of character vectors specifying variables to be included in 
#' the model (columns in the data). If unspecified, then it is assumed to be all 
#' columns in the data besides y and w
#' @param w A character vector specifying the treatment status. 
#' Treatment must be specified as 0 and 1 or TRUE and FALSE.
#' @param trim A list of two numeric elements describing cutoffs at which to 
#' trim propensity scores
#' 
#' @return a vector of propensity scores
#' 
.propensity_logistic <- function(data, x, w, trim){
  #Filter for NAs
  data <- .filter_nas(data)
  
  #Get formula and fit logistic regression
  formula <- as.formula(paste0(w, " ~ ", paste(x, collapse=' + ')))
  fit <- glm(formula, data, family = "binomial")
  
  #Predict probability
  p <- predict(fit, data, type = "response")
  p
}

#' @description Helper function that calculates propensity scores with 
#' Lasso-penalized logistic regression
#'
#' @param data A dataframe object containing the variables and values
#' @param x A list of character vectors specifying variables to be included in 
#' the model (columns in the data). If unspecified, then it is assumed to be all 
#' columns in the data besides y and w
#' @param w A character vector specifying the treatment status. 
#' Treatment must be specified as 0 and 1 or TRUE and FALSE.
#' @param trim A list of two numeric elements describing cutoffs at which to 
#' trim propensity scores
#' @param ... parameters for cv.glmnet
#' 
#' @return a vector of propensity scores
#' 
.propensity_lasso <- function(data, x, w, trim, ...){
  #Filter for NAs
  data <- .filter_nas(data)
  
  #Separate data
  treatment <- as.matrix(data[w])
  X <- model.matrix(~ . , data = data[x])
  
  #Fit LASSO model and predict
  lasso_fit <- glmnet::cv.glmnet(X, data[[w]], alpha = 1, ...)
  theta_hat <- predict(lasso_fit, X, s = "lambda.min")
  
  #Estimate propensity scores, subset to keep overlap assumption
  p <- (1 / (1 + exp(-theta_hat)))
  p
}

#' @description Helper function that calculates propensity scores with random 
#' forests
#'
#' @param data A dataframe object containing the variables and values
#' @param x A list of character vectors specifying variables to be included in 
#' the model (columns in the data). If unspecified, then it is assumed to be all 
#' columns in the data besides y and w
#' @param w A character vector specifying the treatment status. 
#' Treatment must be specified as 0 and 1 or TRUE and FALSE.
#' @param trim A list of two numeric elements describing cutoffs at which to 
#' trim propensity scores
#' @param ... Additional arguments for causal_forest
#' 
#' @return a vector of propensity scores.
#' 
.propensity_cf <- function(data, x, y, w, trim, ...) {
  #Filter for NAs
  data <- .filter_nas(data)
  
  #Separate into treatment, outcome, and covariates
  treatment <- data[[w]]
  X <- data[x]
  Y <- data[[y]]
  
  #Fit causal forest and estimate propensity scores
  cf_fit <- grf::causal_forest(X, Y, treatment, ...) 
  p <- cf_fit$W.hat
  names(p) <- 1:length(p)
  p
}

#' Check calibration and overlap assumptions via plot.
#' 
#' @description Creates a Q-Q plot and a histogram to check 
#' calibration and overlap assumptions, respectively.
#'
#' @param p a vector of propensity scores.
#' @param w a vector of treatment indicators. 
#' Treatment must be specified as 0 and 1 or TRUE and FALSE.
#' @param breaks an integer number of breaks for the histogram, default to 20.
#' 
#' @return a plot object that graphically tests overlap assumption and calibration.
#' 
#' @export
#' 
#' @examples 
#' data(lalonde)
#' p <- propensity_score(lalonde, y = "re78", w = "treat", model = "logistic", plot = FALSE)
#' w <- lalonde$treat
#' check_propensity(p = p, w = w)
#' 
#' @import graphics 
check_propensity <- function(p, w, breaks = 20){
  par(mfrow = c(1,2))
  hist(p, breaks = breaks, 
       main = "Propensity Scores Frequency", 
       xlab = "Propensity Score")
  {plot(smooth.spline(p, w, df = 4),
        main = "Propensity Score vs. Treatment", 
        xlab = "Propensity Score", ylab = "W",
        xlim = c(0, 1), ylim = c(0, 1))
    abline(0, 1)}
}

#' Calculate propensity scores.
#' 
#' @description Function that calculates propensity scores with an option to
#' check calibration and overlap assumptions
#'
#' @param data a dataframe object containing the variables and values.
#' @param x a list of character vectors specifying variables to be included in 
#' the model (columns in the data). If unspecified, then it is assumed to be all 
#' columns in the data besides y and w.
#' @param y a character vector specifying the response variable.
#' @param w a character vector specifying the treatment variable. 
#' Treatment must be specified as 0 and 1 or TRUE and FALSE.
#' @param subset an optional vector specifying a subset of observations to be
#' used in the fitting process.
#' @param model a character string naming the model used (one of "logistic", 
#' "lasso", or "causal forest").
#' @param trim a list of two numeric elements describing cutoffs at which to 
#' trim propensity scores.
#' @param plot logical; if TRUE then produces histogram of propensity scores to
#' check that the overlap assumption is satisfied and a Q-Q plot to check
#' the calibration.
#' @param ... additional arguments to be passed to the model fitting functions.
#' 
#' @return a vector of propensity scores where the indices of the vector correspond to the indices in \code{data}.
#' 
#' @references 
#' Imbens, Guido, and Donald Rubin. \emph{Causal Inference for Statistics, Social, and Biomedical Sciences: An Introduction}. 
#' Cambridge University Press, 2015.
#' 
#' @export
#' 
#' @examples
#' data(lalonde)
#' logit_propensities <- propensity_score(lalonde, y = "re78", w = "treat", model = "logistic")
#' lasso_propensities <- propensity_score(lalonde, y = "re78", w = "treat", model = "lasso")
#' cf_propensities <- propensity_score(lalonde, y = "re78", w = "treat", model = "causal forest")
#' 
#' @import stats
propensity_score <- function(data, 
                             x, 
                             y, 
                             w, 
                             subset = FALSE, 
                             model = "logistic", 
                             trim = FALSE, 
                             plot = TRUE, 
                             ...){
  #Filter for NAs
  data <- .filter_nas(data)
  
  #Check if missing
  if (missing(x)) {
    x <- setdiff(names(data), c(w, y))
  }
  
  #Check model type
  if (model == "logistic") {
    p <- .propensity_logistic(data, x, w, ...)
  }
  else if (model == "lasso") {
    p <- .propensity_lasso(data, x, w, ...)    
  }
  else if (model == "causal forest") {
    p <- .propensity_cf(data, x, y, w, ...)
  }
  else stop("model must be one of 'logistic', 'lasso', or 'causal forest'.")
  
  if (sum(trim) <= 0) {
    trim <- c(0, 1)
  }
  ind <- which(p > trim[1] & p < trim[2])
  
  #Check parameters
  if (plot == T) {
    #Check calibration/overlap
    check_propensity(p[ind], data[ind, w])
  }
  
  # return vector of propensities (indexed by index in `data`)
  p[ind]
}

#' @description Helper function to perform score matching
#'
#' @param diffs_set A matrix or dataframe giving some measure of distance between treatment and control points
#' @param dist The maximum distance between points to identify if they are matches.
#' @param repl logical; whether to match with or without replacement. 
#' 
#' @return a dataframe with two columns containing matched indexes. 
.matcher <- function(diffs_set, dist, repl){
  treat_indices <- as.numeric(colnames(diffs_set))
  control_indices <- vector()
  i <- 1
  
  for (treated in treat_indices) {
    if (min(diffs_set[, treated]) <= dist) {
      control_indices[i] <- as.numeric(names(which.min(diffs_set[, treated]))[1])
      # greedy approach
      if (repl == F) {
        diffs <- diffs[-which(rownames(diffs_set) == control_indices[i]), ]
      }
    } else {
      control_indices[i] <- NA
    }
    i <- i + 1
  }
  
  out <- data.frame(treat_indices, control_indices)
  out <- na.omit(out)
  out
}

#' Performs propensity score matching.
#'
#' @param data a dataframe object containing the variables and values.
#' @param w a character vector describing the treatment variable.
#' @param p a numeric vector of propensity scores.
#' @param max_distance the maximum distance between propensity scores to form a match. 
#' Default is 1.
#' @param replacement logical; if FALSE then treated and control units 
#' are matched without replacement of control units (a control unit can
#' only be matched once).
#'
#' @param type a string specifying which type of propensity score matching to perform. 
#' Choices are 'default', and 'linear'. Linear takes the logit of the propensity score before taking differences.
#' If not specified as linear, the default score matching will be used, which subtracts the scores as they are. 
#' 
#' @return a dataframe with two columns, the first corresponding to indices for treated observations, 
#' and the second corresponding to indices for matched control observations. 
#' 
#' @export
#' 
#' @examples 
#' data(lalonde)
#' 
#' p <- propensity_score(lalonde, y = "re78", w = "treat")
#' propensity_match(lalonde, w = "treat", p = p, max_distance = 1)

propensity_match <- function(data, w, p, max_distance = 1, replacement = TRUE, type = 'default'){
  # FILLER: UNTIL MAX_DISTANCE IS USED. HERE IN ORDER TO AVOID ERROR FOR 
  # UNUSED ARGUMENT max_distance ON FUNCITON CALL
  useless <- max_distance
  
  #Filter for NAs
  data <- .filter_nas(data)
  
  if (nrow(data) != length(p))
    stop("Data and propensity scores do not have the same length")
  
  if (length(setdiff(unique(data[[w]]), c(0, 1))) != 0)
    stop("Treatment is not 0 or 1.")
  
  if (type == 'linear'){
    p <- log(p/(1-p))
  }
  
  else if (type != 'default'){
    stop("propensity score matching type must be either default or linear")
  }
  
  p_treated <- p[data[w] == 1]
  p_control <- p[data[w] == 0]
  
  if (replacement == F & (length(p_treated) > length(p_control)))
    stop("When matching without replacement, number of treated units needs to be equal to or less than number of control units.")
  
  diffs <- abs(outer(p_control, p_treated,`-`))
  
  .matcher(diffs, max_distance, replacement)
}

#' Performs general matching with exact and mahalanobis distnace.
#'
#' @param data a dataframe object containing the variables and values.
#' @param w a character vector describing the treatment variable.
#' @param max_distance the maximum distance between propensity scores to form a match. 
#' Default is 0, as this is the distance for exact score matching
#' @param replacement logical; if FALSE then treated and control units 
#' are matched without replacement of control units (a control unit can
#' only be matched once).
#' @param type a string specifying which type of score matching to perform. 
#' Choices are 'exact', which is the default, and 'mahalanobis', which uses the Mahalanobis distance. 
#' If a different method is mentioned, then exact is performed by default. 
#' 
#' @return a dataframe with two columns, the first corresponding to indices for treated observations, 
#' and the second corresponding to indices for matched control observations. 
#' 
#' @export
#' 
#' @examples 
#' data(lalonde)
#' general_match(lalonde, w = "treat", max_distance = 10, type = 'mahalanobis')
general_match <- function(data, w, max_distance = 0, replacement = TRUE, type = 'exact'){
  
  #Helper function for mahalanobis distance"
  mlnobis <- function(x, cv, num_cols){
    vec <- x[1:num_cols] - x[(num_cols+1):(num_cols*2)]
    inv <- solve(cv)
    vec <- as.matrix(vec)
    t(vec) %*% inv %*% vec
  }
  
  if (length(setdiff(unique(data[[w]]), c(0, 1))) != 0)
    stop("Treatment is not 0 or 1.")
  
  treat <- data[data[w] == 1,]
  treat <- treat[,!colnames(treat) == w]
  cont <- data[data[w] == 0,]
  cont <- cont[,!colnames(cont) == w]
  
  comb <- merge(treat, cont, by = NULL)
  
  if (replacement == F & (length(treat) > length(cont)))
    stop("When matching without replacement, number of treated units needs to be equal to or less than number of control units.")
  
  col_cnt <-ncol(data) - 1
  diffs <- 0
  if (type == 'mahalanobis'){
    cv <- cov(cont)
    diffs <- apply(comb, 1, function(x)mlnobis(x, cv, col_cnt))
  }
  
  else if (type == 'exact'){
    diffs <- apply(comb, 1, function(x)as.numeric(!all(x[1:col_cnt]==x[(col_cnt+1):(col_cnt*2)])))
  }
  
  else{
    stop("type must be either exact or mahalanobis")
  }
  
  diffs <- matrix(diffs, nrow = nrow(cont))
  rownames(diffs) <- seq(1, nrow(diffs))
  colnames(diffs) <- seq(1, ncol(diffs))

  if (type == 'exact'){
    max_distance <-  0
  }
  .matcher(diffs, max_distance, replacement)
}

#' Assess balance in multivariate covariate distributions directly or by using 
#' propensity scores.
#'
#' @param data a dataframe object containing the variables and values.
#' @param x a list of character vectors specifying the covariates.
#' @param w a character vector specifying the treatment variable. 
#' Treatment must be specified as 0 and 1 or TRUE and FALSE.
#' @param p a numeric vector of propensity scores. Only specified if \code{method = "NLPD"}.
#' @param method a character vector specifying the which method to use to assess balance 
#' (one of "Mahalanobis distance" or "NLPD")
#' 
#' @return a numeric value for difference between the covariate distributions for treated and control groups.
#' 
#' @export
#' 
#' @examples
#' data(lalonde)
#' 
#' vars <- names(lalonde)
#' covariates <- vars[!vars %in% c("re78", "treat")]
#' 
#' assess_covariate_balance(lalonde, x = covariates, w = "treat")
#' p <- propensity_score(lalonde, y = "re78", w = "treat")
#' 
#' assess_covariate_balance(lalonde, x = covariates, w = "treat", p = p, method = "NLPD")
assess_covariate_balance <- function(data, x, w, p = NULL, method = "Mahalanobis distance"){
  #Filter for NAs
  data <- .filter_nas(data)
  
  if (method == "Mahalanobis distance") {
    X <- data[x]
    control <- X[data[w] == 1, ]
    treated <- X[data[w] == 0, ]
    
    cov_c <- cov(control, use = 'pairwise.complete')
    cov_t <- cov(treated, use = 'pairwise.complete')
    
    diff_in_means <- colMeans(treated) - colMeans(control)
    diff <- sqrt( abs( t(diff_in_means) %*% ((cov_t + cov_c) / 2)^(-1) %*% diff_in_means ) )
    diff <- diff[[1, 1]]
  } else if (method == "NLPD") {
    #NLPD = Normalized linearized propensity difference
    linear_p <- log( p / (1 - p) )
    
    linear_p_treated <- linear_p[data[w] == 1]
    linear_p_control <- linear_p[data[w] == 0]
    
    var_t <- var(linear_p_treated)
    var_p <- var(linear_p_control)
    diff <- (mean(linear_p_treated) - mean(linear_p_control)) / 
      sqrt( (var(linear_p_treated) + var(linear_p_control)) / 2 )
  } else{
    stop("method not valid")
  }
  
  return(diff)
}

#'
#' @param data a dataframe object containing the variables and values.
#' @param x a list of character vectors specifying the covariates.
#' @param w a character vector describing the treatment variable. 
#' Treatment must be specified as 0 and 1 or TRUE and FALSE.
#' @param method a character vector specifying the method for which to assess balance in 
#' covariate distributions across treatment and control groups. 
#' (one of "Normalized difference" or "Dispersion")
#' 
#' @return a dataframe object containing the difference (measured differenty with different methods)
#' across covariates with their respective names
#' 
.get_diff_by_covar <- function(data, x, w, method){
  
  X <- data[x]
  control <- X[data[w] == 0, ]
  treated <- X[data[w] == 1, ]
  
  var_c <- apply(control, 2, var)
  var_t <- apply(treated, 2, var)
  
  if (method == "Normalized difference") {
    means_c <- colMeans(control)
    means_t <- colMeans(treated)
    
    diff <- (means_t - means_c) / sqrt( ( var_t + var_c) / 2 )
  } else if (method == "Dispersion") {
    diff <- log(sqrt(var_t)) - log(sqrt(var_c))
  }
  Covariates <- names(diff)
  data.frame(Covariates, diff)
}


#' @param data a dataframe object containing the variables and values.
#' @param x a list of character vectors specifying the covariates.
#' @param w a character vector describing the treatment variable.
#' Treatment must be specified as 0 and 1 or TRUE and FALSE.
#' @param method a character vector specifying the method for which to assess balance in 
#' covariate distributions across treatment and control groups. 
#' (one of "Normalized difference" or "Dispersion")
#' @param threshold the threshold to plot. Default is 0.2.
#' 
#' @return a \code{ggplot} object that shows balance by covariate
#' 
.balance_wo_matching <- function(data, x, w, method, threshold){
  
  df_plot <- .get_diff_by_covar(data, x, w, method)
  
  Covariates <- NULL
  g <- ggplot(df_plot, aes(x = diff, y = Covariates)) +
    geom_point() +
    xlab(method) +
    geom_vline(xintercept = 0) +
    geom_vline(xintercept = threshold, lty = "dashed") +
    geom_vline(xintercept = -threshold, lty = "dashed") 
  g
}

#' @description Assess balance in covariate distribution with propensity score matching 
#' and assess ability to adjust for differences in covariates by treatment status.
#'
#' @param data a dataframe object containing the variables and values.
#' @param x a list of character vectors specifying the covariates.
#' @param w a character vector describing the treatment variable.
#' Treatment must be specified as 0 and 1 or TRUE and FALSE.
#' @param method a character vector specifying the method for which to assess balance in 
#' covariate distributions across treatment and control groups.
#' (one of "Normalized difference" or "Dispersion")
#' @param threshold the threshold to plot. Default is 0.2.
#' 
#' @return a \code{ggplot} object that shows balance by covariate with matching
#' 
.balance_w_matching <- function(data, x, w, matched_indices, method, threshold, colors){
  
  matched_indices <- c(matched_indices[[1]], matched_indices[[2]])
  matched <- ifelse(as.numeric(rownames(data)) %in% matched_indices, 1, 0)
  
  diff_matched <- .get_diff_by_covar(data[matched == 1, ], x, w, method)
  diff_matched$matched <- 1
  diff_unmatched <- .get_diff_by_covar(data[matched == 0, ], x, w, method)
  diff_unmatched$matched <- 0
  
  df_plot <- rbind(diff_matched, diff_unmatched)
  
  if (is.null(colors)) colors <- c("#ff4d47", "#74caff")
  
  Covariates <- NULL
  g <- ggplot(df_plot,
              aes(
                x = diff,
                y = Covariates,
                group = as.factor(matched),
                col = as.factor(matched)
              )) +
    geom_point() +
    xlab(method) +
    geom_vline(xintercept = 0) +
    geom_vline(xintercept = threshold, lty = "dashed") +
    geom_vline(xintercept = -threshold, lty = "dashed") +
    guides(col = guide_legend(title = "Matched")) +
    scale_fill_manual(values = colors) +
    theme_minimal()
  g
}

#' Assess balance in covariate distribution graphically.
#' 
#' @description Assess balance in covariate distributio and assess ability to adjust 
#' for differences in covariates by treatment status.
#'
#' @param data a dataframe object containing the variables and values.
#' @param x a list of character vectors specifying the covariates.
#' @param w a character vector describing the treatment variable.
#' @param matched_indices a data frame of two columns, one of treatment indices and one of 
#' matched (produced by \code{propensity_match}).
#' @param method a character vector specifying the method for which to assess balance in 
#' covariate distributions across treatment and control groups.
#' @param threshold the threshold to plot. Default is 0.2.
#' @param colors an optional vector of two specifying the colors for the control and
#' treatment groups, respectively.
#' @param title an optional plot title.
#' 
#' @return a \code{ggplot} object.
#' 
#' @export
#' 
#' @examples 
#' data(lalonde)
#' 
#' p <- propensity_score(lalonde, y = "re78", w = "treat")
#' matched_indices <- propensity_match(lalonde, w = "treat", p = p, max_distance = .00001)
#' 
#' vars <- names(lalonde)
#' covariates <- vars[!vars %in% c("re78", "treat")]
#' 
#' balance_plot(lalonde, x = covariates, w = "treat", matched_indices = matched_indices)
#' 
#' @importFrom ggplot2 ggplot aes geom_point xlab ylab guides guide_legend scale_fill_manual theme_minimal
#' geom_density ggtitle geom_vline
balance_plot <- function(data, 
                         x, 
                         w, 
                         matched_indices = NULL, 
                         method = "Normalized difference", 
                         threshold = 0.2, 
                         colors = NULL, 
                         title = NULL){
  #Filter for NAs
  data <- .filter_nas(data)
  
  if (is.null(matched_indices)) {
    g <- .balance_wo_matching(data, x, w, method, threshold)
  } else {
    g <- .balance_w_matching(data, x, w, matched_indices, method, threshold, colors)
  }
  
  if (!is.null(title)) g <- g + ggtitle(title)
  g
}

#' Assess balance for a single covairate.
#'
#' @param data a dataframe object containing the variables and values.
#' @param covariate a character vectors specifying the covariate
#' @param w a character vector specifying the treatment variable.
#' Treatment must be specified as 0 and 1 or TRUE and FALSE.
#' @param colors an optional vector of two specifying the colors for the control and
#' treatment groups, respectively.
#' @param title an optional plot title.
#' 
#' @return a \code{ggplot} object.
#' 
#' @export
#' 
#' @examples
#' data(lalonde)
#' 
#' # histogram for continues variable: 
#' univariate_balance_plot(lalonde, "age", "treat")
#' # bar plot for discrete factor variable:
#' lalonde$hisp <- factor(lalonde$hisp)
#' univariate_balance_plot(lalonde, "hisp", "treat")
#' 
#' @importFrom ggplot2 ggplot aes geom_bar xlab ylab guides guide_legend scale_fill_manual theme_minimal
#' geom_density ggtitle stat
univariate_balance_plot <- function(data, covariate, w, colors = NULL, title = NULL){
  #Filter for NAs
  data <- .filter_nas(data)
  
  if (length(setdiff(unique(data[[w]]), c(0, 1))) != 0)
    stop("Treatment is not 0 or 1.")
  
  if (!is.factor(data[[covariate]]) & length(unique(data[[covariate]])) < 5) {
    warning("Covariate has fewer than 5 unique values. Consider turning variable into a factor to construct a bar plot.")
  }
  
  if (is.null(colors)) colors <- c("#ff4d47", "#74caff")
  
  if (is.factor(data[[covariate]])) {
    prop <- NULL
    g <- ggplot(data,
                aes(
                  x = !!as.symbol(covariate),
                  group = as.factor(!!as.symbol(w)),
                  fill = as.factor(!!as.symbol(w))
                )) +
      geom_bar(aes(y = stat(prop)), position = "dodge") +
      xlab(covariate) +
      ylab("Percent") +
      guides(fill = guide_legend(title = "Treatment")) +
      scale_fill_manual(values = colors) +
      theme_minimal()
  } else {
    g <- ggplot(data,
                aes(
                  x = !!as.symbol(covariate),
                  group = as.factor(!!as.symbol(w)),
                  fill = as.factor(!!as.symbol(w))
                )) +
      geom_density(alpha = .5) +
      ylab("Density") +
      guides(fill = guide_legend(title = "Treatment")) +
      scale_fill_manual(values = colors) +
      theme_minimal()
  }
  
  if (!is.null(title)) g <- g + ggtitle(title)
  
  g
}
