% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ate.R
\name{aipw_ate}
\alias{aipw_ate}
\title{Estimate average treatment effect (ATE) with
augmented inverse propensity score weighting.}
\usage{
aipw_ate(data, x, y, w, p, cf = TRUE)
}
\arguments{
\item{data}{a dataframe object containing the variables and values.}

\item{x}{a list of character vectors specifying variables to be included in 
the model (columns in the data). If unspecified, then it is assumed to be all 
columns in the data besides y and w.}

\item{y}{a character vector specifying the response variable.}

\item{w}{a character vector specifying the treatment status.}

\item{p}{a vector containing propensity score values.}

\item{cf}{logical; if TRUE then includes confidence interval on ATE.}
}
\value{
a list of ATE, 95 percent confidence interval upperbound and lowerbound or just ATE, depending on user input of \code{cf}.
}
\description{
Estimate average treatment effect (ATE) with
augmented inverse propensity score weighting.
}
\details{
Estimates the ATE \eqn{\tau} using the doubly robust method described in Scharfstein, Robins and Rotznitzky (1998) that combines
both regression and propensity score weighting. 
\deqn{\tau = E \left[  W_i \frac{Y_i-\tau(1,X_i)}{e(X_i)} + (1-W_i) \frac{Y_i-\tau(0,X_i)}{1-e(X_i)} + \tau(1,X_i) - \tau(0,X_i)\right]}
where \eqn{e(X_i) = P(W_i = 1 | X_i = x)} is the propensity score and \eqn{\tau(1, X_i)} and \eqn{\tau{0, X_i}} are estimated in
the first stage via OLS regression.
}
\examples{
data("lalonde")

# calculate propensity scores
p <- propensity_score(lalonde, y = "re78", w = "treat", model = "logistic")
ate <- aipw_ate(data = lalonde, y = "re78", w = "treat", p = p)

}
\references{
Rotnitzky, Andrea, James M. Robins, and Daniel O. Scharfstein. 1998. ``Semiparametric 
Regression for Repeated Outcomes with Nonignorable Nonresponse." \emph{Journal of the American Statistical Association}. 
Vol. 93, No. 444, Dec. pgs. 1321-1339.

Cao, Weihua, Anastasios A. Tsiatis, and Marie Davidian. 2009. ``Improving efficiency and robustness of the doubly 
robust estimator for a population mean with incomplete data". \emph{Biometrika}. Vol. 96(3). pgs. 723–734.
}
