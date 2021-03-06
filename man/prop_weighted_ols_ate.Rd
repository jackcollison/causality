% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ate.R
\name{prop_weighted_ols_ate}
\alias{prop_weighted_ols_ate}
\title{Estimate average treatment effect (ATE) with
inverse propensity score weighting using OLS.}
\usage{
prop_weighted_ols_ate(data, y, w, p, cf = TRUE)
}
\arguments{
\item{data}{a dataframe object containing the variables and values.}

\item{y}{a character vector specifying the response variable.}

\item{w}{a character vector specifying the treatment status.}

\item{p}{a vector containing propensity score values.}

\item{cf}{logical; if TRUE then includes confidence interval on ATE.}
}
\value{
a list of ATE, 95 percent confidence interval upperbound and lowerbound or just ATE, depending on user input of \code{cf}
}
\description{
Estimate average treatment effect (ATE) with
inverse propensity score weighting using OLS.
}
\details{
Computes an estimate of the ATE \eqn{\tau} via weighted OLS regression using propensity score weighting
The weights are given by:
\deqn{w_i = \frac{W_i}{e(X_i)} - \frac{1 - W_i}{1 - e(X_i)}}
where \eqn{e(X_i) = P(W_i = 1 | X_i = x)} is the propensity score.
}
\examples{
data("lalonde")

# calculate propensity scores
p <- propensity_score(lalonde, y = "re78", w = "treat", model = "logistic")
ate <- prop_weighted_ols_ate(data = lalonde, y = "re78", w = "treat", p = p)

}
\references{
Aronow, Peter M.; Samii, Cyrus. ``Estimating average causal effects under general interference, with application to a social network experiment." 
\emph{Annals of Applied Statistics} 11 (2017), no. 4, 1912--1947
\url{https://arxiv.org/pdf/1305.6156.pdf}
}
