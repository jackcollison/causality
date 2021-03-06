% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/hte.R
\name{s_learner}
\alias{s_learner}
\title{Estimate heterogeneous treatment effects (HTEs)
using the S-Learner strategy.}
\usage{
s_learner(data, x, y, w, base_learner = "regression forest", plot = TRUE, ...)
}
\arguments{
\item{data}{a dataframe object containing the variables and values.}

\item{x}{a list of character vectors specifying variables to be included in 
the model (columns in the data). If unspecified, then it is assumed to be all 
columns in the data besides y and w.}

\item{y}{a character vector specifying the response variable.}

\item{w}{a character vector specifying the treatment status.}

\item{base_learner}{a character vector specifying the base learner. One of "regression forest"
or "OLS". Default is "regression forest".}

\item{plot}{logical; if TRUE, then plots a histogram of treatment effects.}

\item{...}{additional arguments to the base learner.}
}
\value{
a list of two. The first element is a vector of conditional average treatment effect for each observation. 
The second element is the estimated average treatment effect.
}
\description{
Estimate heterogeneous treatment effects (HTEs)
using the S-Learner strategy.
}
\details{
Implements the S-learner algorithm described in Künzel et al. (2019) for estimating
conditional average treatment effects (CATE). In the S-learner algorithm, the treatment \eqn{W}
is included as a feature similar to all of the other covariates without the indicator being given any special role.
The combined response function 
\deqn{\mu(x, w) = E \left[ Y^{obs} | X = x, W = w \right]}
can then be estimated using any base learner (supervised machine learning or regression algorithm) on the entire dataset. 
Here we implement the S-learner with the option forA linear regression or a regression forest (see Athey, Tibshirani, and Wager (2016)) 
as the base learner.

The CATE estimator is then given by 
\deqn{\hat{\tau}(x) = \hat{\mu}(x, 1) - \hat{\mu}(x, 0).}
}
\examples{
data("lalonde")
hte <- s_learner(data = lalonde, y = "re78", w = "treat", num.trees = 100, mtry = 3)

}
\references{
Künzel, Sören R., Jasjeet S. Sekhon, Peter J. Bickel, and Bin Yu. 2019. 
``Metalearners for estimating heterogeneous treatment effects using machine learning." 
\emph{Proceedings of the National Academy of Sciences of the United States of America.}
Mar. 116(10): 4156–4165. \url{https://doi.org/10.1073/pnas.1804597116}

Athey, Susan, Julie Tibshirani, and Stefan Wager. 2016. ``Generalized Random Forests."
Working paper; Forthcoming in the Annals of Statistics. \url{https://arxiv.org/abs/1610.01271}
}
\seealso{
\code{\link{t_learner}}, \code{\link{x_learner}}
}
