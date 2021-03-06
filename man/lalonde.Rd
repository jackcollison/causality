% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/lalonde.R
\docType{data}
\name{lalonde}
\alias{lalonde}
\title{Lalonde Data}
\format{A data frame with 445 observations and 12 variables, described below:
\describe{
  \item{age}{age in years}
  \item{educ}{years of schooling}
  \item{black}{indicator variable for Black}
  \item{hisp}{indicitaor variable for Hispanic}
  \item{married}{indicator variable for marital status}
  \item{nodegr}{indicator variable for high school diploma}
  \item{re74}{real earnings in 1974}
  \item{re75}{real earnings in 1975}
  \item{re78}{real earnings in 1978}
  \item{u74}{indicator variable for 1974 earnings being zero}
  \item{u75}{indicator variable for 1975 earnings being zero}
  \item{treat}{indicator variable for treatment status}

}}
\source{
\href{https://users.nber.org/~rdehejia/data/.nswdata2.html}{NBER}
}
\usage{
data(lalonde)
}
\description{
Data from a paper by Robert Lalonde, "Evaluating the Econometric 
Evaluations of Training Programs," \emph{American Economic Review}, Vol. 76, 
pp. 604-620.
}
\details{
We use data from Lalonde's economic examination of training programs in which he estimates the effect of
an employment program designed to help disadvantaged workers enter the labor market (Lalonde, 1986). The
data comes from the National Supported Work Demonstration (NSW), a temporary employment program that aimed
to help unskilled workers find jobs. At the time of the study, the NSW program was operating in ten cities
across the country, admitting women in families under Aid to Families with Dependent Children (AFDC),
ex-drug addicts, ex-criminal offenders, and high school dropouts. Individuals were randomly assigned to a
treatment or control group where the treated were guaranteed jobs for 9-18 months (depending on the target
group and site). The treatment group was divided into sub-groups of three to five individuals who typically
worked together and met with an NSW representative to discuss work performance and issues. There is within-site
and between-site variation (e.g. male and female participants usually had different types of work and the
different sites had different base wages). Earnings and demographic data was collected for both the treatment
and control groups, both during and at nine-month intervals after the conclusion of the program. Although there
was a high attrition rate, it did not affect the experimental design. Lalonde provides experimental and non-
experimental examinations of the data. For our purposes, we use the data in a non-experimental context.
}
\examples{
data(lalonde)
treat <- attr(lalonde, "treat")
age <- lalonde$age
}
\references{
Robert Lalonde, "Evaluating the Econometric Evaluations of Training Programs," \emph{American Economic Review}, Vol. 76, pp. 604-620.
}
\keyword{datasets}
