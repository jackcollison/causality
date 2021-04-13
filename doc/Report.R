## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)
library(causality)

## -----------------------------------------------------------------------------
#Import data
data(lalonde)

#Continuous variable
univariate_balance_plot(lalonde, "age", "treat")

#Factor
lalonde$hisp <- factor(lalonde$hisp)
univariate_balance_plot(lalonde, "hisp", "treat")

## -----------------------------------------------------------------------------
#Import data
data(lalonde)

#Get covariates of interest
vars <- names(lalonde)
covariates <- vars[!vars %in% c("re78", "treat")]

#Assess balance
assess_covariate_balance(lalonde, x = covariates, w = "treat")

## -----------------------------------------------------------------------------
#Calculate propensity scores
p <- propensity_score(lalonde, y = "re78", w = "treat", plot = F)

#Assess balance
assess_covariate_balance(lalonde, x = covariates, w = "treat", p = p, method = "NLPD")

## -----------------------------------------------------------------------------
#Balance plot without matching
balance_plot(lalonde, x = covariates, w = "treat")

#Matched indicies
matched_indices <- propensity_match(lalonde, w = "treat", p = p, max_distance = .00001)

#Balance plot with matching
balance_plot(lalonde, x = covariates, w = "treat", matched_indices = matched_indices)

## -----------------------------------------------------------------------------
#Estimate propensity scores
p <- propensity_score(lalonde, y = "re78", w = "treat", model = "logistic", plot = F)
head(p)

## -----------------------------------------------------------------------------
#Estimate propensity scores
p <- propensity_score(lalonde, y = "re78", w = "treat", model = "lasso", plot = F)
head(p)

## -----------------------------------------------------------------------------
#Estimate propensity scores
p <- propensity_score(lalonde, y = "re78", w = "treat", model = "causal forest", plot = F)
head(p)

## ---- fig.height = 4, fig.width = 7-------------------------------------------
#Calculate propensity scores, load treatment
p <- propensity_score(lalonde, y = "re78", w = "treat", model = "logistic", plot = FALSE)
w <- lalonde$treat

#Check assumptions
check_propensity(p = p, w = w)

## -----------------------------------------------------------------------------
#Calculate propensity scores
p <- propensity_score(lalonde, y = "re78", w = "treat", plot = F)

#Match propensity scores
head(propensity_match(lalonde, w = "treat", p = p, max_distance = .2))

## -----------------------------------------------------------------------------
#General matching (FIX diffs2 bug)
head(general_match(lalonde, w = "treat", max_distance = 10, type = "mahalanobis"))

## -----------------------------------------------------------------------------
#Estimate ATE
naive_ate(data = lalonde, y = "re78", w = "treat")

## -----------------------------------------------------------------------------
#Calculate propensity scores and IPW ATE
p <- propensity_score(data = lalonde, y = "re78", w = "treat", model = "logistic", plot = F)
ipw_ate(data = lalonde, y = "re78", w = "treat", p = p)

## -----------------------------------------------------------------------------
#Calculate propensity scores and IPW OLS
p <- propensity_score(lalonde, y = "re78", w = "treat", model = "logistic", plot = F)
prop_weighted_ols_ate(data = lalonde, y = "re78", w = "treat", p = p)

## -----------------------------------------------------------------------------
#Calculate propensity scores and AIPW ATE
p <- propensity_score(lalonde, y = "re78", w = "treat", model = "logistic", plot = F)
aipw_ate(data = lalonde, y = "re78", w = "treat", p = p)

## -----------------------------------------------------------------------------
#Calculate propensity scores and ATE
p <- propensity_score(lalonde, y = "re78", w = "treat", model = "logistic", plot = F)
propensity_strat_ate(data = lalonde, 
                     y = "re78", 
                     w = "treat",
                     p = p, 
                     n_strata = 4, 
                     min_obs = 1)

## -----------------------------------------------------------------------------
#Estimate ATE
double_selection_ate(data = lalonde, y = "re78", w = "treat")

## -----------------------------------------------------------------------------
#Estimate ATT
rf_att(data = lalonde, y = "re78", w = "treat", num.trees = 100, mtry = 3)

## -----------------------------------------------------------------------------
#Calculate propensity scores, estimate ATT
p <- propensity_score(lalonde, y = "re78", w = "treat", model = "logistic", plot = F)
aipw_rf_att(data = lalonde, y = "re78", w = "treat", p = p, num.trees = 100, mtry = 3)

## ---- fig.height = 4, fig.width = 7-------------------------------------------
#Estimate HTE
head(s_learner(data = lalonde, y = "re78", w = "treat", num.trees = 100, mtry = 3)$CATE)

## ---- fig.height = 4, fig.width = 7-------------------------------------------
#Estimate HTE
head(t_learner(data = lalonde, y = "re78", w = "treat", num.trees = 100, mtry = 3)$CATE)

## ---- fig.height = 4, fig.width = 7-------------------------------------------
#Estimate HTE
head(x_learner(data = lalonde, y = "re78", w = "treat", num.trees = 100, mtry = 3)$CATE)

