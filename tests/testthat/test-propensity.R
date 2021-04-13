library(testthat)
library(ggplot2)
library(grf)
context("propensity")
source('../../R/propensity.R')
load("../../data/lalonde.RData")

test_that("propensity_score returns error when model is not one of lasso, logistic, causal forrest", {
  expect_error(propensity_score(lalonde, y = "re78", w = "treat", model = "test"))
})

test_that("propensity_score runs without any errors for all three model types", {
  expect_silent(propensity_score(lalonde, y = "re78", w = "treat", model = "logistic"))
  expect_silent(propensity_score(lalonde, y = "re78", w = "treat", model = "lasso"))
  propensity_score(lalonde, y = "re78", w = "treat", model = "causal forest")
})

test_that("propensity_match returns error when data and scores have size mismatch", {
  p <- propensity_score(lalonde, y = "re78", w = "treat")
  expect_error(propensity_match(lalonde, w = "treat", p = seq(1,nrow(lalonde)+5), max_distance = 1))
})

test_that("propensity_match returns error when treatment column is wrong", {
  p <- propensity_score(lalonde, y = "re78", w = "treat")
  expect_error(propensity_match(lalonde, w = "age", p = p, max_distance = 1))
})

test_that("propensity_match returns error when replacement is false and more control than treated units", {
  p <- propensity_score(lalonde, y = "re78", w = "treat")
  d2 = lalonde
  d2$treat = as.numeric(!d2$treat)
  expect_error(propensity_match(d2, w = "treat", p = p, max_distance = 1, replacement = FALSE))
})

test_that("propensity_match returns error when something other than default or linear is specified", {
  p <- propensity_score(lalonde, y = "re78", w = "treat")
  expect_error(propensity_match(lalonde, w = "treat", p = p, max_distance = 1, type = 'test'))
})

test_that("propensity_match returns error when treatment column is wrong", {
  expect_error(general_match(lalonde, w = "age", max_distance = 1))
})

test_that("propensity_match returns error when replacement is false and more control than treated units", {
  d2 = lalonde
  d2$treat = as.numeric(!d2$treat)
  expect_error(general_match(d2, w = "treat", max_distance = 1, replacement = FALSE))
})

test_that("propensity_match runs for both default and linear matching", {
  p <- propensity_score(lalonde, y = "re78", w = "treat")
  expect_silent(propensity_match(lalonde, w = "treat", p = p, max_distance = 1, type = 'linear'))
  expect_silent(propensity_match(lalonde, w = "treat", p = p, max_distance = 1, type = 'default'))
})

test_that("general_match returns error when something other than default or exact or mahalanobis", {
  expect_error(general_match(lalonde, w = "treat", max_distance = 1, type = 'test'))
})

test_that("general_match runs for both exact and mahalanobis", {
  expect_silent(general_match(lalonde, w = "treat", max_distance = 1, type = 'exact'))
  expect_silent(general_match(lalonde, w = "treat", max_distance = 1, type = 'mahalanobis'))
})

test_that("univariate_variance_plot returns error when treatment is not 0 or 1", {
  expect_error(univariate_balance_plot(lalonde, "age", "age"))
})

test_that("univariate variance plot shows warning when covariate is not factor and has fewer than 5 unique values", {
  d2 = lalonde
  d2$dummy = rep(1, nrow(d2))
  expect_warning(univariate_balance_plot(d2, "dummy", "treat"))
})

test_that("assess covariate balance returns error for invalid method", {
  vars <- names(lalonde)
  covariates <- vars[!vars %in% c("re78", "treat")]
  expect_error(assess_covariate_balance(lalonde, x = covariates, w = "treat", method = "test"))
})

test_that("assess_covariate_balance works for both methods",{
  p <- propensity_score(lalonde, y = "re78", w = "treat")
  vars <- names(lalonde)
  covariates <- vars[!vars %in% c("re78", "treat")]
  expect_silent(assess_covariate_balance(lalonde, x = covariates, w = "treat", p = p, method = "NLPD"))
  expect_silent(assess_covariate_balance(lalonde, x = covariates, w = "treat", method = "Mahalanobis distance"))
})