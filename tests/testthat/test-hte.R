context("hte")
library(testthat)
source('../../R/hte.R')
load("../../data/lalonde.RData")

test_that("s_learner returns error if base learner invalid", {
  expect_error(s_learner(data = lalonde, y = "re78", w = "treat", base_learner = "test"))
})

test_that("s_learner runs for both base learners", {
  expect_silent(s_learner(data = lalonde, y = "re78", w = "treat", base_learner = "regression forest", num.trees = 100, mtry = 3))
  expect_silent(s_learner(data = lalonde, y = "re78", w = "treat", base_learner = "OLS"))
})

test_that("t_learner returns error if base learner invalid", {
  expect_error(t_learner(data = lalonde, y = "re78", w = "treat", base_learner = "test"))
})

test_that("t_learner runs for both base learners", {
  expect_silent(t_learner(data = lalonde, y = "re78", w = "treat", base_learner = "regression forest", num.trees = 100, mtry = 3))
  expect_silent(t_learner(data = lalonde, y = "re78", w = "treat", base_learner = "OLS"))
})

test_that("x_learner returns error if base learner invalid", {
  expect_error(x_learner(data = lalonde, y = "re78", w = "treat", base_learner = "test"))
})

test_that("x_learner runs without errors for both base learners", {
  expect_silent(x_learner(data = lalonde, y = "re78", w = "treat", base_learner = "regression forest", num.trees = 100, mtry = 3))
  expect_warning(x_learner(data = lalonde, y = "re78", w = "treat", base_learner = "OLS"))
})

test_that("x_learner returns error if base cate model invalid", {
  expect_error(x_learner(data = lalonde, y = "re78", w = "treat", base_learner = "OLS", cate_model = "test"))
})

test_that("x_learner runs without errors for both cate models", {
  expect_silent(x_learner(data = lalonde, y = "re78", w = "treat", base_learner = "regression forest", cate_model = "regression forest"))
  expect_warning(x_learner(data = lalonde, y = "re78", w = "treat", base_learner = "regression forest", cate_model = "OLS"))
})

test_that("x_learner returns error if propensity model invalid", {
  expect_error(x_learner(data = lalonde, y = "re78", w = "treat", base_learner = "OLS", cate_model = "OLS", propensity_model = "test"))
})

test_that("x_learner runs without errors for all propensity_models", {
  expect_silent(x_learner(data = lalonde, y = "re78", w = "treat", base_learner = "regression forest", cate_model = "regression forest", propensity_model = "logistic"))
  expect_silent(x_learner(data = lalonde, y = "re78", w = "treat", base_learner = "regression forest", cate_model = "regression forest", propensity_model = "lasso"))
  expect_silent(x_learner(data = lalonde, y = "re78", w = "treat", base_learner = "regression forest", cate_model = "regression forest", propensity_model = "causal forest"))
})