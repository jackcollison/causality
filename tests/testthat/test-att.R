context("att")
library(testthat)
source('../../R/att.R')
source('../../R/propensity.R')
load("../../data/lalonde.RData")

test_that("aipw_rf_att check if error is thrown when dimension of data and propensity scores don't match", {
  expect_error(aipw_rf_att(data = lalonde, y = "re78", w = "treat", p = seq(1,nrow(lalonde)+5), num.trees = 100, mtry = 3))
})

test_that("rf_att runs with and without cf", {
  expect_silent(rf_att(data = lalonde, y = "re78", w = "treat", cf = TRUE, num.trees = 100, mtry = 3))
  expect_silent(rf_att(data = lalonde, y = "re78", w = "treat", cf = FALSE, num.trees = 100, mtry = 3))
})


test_that("aipw_rf_att runs with and without cf", {
  p <- propensity_score(lalonde, y = "re78", w = "treat", model = "logistic")
  expect_silent(aipw_rf_att(data = lalonde, y = "re78", w = "treat", p = p, cf = TRUE, num.trees = 100, mtry = 3))
  expect_silent(aipw_rf_att(data = lalonde, y = "re78", w = "treat", p = p, cf = FALSE, num.trees = 100, mtry = 3))
})