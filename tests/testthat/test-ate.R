context("ate")
source('../../R/ate.R')
source('../../R/propensity.R')
load("../../data/lalonde.RData")

test_that("naive_ate runs with and without cf", {
  expect_silent(naive_ate(data = lalonde, y = "re78", w = "treat", cf = TRUE))
  expect_silent(naive_ate(data = lalonde, y = "re78", w = "treat", cf = FALSE))
})

test_that("ipw_ate check that error is thrown when propensity and data dimensions don't match", {
  expect_error(imw_ate(data = lalonde, y = "re78", w = "treat", p = seq(1,nrow(lalonde)+5)))
})

test_that("ipw_ate runs with and without cf", {
  p <- propensity_score(data = lalonde, y = "re78", w = "treat", model = "logistic")
  expect_silent(ipw_ate(data = lalonde, y = "re78", w = "treat", p = p, cf = TRUE))
  expect_silent(ipw_ate(data = lalonde, y = "re78", w = "treat", p = p, cf = FALSE))
})

test_that("prop_weighted_ols_ate check that error is thrown when propensity and data dimensions don't match", {
  expect_error(prop_weighted_ols_ate(data = lalonde, y = "re78", w = "treat", p = seq(1,nrow(lalonde)+5)))
})

test_that("prop_weighted_ols_ate runs with and without cf", {
  p <- propensity_score(data = lalonde, y = "re78", w = "treat", model = "logistic")
  expect_silent(prop_weighted_ols_ate(data = lalonde, y = "re78", w = "treat", p = p, cf = TRUE))
  expect_silent(prop_weighted_ols_ate(data = lalonde, y = "re78", w = "treat", p = p, cf = FALSE))
})

test_that("aipw_ate check that error is thrown when propensity and data dimensions don't match", {
  expect_error(aipw_ate(data = lalonde, y = "re78", w = "treat", p = seq(1,nrow(lalonde)+5)))
})

test_that("aipw_ate runs with and without cf", {
  p <- propensity_score(data = lalonde, y = "re78", w = "treat", model = "logistic")
  expect_silent(aipw_ate(data = lalonde, y = "re78", w = "treat", p = p, cf = TRUE))
  expect_silent(aipw_ate(data = lalonde, y = "re78", w = "treat", p = p, cf = FALSE))
})

test_that("aipw_ate runs with and without x specified", {
  p <- propensity_score(data = lalonde, y = "re78", w = "treat", model = "logistic")
  expect_silent(aipw_ate(data = lalonde, y = "re78", w = "treat", p = p, cf = TRUE))
  expect_silent(aipw_ate(data = lalonde, x = c("age"), y = "re78", w = "treat", p = p, cf = FALSE))
})

test_that("propensity_strat_ate check that error is thrown when propensity and data dimensions don't match", {
  expect_error(propensity_strat_ate(data = lalonde, y, w, p=seq(1,nrow(lalonde)+5), n_strata = 4, min_obs = 1, model = "naive"))
})

test_that("propensity_strat_ate check that error is thrown when ate model is not valid", {
  p <- propensity_score(lalonde, y = "re78", w = "treat", model = "logistic")
  expect_error(propensity_strat_ate(data = lalonde, y, w, p=p, n_strata = 4, min_obs = 1, model = "test"))
})

test_that("propensity_strat_ate check error is thrown if min_obs and/or n_strata not valid integers", {
  p <- propensity_score(lalonde, y = "re78", w = "treat", model = "logistic")
  expect_error(propensity_strat_ate(data = lalonde, y, w, p=p, n_strata = 4.5, min_obs = 1, model = "naive"))
  expect_error(propensity_strat_ate(data = lalonde, y, w, p=p, n_strata = 4, min_obs = 1.5, model = "naive"))
  expect_error(propensity_strat_ate(data = lalonde, y, w, p=p, n_strata = 4.5, min_obs = 1.5, model = "naive"))
  expect_error(propensity_strat_ate(data = lalonde, y, w, p=p, n_strata = -2, min_obs = 1, model = "naive"))
  expect_error(propensity_strat_ate(data = lalonde, y, w, p=p, n_strata = 4, min_obs = -2, model = "naive"))
  expect_error(propensity_strat_ate(data = lalonde, y, w, p=p, n_strata = -4, min_obs = -2, model = "naive"))
  expect_error(propensity_strat_ate(data = lalonde, y, w, p=p, n_strata = 0, min_obs = 0, model = "naive"))
})

test_that("propensity_strat_ate check error thrown if n_strata larger than data size", {
  p <- propensity_score(lalonde, y = "re78", w = "treat", model = "logistic")
  expect_error(propensity_strat_ate(data = lalonde, y, w, p=p, n_strata = nrow(lalonde)+5, min_obs = 1, model = "naive"))
})

test_that("propensity_strat_ate check error thrown if min_obs larger than data size", {
  p <- propensity_score(lalonde, y = "re78", w = "treat", model = "logistic")
  expect_error(propensity_strat_ate(data = lalonde, y, w, p=p, n_strata = 4, min_obs = nrow(lalonde)+5, model = "naive"))
})

test_that("propensity_strat_ate check that n_strata bounding works properly", {
  p <- propensity_score(lalonde, y = "re78", w = "treat", model = "logistic")
  output_1 = propensity_strat_ate(data = lalonde, y = "re78", w = "treat", p=p, n_strata = 11, min_obs = 40, model = "naive")
  output_2 = propensity_strat_ate(data = lalonde, y = "re78", w = "treat", p=p, n_strata = 20, min_obs = 40, model = "naive")
  expect_equal(output_1, output_2)
})

test_that("propensity_strat_ate works for different ate models", {
  p <- propensity_score(lalonde, y = "re78", w = "treat", model = "logistic")
  expect_silent(propensity_strat_ate(data = lalonde, y = "re78", w = "treat", p=p, n_strata = 4, min_obs = 1, model = "naive"))
  expect_silent(propensity_strat_ate(data = lalonde, y = "re78", w = "treat", p=p, n_strata = 4, min_obs = 1, model = "ipw"))
  expect_silent(propensity_strat_ate(data = lalonde, y = "re78", w = "treat", p=p, n_strata = 4, min_obs = 1, model = "weighted ols"))
  expect_silent(propensity_strat_ate(data = lalonde, y = "re78", w = "treat", p=p, n_strata = 4, min_obs = 1, model = "aipw"))
})