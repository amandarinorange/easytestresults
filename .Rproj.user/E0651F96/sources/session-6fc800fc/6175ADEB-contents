context("check-easyresults")
library(testthat)        # load testthat package
library(lmtest)
library(tibble)
library(stats)
library(emmeans)
library(interactions)
library(dplyr)
library(tidyr)
library(purrr)
library(easytestresults)       # load our package

# Test whether the output is a list
test_that("easyresults() returns a list", {
  output <- easy_results(outcome = "outcome_norm", treatment_name = "condition", df = toydata, family = "lm")
  expect_type(output, "list")
})


colnames(toydata)
