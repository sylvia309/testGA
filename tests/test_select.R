# TESTS
# ------------------------------------------------------------------------------
# Organization: Each set of test corresponding to each fun or subfun called in
#               select() is presented in a separate context
# Description:  1. Tests for functions initialization, fit and select
#               2. Handles problematic cases inputed by user and
#               stochasic components of the algorithm

source("../R/selectFinal.R", chdir = TRUE)
library(testthat)

#------------------------------------------------------------------------------

context("Initialization")
# Test for the function of initialization of matrix

test_that("initialization has a non admitted input (NULL or NA)", {
  x <- NULL
  z <- NA

  expect_length(initialization(x), 0)
  expect_error(initialization(z), "incorrect number of dimensions")
})

#------------------------------------------------------------------------------

context("Fit")
# Test for the function of fit

# Test for checking dataset input is not NA
test_that("fit of has a non admitted dataset: NA", {
  z <- NA

  # Parameters
  f <- AIC         # Default: AIC for fitness function
  model <- lm      # Default: fitting linear model
  generation <- 20 # Default: generation size
  size <- 50       # Default: population size

  matrix <- initialization(mtcars,size)

  expect_error(fit(z, matrix, f, model), "incorrect number of dimensions")
})

# Test for checking function works with all types of dataset (including non numeric vars only)
test_that("fit is the correct type for a sample R dataset", {

  # Parameters
  f <- AIC         # Default: AIC for fitness function
  model <- lm      # Default: fitting linear model
  generation <- 20 # Default: generation size
  size <- 50       # Default: population size

  matrix <- initialization(mtcars, size)
  expect_type(fit(mtcars, matrix, f, model), "list")

  matrix <- initialization(iris, size)
  expect_type(fit(iris, matrix, f, model), "list")

  matrix <- initialization(ToothGrowth, size)
  expect_type(fit(ToothGrowth, matrix, f, model), "list")

  matrix <- initialization(USArrests, size)
  expect_type(fit(USArrests , matrix, f, model), "list")
})

#------------------------------------------------------------------------------

context("Selection")

test_that("parent matrix has the same size with the initial matrix",{
  test_matrix <- matrix(0, 5, 10)
  test_matrix <- apply(test_matrix, c(1,2), function(x) sample(c(TRUE, FALSE),1))
  test_fitness <- sample(10)
  expect_true(identical(dim(selection(test_fitness, test_matrix)), dim(test_matrix)))
})

test_that("parent matrix is different from the initial matrix",{
  test_matrix <- matrix(0, 12, 26)
  test_matrix <- apply(test_matrix, c(1,2), function(x) sample(c(TRUE, FALSE),1))
  test_fitness <- sample(26)
  expect_false(identical(selection(test_fitness, test_matrix), test_matrix))
})

#------------------------------------------------------------------------------

context("Crossover")

# Test that we get different results when calling crossover function
test_that("the result of crossover is random",{
  test_indiv1 <- sample(c(TRUE,FALSE), 40, TRUE)
  test_indiv2 <- sample(c(TRUE,FALSE), 40, TRUE)
  expect_false(identical(crossover(test_indiv1, test_indiv2), crossover(test_indiv1, test_indiv2)))
})

# Test that the new individuals have the same length as their parents
test_that("length of the individual does not change",{
  test_indiv1 <- sample(c(TRUE,FALSE), 12, TRUE)
  test_indiv2 <- sample(c(TRUE,FALSE), 12, TRUE)
  expect_true(length(crossover(test_indiv1, test_indiv2)['newIndiv1',]) == length(crossover(test_indiv1, test_indiv2)['newIndiv2',]))
  expect_true(length(crossover(test_indiv1, test_indiv2)['newIndiv1',]) == length(test_indiv1))
})

#------------------------------------------------------------------------------

context("Reproduction")

test_that("the crossed over matrix has the same size with parent matrix",{
  test_matrix <- matrix(0, 5, 10)
  test_matrix <- apply(test_matrix, c(1,2), function(x) sample(c(TRUE, FALSE),1))
  expect_true(identical(dim(reproduction(test_matrix)), dim(test_matrix)))

  test_matrix2 <- matrix(0, 16, 20)
  test_matrix2 <- apply(test_matrix2, c(1,2), function(x) sample(c(TRUE, FALSE),1))
  expect_true(identical(dim(reproduction(test_matrix2)), dim(test_matrix2)))
})


#------------------------------------------------------------------------------
context("Select")

# Test that the output is of class 'list'
test_that("select has output of incorrect class: list", {

  test_result <- select(mtcars, AIC, lm)
  expect_is(test_result, "list")
  expect_is(test_result$finalPredictors, "character")
  expect_is(test_result$finalCoeff, "numeric")

  test_result2 <- select(mtcars, AIC, glm)
  expect_is(test_result2, "list")
  expect_is(test_result2$finalPredictors, "character")
  expect_is(test_result2$finalCriterion, "numeric")
  expect_is(test_result2$finalCoeff, "numeric")

  test_result3 <- select(iris, AIC, glm)
  expect_is(test_result3, "list")
  expect_is(test_result3$finalPredictors, "character")
  expect_is(test_result3$finalCriterion, "numeric")
  expect_is(test_result3$finalCoeff, "numeric")

  test_result4 <- select(USArrests, AIC, lm)
  expect_is(test_result4, "list")
  expect_is(test_result4$finalPredictors, "character")
  expect_is(test_result4$finalCriterion, "numeric")
  expect_is(test_result4$finalCoeff, "numeric")

  test_result5 <- select(ToothGrowth, AIC, lm)
  expect_is(test_result5, "list")
  expect_is(test_result5$finalPredictors, "character")
  expect_is(test_result5$finalCriterion, "numeric")
  expect_is(test_result5$finalCoeff, "numeric")

})

# Test-version of initialization function: binomial
# Input:  initial dataset
# Output: (nrow,50) matrix with initial values for population

# We know that final criterion function is 154.2075 for this dataset

#Using sample dataset initialized by uniform
test_that("select output fitness fn converges to a final criterion in neighborhood 154.2075 +/-
          alpha* Machine Epsilon", {

expect_lte(
  #I create a 'mock' function 'initialization' to test desired feature
  with_mock(
    initialization = function(dataset, size){

      featuresLength <- length(dataset[,-1])  #

      # Create randomized matrix
      set.seed(1)
      individualsMatrix <- t(replicate(featuresLength,as.logical(round(runif(size, 0,1)))))

      return(individualsMatrix)
    },
    select(mtcars, AIC, lm, 20, 50)[2][[1]]
    ),
  154.2075+10^(20)*.Machine$double.eps
  )

expect_gte(
  #I create a 'mock' function 'initialization' to test desired feature
  with_mock(
    initialization = function(dataset, size){

      featuresLength <- length(dataset[,-1])  #

      # Create randomized matrix
      set.seed(1)
      individualsMatrix <- t(replicate(featuresLength,as.logical(round(runif(size, 0,1)))))

      return(individualsMatrix)
    },
    select(mtcars, AIC, lm, 20, 50)[2][[1]]
  ),
  154.2075-10^(20)*.Machine$double.eps
)
})


#-----------------------------
# The following two tests aim to create upper/lower bound criterium function with a known dataset (mtcars)

# Test-version of initialization function: binomial
# Input:  initial dataset
# Output: (nrow,50) matrix with initial values for population
# We know that final criterion function is 154.2075 in using the sample dataset
# and initialized by a uniform

test_that("select output converges to a final criterion in neghborhood 154.2075 +/- alpha* Machine Epsilon", {
    expect_lte(
      #I create a 'mock' function 'initialization' to test desired feature
      with_mock(
        initialization = function(dataset,size){

          featuresLength <- length(dataset[,-1])

          # Simulated initialized mtx: uniform
          set.seed(1)
          individualsMatrix <- t(replicate(featuresLength,as.logical(round(rbinom(size, 1,1/2)))))
          return(individualsMatrix)
        },
        select(mtcars, AIC, lm, 20, 50)[2][[1]]
      ),
      154.2075+10^(20)*.Machine$double.eps
    )

    expect_gte(
      #I create a 'mock' function 'initialization' to test desired feature
      with_mock(
        initialization = function(dataset,size){

          featuresLength <- length(dataset[,-1])

          # Simulated initialized mtx: uniform
          set.seed(1)
          individualsMatrix <- t(replicate(featuresLength,as.logical(round(rbinom(size, 1,1/2)))))
          return(individualsMatrix)
        },
        select(mtcars, AIC, lm, 20, 50)[2][[1]]
      ),
      154.2075-10^(20)*.Machine$double.eps
    )
})
