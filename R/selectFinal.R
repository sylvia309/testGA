

#Customizable parameters of out fitness function
f <- AIC         # Default: AIC for fitness function
model <- lm      # Default: fitting linear model
generation <- 20 # Default: generation size
size <- 50       # Default: population size



initialization <- function(dataset, size){

  featuresLength <- length(dataset[,-1])  #

  # Create randomized matrix
  individualsMatrix <- t(replicate(featuresLength,sample(x = c(FALSE,TRUE), size = size, replace = TRUE)))

  return(individualsMatrix)
}

# parameters: original data and the pop matrix
# returns a list of the selection error vector and the fitness vector of one generation and the fittest chromosome
# the length of the two vectors equals to the number of individuals/chromosomes

fit <- function(data, matrix, f, model){

  selectionError <- c()

  y <- data[,1]

  chrCount <- ncol(matrix)

  featuresLength <- nrow(matrix)

  #stopifnot((ncol(data)-1) == featuresLength)
  for(i in seq(1,chrCount)){
    if(sum(matrix[,i]) > 0){
      y <- data[,1]
      feature_index <- seq(2, featuresLength+1)
      selected_index <- subset(feature_index, matrix[,i])
      variables<- colnames(data)[selected_index]
      eq <- as.formula(paste("y", paste(variables, collapse = " + "),sep = " ~ "))
      selectionError <- c(selectionError, f(model(eq, data=data)))      # AIC Lin. model with no intercept
      modelRes <- model(eq, data=data)

    }
    else{

      selectionError <- c(selectionError, Inf)
    }
  }
  fitness <- 1.5 * rank(-selectionError) #


  #Return a list with tree components
  result <- list()
  result$error <- selectionError                 # Selection errors
  result$fitness <- fitness                      # Fitness vals
  result$fittest <- matrix[, which.max(fitness)] # Fit tests vals (T or F)
  return(result)

}
# parameter: fitness vector and the individual matrix
# returns a list; parent chromosomes - same length of chromosomes & fittest
selection <- function(fitness, matrix){
  chrCount <- dim(matrix)[2]

  stopifnot(!chrCount%%2)   #If not
  parent_index <- sample(x=1:chrCount, size=chrCount, replace=TRUE, prob=unlist(fitness))  #Unlisted prob. weights
  parent <- matrix[, parent_index]
  return(parent)
}

#temp2 <- selection(fitness, matrix)


# If number is lower than mutation rate then we change the variable
# According to literature the mutation rate is 1/(# of features)

mutation <- function(indiv){

  m <- length(indiv)
  doMutation <- sample(c(1,0),length(indiv),prob=c(1/m,1-(1/m)),replace=T)

  for(i in 1:length(doMutation)) {

    if (doMutation[i] == 1 && indiv[i] == 0) {
      indiv[i] = TRUE
    } else if (doMutation[i] == 1 && indiv[i] == 1) {
      indiv[i] = FALSE
    } else {

    }
  }

  return(indiv)
}


# This crossover function is going to take in 2 individuals
# This does a uniform crossover
# Also known as the "coin flipping" method in which each gene has a 50% chance of coming from individual 1 and a 50% chance of coming from individual 2

crossover <- function(indiv1, indiv2){

  newIndiv1 <- rep(FALSE,length(indiv1))
  newIndiv2 <- rep(FALSE,length(indiv1))

  for(i in 1:length(indiv1)){

    newIndiv1[i] <- sample(x = c(indiv1[i],indiv2[i]),1)
    newIndiv2[i] <- sample(x = c(indiv1[i],indiv2[i]),1)

  }

  # Call the mutate function on crossed individuals
  newIndiv1 <- mutation(newIndiv1)
  newIndiv2 <- mutation(newIndiv2)

  return(rbind(newIndiv1,newIndiv2))
}

# This takes in the parent chromosomes matrix that we got from the selection function
# Note this has the same dimension as the individuals matri
reproduction <- function(parentMatrix) {

  featuresLength <- nrow(parentMatrix)
  size <- ncol(parentMatrix)
  # Initialize blank matrix
  crossedMatrix <- matrix(FALSE,featuresLength, size)

  for (i in seq(1, size, 2)) {

    intermediateCross <- crossover(parentMatrix[,i],parentMatrix[,i+1])
    crossedMatrix[,i] <- intermediateCross[1,]
    crossedMatrix[,i+1] <- intermediateCross[2,]

  }

  return(crossedMatrix)
}

#' Implements a Genetic Algorithm for Variable Selection
#'
#' @description This R package implements a genetic algorithm for variable selection in regression problems. The functionality of this package extends to linear regressions or generalized linear models.
#'
#' \code{select} returns a list that includes the best fitted predictors, fitness score, and coefficient values for the predictors.
#' @param dataset dataset that contains the predictors and the response variables. Note: The first column must be the response variable.
#' @param criterion function that specifies which objective criterion/fitness function will be used to calculate fitness scores. Note: The default is AIC
#' @param model model that user wishes to fit to the data either using \code{lm} or \code{glm}. Note: The default is \code{lm}
#' @param generation total number of generations that the user wants through iterate through. Note: The default is 20 generations
#' @param size population size Note: The number entered must be even. The default is set to 50 chromosomes
#'
#' @details The goal of the package is to utilize the following components of genetic algortihms:
#' selection, crossover and mutation in order to improve feature selection.
#' This type of algorithm mimics the process of natural selection, hence it favors the fittest individuals over all the generations.
#'
#' This package utilizes crossover and mutation processes adapted to a problem of variable selection. In this package, the population is represented by a matrix and the length (of chromosomes in traditional genetic algorithms) is equivalent to the number of columns of this population matrix. Preceeding the three building blocks we have mentioned, the algorithm requires two stages: initialization and fitness assignment.
#' @return \code{select} returns a list the best fitted predictors, fitness score, and coefficient values for the predictors.
#' @export
#'
#' @usage select(dataset, criterion, model, generation, size)
#' @examples select(dataset = mtcars, criterion = AIC, model = lm, generation = 50, size = 10)
#' @examples select(dataset = iris, criterion = BIC, model = glm)
#' @examples select(iris)
#'
select <- function(dataset, criterion = AIC, model = lm, generation = 20, size = 50){


  individualMatrix <- initialization(dataset, size)

  for(i in 1:generation) {

    fitnessFuncResults <- fit(dataset, individualMatrix, criterion, model)
    parentChromMatrix <- selection(fitnessFuncResults$fitness, individualMatrix)
    individualMatrix <- reproduction(parentChromMatrix)


  }

  finalPredictors <- colnames(dataset)[-1][c(fitnessFuncResults$fittest)]
  finalCriterion <- fitnessFuncResults$error[which.max(fitnessFuncResults$fitness)]

  y <- dataset[,1]
  feature_index <- seq(2, ncol(dataset))
  selected_index <- subset(feature_index, fitnessFuncResults$fittest)
  variables<- colnames(dataset)[selected_index]
  eq <- as.formula(paste("y", paste(variables, collapse = " + "),sep = " ~ "))
  finalCoeff <- model(eq, data=dataset)$coefficients

  finalModel <- list()
  finalModel$finalPredictors <- finalPredictors
  finalModel$finalCriterion <- finalCriterion
  finalModel$finalCoeff <- finalCoeff

  return(finalModel)
}
