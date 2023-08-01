#'
#' Small, artificially-generated toy data designed for testing the easytestresults package
#'
#' @docType data
#'
#' @usage data(toydata)
#'
#' @format An object of class \code{"data.frame"}
#' \describe{
#'  \item{condition}{sample experimental condition variable)}
#'  \item{outcome_nb}{A sample outcome variable of counts (integer format) that follows a negative binomial distribution
#'  \item{outcome_binom}{A sample outcome variable of 0s or 1s to simulate a binary outcome for logistic regression
#'  \item{outcome_norm}{A sample continuous, normally distributed outcome variable to simulate a DV for OLS regression
#'  \item{segment}{A categorical/factor variable to demonstrate a user segment}}
#' }
#' @references This data set was artificially created for the easytestresults package.
#' @keywords datasets
#' @examples
#'
#' data(toydata)
#' head(toydata)
#'
"toydata"

# Set the number of observations
n <- 500

# Simulate the continuous outcome variables
outcome_norm <- rnorm(n, mean = 0, sd = 1)
outcome_binom <- rbinom(n, size = 1, prob = .6)
outcome_nbinom <- rnbinom(n, mu = 6, size = .2)

# Define the three levels of the categorical predictor (condition)
condition_levels <- c("control", "high", "low")

# Simulate the categorical predictor (condition)
condition <- sample(condition_levels, n, replace = TRUE)

# Define the three levels of the factor predictor (segment)
segment_levels <- c("male", "female", "unknown")

# Simulate the factor predictor (segment)
segment <- sample(segment_levels, n, replace = TRUE)

# Create the data frame
toydata <- data.frame(outcome_norm, outcome_binom, outcome_nbinom, condition, segment)


save(toydata, file = "~/Documents/easytestresults//data/toydata.RData")
df <- toydata
colnames(df)
