# 00_multivariate_offsets.R
# example code for using offsets in the analysis of multivariate abundance data
library(mvabund)

# Example data
response <- mvabund(cbind(species1 = rpois(20, 5), 
                          species2 = rpois(20, 10)))  # Response matrix
predictors <- data.frame(Year = factor(rep(1:2, each = 10)), 
                         Site = rep(1:10, times = 2))
offsets <- log(rep(c(1, 2), each = 10))  # Example offset for sampling effort

# Fit the manyglm model
model <- manyglm(response ~ Year + offset(offsets), 
                 family = "poisson", 
                 data = predictors)

summary(model)
###########

# Example Setup
# 1. Simulated Data:
#   Stations: 5 total.
# Station A, B, C: 1 sample each.
# Station D, E: 3 replicates each.
# Species Abundances: Simulated for two species at each sample.
# Offset: The log of the number of samples (to account for differences in sampling effort).
# 2. Code Implementation:

library(mvabund)

# Step 1: Simulate species abundance data
set.seed(123)
response <- mvabund(cbind(
  species1 = c(rpois(1, 5), rpois(1, 6), rpois(1, 7), rpois(3, 10), rpois(3, 12)),
  species2 = c(rpois(1, 3), rpois(1, 4), rpois(1, 2), rpois(3, 8), rpois(3, 9))
))

# Step 2: Create metadata
stations <- factor(rep(c("A", "B", "C", "D", "E"), times = c(1, 1, 1, 3, 3)))
effort <- rep(c(1, 1, 1, 3, 3), times = c(1, 1, 1, 3, 3))  # Number of samples per station
offsets <- log(effort)  # Log-transformed effort

# Step 3: Fit the manyglm model
model <- manyglm(response ~ stations + offset(offsets), family = "poisson")

# Step 4: Summary of the model
summary(model);plot(model)

##################
# Load necessary libraries
library(mvabund)

# Simulate data
set.seed(123)
data <- data.frame(
  station = rep(1:3, each = 6),  # 3 stations
  year = rep(2010:2011, each = 3, times = 3),  # 2 years
  replicate = rep(1:3, times = 6),  # 3 replicates
  species1 = rnbinom(18, mu = 5, size = 1),
  species2 = rnbinom(18, mu = 3, size = 1)
)

# Calculate offset as log of replicate counts
data$offset <- log(3)  # Log of number of replicates, uniform here for simplicity

# Convert species data to mvabund format
abund_data <- mvabund(data[, c("species1", "species2")])

# Fit the manyglm model with an offset
fit <- manyglm(abund_data ~ data$year + data$station, offset = data$offset, family = "negative.binomial")

# Summary of results
summary(fit)
