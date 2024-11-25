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
