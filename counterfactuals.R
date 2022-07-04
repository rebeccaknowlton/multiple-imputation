library(pdp)
library(randomForest)
library(tidyverse)
library(optimr)

data(pima)
pima$ind <- case_when(pima$diabetes == "pos" ~ 1,
                      pima$diabetes == "neg" ~ 0)

pima_data <- pima[c("pregnant", "glucose", "pressure", "mass", "pedigree", "age", "ind")]

pima_train <- pima_data[1:700,]
pima_test <- pima_data[701:768,]

# random forest to predict diabetes based on pregnant, glucose, pressure, mass, pedigree, and age
rf_model <- randomForest(as.factor(ind) ~ ., 
                         data = pima_train, na.action = na.roughfix)
varImpPlot(rf_model)


# this individual is predicted negative
obs <- pima_test[1,]
predict(rf_model, obs)
obs$ind[1]

# what could we change to get positive prediction?

# distance function
d <- function(x1, x2) {
  sum <- 0
  for(i in 1:(length(x1) - 1)) {
    sum <- sum + (abs(x1[i] - x2[i]) / mad(pima_train[,i], na.rm = TRUE))
  }  
  return(sum)
}

min_max_eq <- function(lam, model, x1, x2, y2) {
  pred <- as.numeric(predict(model, x2, type = "response")) - 1
  return(lam * ((pred - y2) ^ 2) + d(x1, x2))
}

# try just changing glucose value
resolution <- 100
range_vals <- range(pima_train$glucose, na.rm = TRUE)
lam_vals <- seq(from = 0, to = 2, length.out = resolution) 
glucose_vals <- seq(from = range_vals[1], to = range_vals[2], length.out = resolution)


obs_new <- obs
solution_mat <- matrix(nrow = resolution, ncol = 3)
solution_mat[,1] <- lam_vals
for (i in 1:resolution) {
  curr_lam <- lam_vals[i]
  gluc_min <- 1000
  gluc_min_val <- 1000
  for (j in 1:resolution) {
    curr_glucose <- glucose_vals[j] 
    obs_new$glucose <- curr_glucose
    val <- min_max_eq(curr_lam, rf_model, obs, obs_new, 1)
    if (val < gluc_min_val) {
      gluc_min <- curr_glucose
      gluc_min_val <- val
    }
  } 
  solution_mat[i,2] <- gluc_min
  solution_mat[i,3] <- as.numeric(gluc_min_val)
}
colnames(solution_mat) <- c("lam", "glucose", "val")

solution_mat

# we can see that it converges and that when glucose is increased 
# from 122.28 to 155.16, the individual would instead get a positive prediction
predict(rf_model, obs_new)


# now implement multivariate version 

optim_fn <- function(x) {
  min_max_eq()
}

solution_mat_multivar <- matrix(nrow = resolution, ncol = 9)
solution_mat_multivar[,1] <- lam_vals
for(i in 1:resolution) {
  curr_lam <- lam_vals[i]
  optim_fn <- function(x) {
    obs_new <- obs 
    obs_new$pregnant <- x[1]
    obs_new$glucose <- x[2]
    obs_new$pressure <- x[3]
    obs_new$mass <- x[4]
    obs_new$pedigree <- x[5]
    obs_new$age <- x[6]
    return(as.numeric(min_max_eq(curr_lam, rf_model, obs, obs_new, 1)))
  }
  sol <- optimr(obs[1:6], optim_fn, method = 'Nelder-Mead')
  solution_mat_multivar[i, 2:7] <- sol$par
  solution_mat_multivar[i, 8] <- sol$value
  solution_mat_multivar[i, 9] <- as.numeric(predict(rf_model, data.frame(t(sol$par)))) - 1
}
colnames(solution_mat_multivar) <- c("lam", "preg", "gluc", "press", "mass", "pedig", "age", "min max val", "pred")
solution_mat_multivar
# the convergence is messier in the multivariate setting and has multiple possible 
# combinations that can flip the prediction, but we see that increasing the age 
# and possibly increasing the glucose slightly result in a positive prediction
