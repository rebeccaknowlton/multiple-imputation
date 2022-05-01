library(tidyverse)
library(gridExtra)
library(locfit)
library(mice)
library(mitools)
library(quantreg)
library(rpart)
library(caret)
library(rpart.plot)
library(randomForest)

# generate data with n observations
# X1 fully observed, generated from normal distribution
# X2 regressed on X1 with optional quadratic term
# X3 regressed on X1 and X2 with optional quadratic terms and interaction
generate_data <- function(n = 1000, 
                          X1_mean = 0, X1_sd = 4, 
                          X2_noise_mean = 1, X2_noise_sd = 2,
                          X3_noise_mean = 0, X3_noise_sd = 3, 
                          linear_coefs = c(1,2,-1.5), 
                          quadratic_coefs = c(0,0,0), 
                          interaction_coef = 0) {
  X1 <- rnorm(n, X1_mean, X1_sd)
  X2 <- linear_coefs[1] * X1 + quadratic_coefs[1] * X1 ^ 2 +
    rnorm(n, X2_noise_mean, X2_noise_sd)
  X3 <- linear_coefs[2] * X1 + quadratic_coefs[2] * X1 ^ 2 +
    linear_coefs[3] * X2 + quadratic_coefs[3] * X2 ^ 2 + 
    interaction_coef * X1 * X2 + rnorm(n, X3_noise_mean, X3_noise_sd)
  df <- data.frame(X1 = X1,
                   X2 = X2,
                   X3 = X3)
  return(df)
}

# take input dataframe and determine missingness for X2 and X3
induce_missingness <- function(df, X2_constants = c(-0.25, 5, 2),
                               X3_constants = c(-0.5, 2, 3)) {
  # probabilities that X2 and X3 are missing
  X2_prob_missing <- expit(X2_constants[1] - (df$X1 / X2_constants[2]) ^ X2_constants[3])
  X3_prob_missing <- expit(X3_constants[1] + ((df$X1 - X3_constants[2]) / X3_constants[3]))
  
  # indicators for missingness
  R2 <- rbernoulli(nrow(df), X2_prob_missing)
  R3 <- rbernoulli(nrow(df), X3_prob_missing)
  
  # versions of X2 and X3 with blanked out missing values
  X2_missing <- df$X2
  X2_missing[which(R2 == 1)] <- NA
  X3_missing <- df$X3
  X3_missing[which(R3 == 1)] <- NA
  
  updated_df <- data.frame(X1 = df$X1,
                           X2_fully_obs = df$X2,
                           X3_fully_obs = df$X3,
                           X2 = X2_missing,
                           X3 = X3_missing,
                           R2 = R2,
                           R3 = R3)
  return(updated_df)
}

# plot data and report what proportion of values are missing
summarize <- function(df) {
  n <- nrow(df)
  
  # summarize how much data is missing
  cat((sum(df$R2) / n) * 100, "% of X2 data missing \n", sep = "")
  cat((sum(df$R3) / n) * 100, "% of X3 data missing \n", sep = "")
  cat((length(which(df$R2 == 1 | df$R3 == 1)) / n) * 100, "% of cases have either X2 or X3 missing", sep = "")
  # plots
  plot1 <- ggplot(data = df, aes(x = X1, y = X2_fully_obs, color = factor(R2))) + 
    geom_point() + theme_bw() +
    scale_color_manual(values = c("FALSE" = "gray",
                                  "TRUE" = "orange"),
                       name = "X2 missing?")
  plot2 <- ggplot(data = df, aes(x = X1, y = X3_fully_obs, color = factor(R3))) + 
    geom_point() + theme_bw() +
    scale_color_manual(values = c("FALSE" = "gray",
                                  "TRUE" = "orange"),
                       name = "X3 missing?")
  grid.arrange(plot1, plot2, ncol=2)
}

# plot imputed data
plot_imputed <- function(imp, df, k = 1) {
  estimates <- complete(imp, k)
  plot1 <- ggplot(data = df, aes(x = X1, y = X2_fully_obs, color = factor(R2))) + 
    geom_point() + theme_bw() +
    scale_color_manual(values = c("FALSE" = "gray",
                                  "TRUE" = "orange"),
                       name = "X2 missing?") +
    theme(legend.position = "none")
  plot2 <- ggplot(data = df, aes(x = X1, y = X3_fully_obs, color = factor(R3))) + 
    geom_point() + theme_bw() +
    scale_color_manual(values = c("FALSE" = "gray",
                                  "TRUE" = "orange"),
                       name = "X3 missing?") + 
    theme(legend.position = "none")
  plot3 <- ggplot(data = estimates, aes(x = estimates[,1], y = estimates[,2])) +
    geom_point(color = "gray") + theme_bw()
  plot4 <- ggplot(data = estimates, aes(x = estimates[,1], y = estimates[,3])) +
    geom_point(color = "gray") + theme_bw()
  grid.arrange(plot1, plot2, plot3, plot4, ncol = 2)
}

# generate some data
my_data <- induce_missingness(generate_data(linear_coefs = c(0.5, -0.5, 0.2), 
                                            quadratic_coefs = c(-0.4, -0.3, -1), 
                                            interaction_coef = -0.2))
summarize(my_data)


# write function to duplicate and concatenate the data with missingness
duplicate_concatenate <- function(df) {
  n <- nrow(df)
  df2 <- df # copy the original df
  df2$X2 <- rep(NA, n) # blank out all X2 values
  df2$R2 <- rep(TRUE, n)
  df2$X3 <- rep(NA, n) # blank out all X3 values
  df2$R3 <- rep(TRUE, n)
  rbind(df, df2)
}

my_data_dup <- duplicate_concatenate(my_data)

my_data_input <- data.frame(my_data_dup$X1, 
                            my_data_dup$X2,
                            my_data_dup$X3)
imp_my_data <- mice(my_data_input, method = "norm", m = 100, maxit = 5)
plot_imputed(imp_my_data, my_data)

# get the imputed datasets for the original data and the duplicated version
get_imp_datasets_orig <- function(imp) {
  m <- imp$m
  n <- nrow(imp$data)
  imp_datasets <- vector(mode = "list", length = m)
  for (i in 1:m) {
    imp_datasets[[i]] <- complete(imp, i)[1:(n/2),]
    colnames(imp_datasets[[i]]) <- c("X1", "X2", "X3")
  }
  return(imp_datasets)
}
get_imp_datasets_dup <- function(imp) {
  m <- imp$m
  n <- nrow(imp$data)
  imp_datasets <- vector(mode = "list", length = m)
  for (i in 1:m) {
    imp_datasets[[i]] <- complete(imp, i)[(n/2 + 1):n,]
    colnames(imp_datasets[[i]]) <- c("X1", "X2", "X3")
  }
  return(imp_datasets)
}

imp_data_orig <- get_imp_datasets_orig(imp_my_data)
imp_data_dup <- get_imp_datasets_dup(imp_my_data)

# posterior predictive p value
post_pred_p <- function(inputs_df1, inputs_df2) {
  m <- length(inputs_df1$estimates)
  num_estimates <- length(inputs_df1$estimates[[1]])
  p <- rep(NA, num_estimates)
  for (i in 1:num_estimates) {
    tmp_idx <- seq(from = i, to = m * num_estimates, by = num_estimates)
    p[i] <- sum(unlist(inputs_df1$estimates)[tmp_idx] > unlist(inputs_df2$estimates)[tmp_idx]) / m
  }
  return(p)
}

mean_inputs <- function(imp_datasets, col_num) {
  estimates <- lapply(imp_datasets, function(x) {mean(x[,col_num])})
  variances <- lapply(imp_datasets, function(x) {var(x[,col_num]) / nrow(x)})
  df <- nrow(imp_datasets[[1]]) - 1
  inputs <- list(estimates = estimates, 
                 variances = variances,
                 df = df)
  return(inputs)
}

regression_inputs <- function(imp_datasets, formula) {
  estimates <- lapply(imp_datasets, function(x) {coef(lm(formula, data = x))})
  variances <- lapply(imp_datasets, function(x) {vcov(lm(formula, data = x))})
  df <- summary(lm(formula, imp_datasets[[1]]))$df[2]
  inputs <- list(estimates = estimates, 
                 variances = variances,
                 df = df)
  return(inputs)
}


# X2 means
post_pred_p(mean_inputs(imp_data_orig, 2), mean_inputs(imp_data_dup, 2))

# X3 means
post_pred_p(mean_inputs(imp_data_orig, 3), mean_inputs(imp_data_dup, 3))

# regress X2 on X1
post_pred_p(regression_inputs(imp_data_orig, X2 ~ X1), regression_inputs(imp_data_dup, X2 ~ X1))

# regress X3 on X1
post_pred_p(regression_inputs(imp_data_orig, X3 ~ X1), regression_inputs(imp_data_dup, X3 ~ X1))

# regress X2 on X1^2, no intercept
post_pred_p(regression_inputs(imp_data_orig, X2 ~ -1 + I(X1^2)), regression_inputs(imp_data_dup, X2 ~ -1 + I(X1^2)))

# with intercept
post_pred_p(regression_inputs(imp_data_orig, X2 ~ I(X1^2)), regression_inputs(imp_data_dup, X2 ~ I(X1^2)))


### Build classifier to predict orig or dup ###

add_indicator_col <- function(imp_datasets) {
  m <- length(imp_datasets)
  n <- nrow(imp_datasets[[1]])
  ind_col <- append(rep(1, n / 2), rep(0, n / 2))
  for (i in 1:m) {
    imp_datasets[[i]]$ind <- ind_col
  }
  return(imp_datasets)
}

add_R2_R3 <- function(imp_datasets, df_dup) {
  m <- length(imp_datasets)
  n <- nrow(imp_datasets[[1]])
  for (i in 1:m) {
    imp_datasets[[i]]$R2 <- df_dup$R2
    imp_datasets[[i]]$R3 <- df_dup$R3
  }
  return(imp_datasets)
}

get_imp_datasets <- function(imp) {
  m <- imp$m 
  imp_datasets <- vector(mode = "list", length = m)
  for (i in 1:m) {
    imp_datasets[[i]] <- complete(imp, i)
    colnames(imp_datasets[[i]]) <- c("X1", "X2", "X3")
  }
  return(imp_datasets)
}

# these are the m imputed datasets with an ind col to predict  
classifier_data <- add_indicator_col(get_imp_datasets(imp_my_data))
classifier_data <- add_R2_R3(classifier_data, my_data_dup)
stacked_classifier_data <- do.call(rbind, classifier_data)

# example: try cart on one of them
k <- 1 # which imputed dataset to use
model <- rpart(ind ~ X1 + X2 + X3, data = classifier_data[[k]])
plot(model)
text(model)
rpart.plot(model)
plot_imputed(imp_my_data, my_data, k)

# cart tree on stacked data
model2 <- rpart(ind ~ X1 + X2 + X3, data = stacked_classifier_data)
plot(model2)
text(model2)
rpart.plot(model2)

# plot stacked data 
plotX2 <- ggplot(data = stacked_classifier_data, aes(x = X1, y = X2, color = factor(ind))) +
  geom_point(alpha = 0.2) + theme_bw() +
  scale_color_manual(values = c("0" = "gray",
                                "1" = "orange"),
                     name = "Duplicated data?") +
  theme(legend.position = "none")

plotX3 <- ggplot(data = stacked_classifier_data, aes(x = X1, y = X3, color = factor(ind))) +
  geom_point(alpha = 0.2) + theme_bw() + 
  scale_color_manual(values = c("0" = "gray",
                                "1" = "orange"),
                     name = "Duplicated data?") +
  theme(legend.position = "none")
grid.arrange(plotX2, plotX3, ncol = 2)

# write function to plot data colored by classifier output
plot_preds <- function(df, model, prob = FALSE) {
  if (prob) {
    df$pred <- model %>% predict(df, type = "prob")
    df$pred <- df$pred[,2] #take second col because that's the prob = 1
  } else {
    df$pred <- model %>% predict(df)
  }
  plot1 <- ggplot(data = df, aes(x = X1, y = X2)) +
    geom_point(alpha = 0.2, aes(color = pred)) + theme_bw() + 
    scale_color_gradient(low = "gray", high = "dodgerblue4")
  plot2 <- ggplot(data = df, aes(x = X1, y = X3)) + 
    geom_point(alpha = 0.2, aes(color = pred)) + theme_bw() + 
    scale_color_gradient(low = "gray", high = "dodgerblue4")
  grid.arrange(plot1, plot2, ncol = 2)
}

plot_preds(stacked_classifier_data, model2)

# to force tree to be simpler, use only X1 and X2
model3 <- rpart(ind ~ X1 + X2, data = stacked_classifier_data)
rpart.plot(model3)

plot_preds(stacked_classifier_data, model3)

# also try X1 and X3
model4 <- rpart(ind ~ X1 + X3, data = stacked_classifier_data)
rpart.plot(model4)
plot_preds(stacked_classifier_data, model4)

# also try X2 and X3
model5 <- rpart(ind ~ X2 + X3, data = stacked_classifier_data)
rpart.plot(model5)
plot_preds(stacked_classifier_data, model5)

## Logistic Regression (instead of trees)

logistic_model <- glm(ind ~ X1 + X2 + X3, data = stacked_classifier_data, family = "binomial")
summary(logistic_model)

plot_preds(stacked_classifier_data, logistic_model)


## Random Forest
rf_model <- randomForest(as.factor(ind) ~ X1 + X2 + X3, data = stacked_classifier_data)
plot_preds(stacked_classifier_data, rf_model, prob = TRUE)
varImpPlot(rf_model)

## cart column-wise prediction

# predict R3 given X1 and X2
R3_cart_model <- rpart(R3 ~ X1 + X2, data = stacked_classifier_data)
rpart.plot(R3_cart_model)
plot_preds(stacked_classifier_data, R3_cart_model)

# plot R3 in stacked data
plotX2 <- ggplot(data = stacked_classifier_data, aes(x = X1, y = X2, color = factor(R3))) +
  geom_point(alpha = 0.2) + theme_bw() +
  scale_color_manual(values = c("FALSE" = "gray",
                                "TRUE" = "orange"),
                     name = "Duplicated data?") +
  theme(legend.position = "none")

plotX3 <- ggplot(data = stacked_classifier_data, aes(x = X1, y = X3, color = factor(R3))) +
  geom_point(alpha = 0.2) + theme_bw() + 
  scale_color_manual(values = c("FALSE" = "gray",
                                "TRUE" = "orange"),
                     name = "Duplicated data?") +
  theme(legend.position = "none")
grid.arrange(plotX2, plotX3, ncol = 2)

# predict R2 given X1 and X2
R2_cart_model <- rpart(R2 ~ X1 + X2, data = stacked_classifier_data)
rpart.plot(R2_cart_model)
plot_preds(stacked_classifier_data, R2_cart_model)

# plot R2 in stacked data
plotX2 <- ggplot(data = stacked_classifier_data, aes(x = X1, y = X2, color = factor(R2))) +
  geom_point(alpha = 0.2) + theme_bw() +
  scale_color_manual(values = c("FALSE" = "gray",
                                "TRUE" = "orange"),
                     name = "Duplicated data?") +
  theme(legend.position = "none")

plotX3 <- ggplot(data = stacked_classifier_data, aes(x = X1, y = X3, color = factor(R2))) +
  geom_point(alpha = 0.2) + theme_bw() + 
  scale_color_manual(values = c("FALSE" = "gray",
                                "TRUE" = "orange"),
                     name = "Duplicated data?") +
  theme(legend.position = "none")
grid.arrange(plotX2, plotX3, ncol = 2)
