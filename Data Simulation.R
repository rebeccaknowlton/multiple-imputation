library(tidyverse)
library(gridExtra)
library(locfit)
library(mice)
library(mitools)
library(quantreg)

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

### true DGP: linear ###
data_linear <- induce_missingness(generate_data())
summarize(data_linear)
data_linear_input <- data.frame(data_linear$X1,
                                data_linear$X2,
                                data_linear$X3)
# imputation method = "norm"
imp_linear_linear <- mice(data_linear_input, method = "norm", m = 5, maxit = 5)
plot_imputed(imp_linear_linear, data_linear)
# imputation method = "norm" with quadratic term
imp_linear_quadratic <- mice(data_linear_input, method = "norm", 
                             formulas = list(data_linear.X1 = data_linear.X1 ~ data_linear.X2 + data_linear.X3 + I(data_linear.X2**2) + I(data_linear.X3**2),
                                             data_linear.X2 = data_linear.X2 ~ data_linear.X1 + data_linear.X3 + I(data_linear.X1**2) + I(data_linear.X3**2),
                                             data_linear.X3 = data_linear.X3 ~ data_linear.X1 + data_linear.X2 + I(data_linear.X1**2) + I(data_linear.X2**2)),
                             m = 5, maxit = 5)
plot_imputed(imp_linear_quadratic, data_linear)
# imputation method = "cart"
imp_linear_cart <- mice(data_linear_input, method = "cart", m = 5, maxit = 5)
plot_imputed(imp_linear_cart, data_linear)

### true DGP: interaction ###
data_interaction <- induce_missingness(generate_data(interaction_coef = 0.8))
summarize(data_interaction)
data_interaction_input <- data.frame(data_interaction$X1,
                                     data_interaction$X2,
                                     data_interaction$X3)
# imputation method = "norm"
imp_interaction_linear <- mice(data_interaction_input, method = "norm", m = 5, maxit = 5)
plot_imputed(imp_interaction_linear, data_interaction)
# imputation method = "norm" with quadratic term
imp_interaction_quadratic <- mice(data_interaction_input, method = "norm", 
                                  formulas = list(data_interaction.X1 = data_interaction.X1 ~ data_interaction.X2 + data_interaction.X3 + I(data_interaction.X2**2) + I(data_interaction.X3**2),
                                                  data_interaction.X2 = data_interaction.X2 ~ data_interaction.X1 + data_interaction.X3 + I(data_interaction.X1**2) + I(data_interaction.X3**2),
                                                  data_interaction.X3 = data_interaction.X3 ~ data_interaction.X1 + data_interaction.X2 + I(data_interaction.X1**2) + I(data_interaction.X2**2)),
                                  m = 5, maxit = 5)
plot_imputed(imp_interaction_quadratic, data_interaction)
# imputation method = "cart"
imp_interaction_cart <- mice(data_interaction_input, method = "cart", m = 5, maxit = 5)
plot_imputed(imp_interaction_cart, data_interaction)


### true DGP: quadratic ###
data_quadratic <- induce_missingness(generate_data(quadratic_coefs = c(0.3,0.7,-0.1)))
summarize(data_quadratic)
data_quadratic_input <- data.frame(data_quadratic$X1,
                                   data_quadratic$X2,
                                   data_quadratic$X3)
# imputation method = "norm"
imp_quadratic_linear <- mice(data_quadratic_input, method = "norm", m = 5, maxit = 5)
plot_imputed(imp_quadratic_linear, data_quadratic)
# imputation method = "norm" with quadratic term
imp_quadratic_quadratic <- mice(data_quadratic_input, method = "norm", 
                                formulas = list(data_quadratic.X1 = data_quadratic.X1 ~ data_quadratic.X2 + data_quadratic.X3 + I(data_quadratic.X2**2) + I(data_quadratic.X3**2),
                                                data_quadratic.X2 = data_quadratic.X2 ~ data_quadratic.X1 + data_quadratic.X3 + I(data_quadratic.X1**2) + I(data_quadratic.X3**2),
                                                data_quadratic.X3 = data_quadratic.X3 ~ data_quadratic.X1 + data_quadratic.X2 + I(data_quadratic.X1**2) + I(data_quadratic.X2**2)),
                                m = 5, maxit = 5)
plot_imputed(imp_quadratic_quadratic, data_quadratic)
# imputation method = "cart"
imp_quadratic_cart <- mice(data_quadratic_input, method = "cart", m = 5, maxit = 5)
plot_imputed(imp_quadratic_cart, data_quadratic) 

### Compute Estimands ###
get_imp_datasets <- function(imp) {
  m <- imp$m 
  imp_datasets <- vector(mode = "list", length = m)
  for (i in 1:m) {
    imp_datasets[[i]] <- complete(imp, i)
    colnames(imp_datasets[[i]]) <- c("X1", "X2", "X3")
  }
  return(imp_datasets)
}

# write function that adds an indicator column to the imputed datasets for whether X2 or X3 > X1
get_imp_prob_greater_than_X1 <- function(imp, col_num) {
  m <- imp$m 
  imp_datasets <- vector(mode = "list", length = m)
  for (i in 1:m) {
    imp_data <- complete(imp, i)
    imp_data$ind <- as.numeric(imp_data[,col_num] > imp_data[,1])
    imp_datasets[[i]] <- imp_data
    colnames(imp_datasets[[i]]) <- c("X1", "X2", "X3", "ind")
  }
  return(imp_datasets)
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


quantile_inputs <- function(imp_datasets, formula, quant) {
  estimates <- suppressWarnings(lapply(imp_datasets, function(x) {coef(rq(formula, data = x, tau = quant))}))
  variances <- suppressWarnings(lapply(imp_datasets, function(x) {summary.rq(rq(formula, data = x, tau = quant), covariance = TRUE)$coefficients[,2]}))
  df <- nrow(imp_datasets[[1]]) - 1
  inputs <- list(estimates = estimates, 
                 variances = variances,
                 df = df)
  return(inputs)
}

pool_results <- function(inputs_df) {
  estimates <- inputs_df$estimates
  variances <- inputs_df$variances
  pooled <- MIcombine(results = estimates, variances = variances, df.complete = inputs_df$df)
  coef_pool <- pooled["coefficients"]$coefficients
  se_pool <- sqrt(diag(pooled[["variance"]]))
  df_imp <- pooled[["df"]]
  fmi <- pooled[["missinfo"]]
  results <- data.frame("est" = coef_pool,
                        "se" = se_pool,
                        "df" = df_imp,
                        "fmi" = fmi)
  print(results)
}

# X2 mean

mean(data_linear$X2_fully_obs) # true value
pool_results(mean_inputs(get_imp_datasets(imp_linear_linear), 2))
pool_results(mean_inputs(get_imp_datasets(imp_linear_quadratic), 2))
pool_results(mean_inputs(get_imp_datasets(imp_linear_cart), 2))

mean(data_quadratic$X2_fully_obs) # true value
pool_results(mean_inputs(get_imp_datasets(imp_quadratic_linear), 2))
pool_results(mean_inputs(get_imp_datasets(imp_quadratic_quadratic), 2))
pool_results(mean_inputs(get_imp_datasets(imp_quadratic_cart), 2))

mean(data_interaction$X2_fully_obs) # true value
pool_results(mean_inputs(get_imp_datasets(imp_interaction_linear), 2))
pool_results(mean_inputs(get_imp_datasets(imp_interaction_quadratic), 2))
pool_results(mean_inputs(get_imp_datasets(imp_interaction_cart), 2))


# X3 mean

mean(data_linear$X3_fully_obs) # true value
pool_results(mean_inputs(get_imp_datasets(imp_linear_linear), 3))
pool_results(mean_inputs(get_imp_datasets(imp_linear_quadratic), 3))
pool_results(mean_inputs(get_imp_datasets(imp_linear_cart), 3))

mean(data_quadratic$X3_fully_obs) # true value
pool_results(mean_inputs(get_imp_datasets(imp_quadratic_linear), 3))
pool_results(mean_inputs(get_imp_datasets(imp_quadratic_quadratic), 3))
pool_results(mean_inputs(get_imp_datasets(imp_quadratic_cart), 3))

mean(data_interaction$X3_fully_obs) # true value
pool_results(mean_inputs(get_imp_datasets(imp_interaction_linear), 3))
pool_results(mean_inputs(get_imp_datasets(imp_interaction_quadratic), 3))
pool_results(mean_inputs(get_imp_datasets(imp_interaction_cart), 3))


# regress X2 on X1

coef(lm(X2_fully_obs ~ X1, data = data_linear))
pool_results(regression_inputs(get_imp_datasets(imp_linear_linear), X2 ~ X1))
pool_results(regression_inputs(get_imp_datasets(imp_linear_quadratic), X2 ~ X1))
pool_results(regression_inputs(get_imp_datasets(imp_linear_cart), X2 ~ X1))

coef(lm(X2_fully_obs ~ X1, data = data_quadratic))
pool_results(regression_inputs(get_imp_datasets(imp_quadratic_linear), X2 ~ X1))
pool_results(regression_inputs(get_imp_datasets(imp_quadratic_quadratic), X2 ~ X1))
pool_results(regression_inputs(get_imp_datasets(imp_quadratic_cart), X2 ~ X1))

coef(lm(X2_fully_obs ~ X1, data = data_interaction))
pool_results(regression_inputs(get_imp_datasets(imp_interaction_linear), X2 ~ X1))
pool_results(regression_inputs(get_imp_datasets(imp_interaction_quadratic), X2 ~ X1))
pool_results(regression_inputs(get_imp_datasets(imp_interaction_cart), X2 ~ X1))


# regress X3 on X1

coef(lm(X3_fully_obs ~ X1, data = data_linear))
pool_results(regression_inputs(get_imp_datasets(imp_linear_linear), X3 ~ X1))
pool_results(regression_inputs(get_imp_datasets(imp_linear_quadratic), X3 ~ X1))
pool_results(regression_inputs(get_imp_datasets(imp_linear_cart), X3 ~ X1))

coef(lm(X3_fully_obs ~ X1, data = data_quadratic))
pool_results(regression_inputs(get_imp_datasets(imp_quadratic_linear), X3 ~ X1))
pool_results(regression_inputs(get_imp_datasets(imp_quadratic_quadratic), X3 ~ X1))
pool_results(regression_inputs(get_imp_datasets(imp_quadratic_cart), X3 ~ X1))

coef(lm(X3_fully_obs ~ X1, data = data_interaction))
pool_results(regression_inputs(get_imp_datasets(imp_interaction_linear), X3 ~ X1))
pool_results(regression_inputs(get_imp_datasets(imp_interaction_quadratic), X3 ~ X1))
pool_results(regression_inputs(get_imp_datasets(imp_interaction_cart), X3 ~ X1))


# X2 quantiles

q <- 0.5

quantile(data_linear$X2_fully_obs, q)
pool_results(quantile_inputs(get_imp_datasets(imp_linear_linear), X2 ~ 1, q))
pool_results(quantile_inputs(get_imp_datasets(imp_linear_quadratic), X2 ~ 1, q))
pool_results(quantile_inputs(get_imp_datasets(imp_linear_cart), X2 ~ 1, q))

quantile(data_quadratic$X2_fully_obs, q)
pool_results(quantile_inputs(get_imp_datasets(imp_quadratic_linear), X2 ~ 1, q))
pool_results(quantile_inputs(get_imp_datasets(imp_quadratic_quadratic), X2 ~ 1, q))
pool_results(quantile_inputs(get_imp_datasets(imp_quadratic_cart), X2 ~ 1, q))

quantile(data_interaction$X2_fully_obs, q)
pool_results(quantile_inputs(get_imp_datasets(imp_interaction_linear), X2 ~ 1, q))
pool_results(quantile_inputs(get_imp_datasets(imp_interaction_quadratic), X2 ~ 1, q))
pool_results(quantile_inputs(get_imp_datasets(imp_interaction_cart), X2 ~ 1, q))


# X3 quantiles

quantile(data_linear$X3_fully_obs, q)
pool_results(quantile_inputs(get_imp_datasets(imp_linear_linear), X3 ~ 1, q))
pool_results(quantile_inputs(get_imp_datasets(imp_linear_quadratic), X3 ~ 1, q))
pool_results(quantile_inputs(get_imp_datasets(imp_linear_cart), X3 ~ 1, q))

quantile(data_quadratic$X3_fully_obs, q)
pool_results(quantile_inputs(get_imp_datasets(imp_quadratic_linear), X3 ~ 1, q))
pool_results(quantile_inputs(get_imp_datasets(imp_quadratic_quadratic), X3 ~ 1, q))
pool_results(quantile_inputs(get_imp_datasets(imp_quadratic_cart), X3 ~ 1, q))

quantile(data_interaction$X3_fully_obs, q)
pool_results(quantile_inputs(get_imp_datasets(imp_interaction_linear), X3 ~ 1, q))
pool_results(quantile_inputs(get_imp_datasets(imp_interaction_quadratic), X3 ~ 1, q))
pool_results(quantile_inputs(get_imp_datasets(imp_interaction_cart), X3 ~ 1, q))


# Pr(X2 > X1)

mean(data_linear$X2_fully_obs > data_linear$X1) # true value
pool_results(mean_inputs(get_imp_prob_greater_than_X1(imp_linear_linear, 2), 4))
pool_results(mean_inputs(get_imp_prob_greater_than_X1(imp_linear_quadratic, 2), 4))
pool_results(mean_inputs(get_imp_prob_greater_than_X1(imp_linear_cart, 2), 4))

mean(data_quadratic$X2_fully_obs > data_quadratic$X1) # true value
pool_results(mean_inputs(get_imp_prob_greater_than_X1(imp_quadratic_linear, 2), 4))
pool_results(mean_inputs(get_imp_prob_greater_than_X1(imp_quadratic_quadratic, 2), 4))
pool_results(mean_inputs(get_imp_prob_greater_than_X1(imp_quadratic_cart, 2), 4))

mean(data_interaction$X2_fully_obs > data_interaction$X1) # true value
pool_results(mean_inputs(get_imp_prob_greater_than_X1(imp_interaction_linear, 2), 4))
pool_results(mean_inputs(get_imp_prob_greater_than_X1(imp_interaction_quadratic, 2), 4))
pool_results(mean_inputs(get_imp_prob_greater_than_X1(imp_interaction_cart, 2), 4))


# Pr(X3 > X1)

mean(data_linear$X3_fully_obs > data_linear$X1) # true value
pool_results(mean_inputs(get_imp_prob_greater_than_X1(imp_linear_linear, 3), 4))
pool_results(mean_inputs(get_imp_prob_greater_than_X1(imp_linear_quadratic, 3), 4))
pool_results(mean_inputs(get_imp_prob_greater_than_X1(imp_linear_cart, 3), 4))

mean(data_quadratic$X3_fully_obs > data_quadratic$X1) # true value
pool_results(mean_inputs(get_imp_prob_greater_than_X1(imp_quadratic_linear, 3), 4))
pool_results(mean_inputs(get_imp_prob_greater_than_X1(imp_quadratic_quadratic, 3), 4))
pool_results(mean_inputs(get_imp_prob_greater_than_X1(imp_quadratic_cart, 3), 4))

mean(data_interaction$X3_fully_obs > data_interaction$X1) # true value
pool_results(mean_inputs(get_imp_prob_greater_than_X1(imp_interaction_linear, 3), 4))
pool_results(mean_inputs(get_imp_prob_greater_than_X1(imp_interaction_quadratic, 3), 4))
pool_results(mean_inputs(get_imp_prob_greater_than_X1(imp_interaction_cart, 3), 4))
