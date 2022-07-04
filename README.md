# multiple-imputation
### Fall 2021 &amp; Spring 2022 Research Elective with Jared Murray

The goal of this project was to learn about multiple imputation for missing data and explore potential diagnostic tools for imputation models.

## Data Simulation

This script contains code to generate n observations for X1, X2, and X3. X1 is fully observered and generated from a normal distribution. X2 is regressed on X1 with an optional quadratic term, and X3 is regressed on X1 and X2 with optional quadratic terms and an interaction. A function then induces missingness on X2 and X3. 

The code includes examples where the true data generating process is linear, has an interaction term, and has a quadratic term. For each of these, we generate multiply imputed datasets using multiple imputation by chained equations (MICE). We experiment with imputation models based on linear regression, including an interaction term, and using CART, to see the results when imputation models are correctly or incorrectly specified for the three different data generating processes.

Lastly, we computed pooled estimands for each of the different combinations and compare them against the true values. Estimands include univariate means, regression coefficients, quantiles, and the probability that X2 and X3 are bigger than X1.

## Posterior Predictive Checks

This script implements the posterior predictive checks for an imputation model described in He and Zaslavsky (2012). After generating some data with quadratic and interaction terms and then inducing missingness on X2 and X3, the data is duplicated and concatenated but with X2 and X3 completely missing. 

100 multiply imputed datasets are generated for the stacked dataset using MICE. Then, posterior predictie p-values of an estimand are calculated by finding the proportion of estimands that are larger in the original data with imputations than the fully imputed duplicated data. We calculate p-values for univariate means and regression coefficients.

Since posterior predictive checks rely on the particular estimand one chosen, they may fail to reveal the shortcomings of an imputation model. For example, quadratic data with a linear imputation model may produce reasonable p-values for univariate means and basic linear regression, but the same data results in extreme p-values when regressing on quadratic terms due to the misspecification of the imputation model. Posterior predictive checks have limitations as a diagnostic tool since the imputer doesn't know what estimands the analyst will choose to compute.

For a more flexible diagnistic tool, we next wanted to train a classifier on the duplicated and concatenated data used in the posterior predictive checks. For the 100 multiply imputed datasets, we added an indicator column where 1 = original data and 0 = duplicated data, then we combined all of the MI datasets into a single dataframe and trained a classifier to predict whether a data point is original or from the duplicated data. If a classifier can accurately predict whether the data was original versus duplicated, then this indicates where the imputation model may be misspecified.

Classifiers used include CART, logistic regression, and random forests. Various plots illustrate how the classifiers can identify original data compared to data imputed by a misspecified model. In addition to classifiers trained on original versus duplicated, we also include classifiers to predict missingness in the original dataset for X2 and X3. In general, the flexibility of the random forest model can make quite accurate predictions and visually reveal the failings of the misspecified model, compared to the posterior predictive checks that depended on chosen estimands.

## Counterfactuals

As I was fitting classifiers to imputed data, I had the idea that counterfactuals could be a useful tool for diagnosing an imputation model. Counterfactual explanations are a method in interpretable ML that ask how the inputs would have had to vary in order to receive a different prediction from the model. Applying this idea to an imputation model, if we can generate a counterfactual explanation that says which variables should change in order to predict that the data point was imputed, that would reveal where the imputation model is possibly misspecified.

The attached slides describe in more detail the method for generating counterfactuals from "Counterfactual Explanations Without Opening the Black Box: Automated Decisions and the GDPR" (Wachter, Mittelstadt & Russell, 2018). The R script implements this method on the Pima dataset. While I didn't have time during the spring 2022 semester to extend this method to the classifiers fit to imputed data, that would be a possible area for future research.

