# Header #############################################################
#
# Author: Lisa Nicvert
# Email:  lisa.nicvert@univ-lyon1.fr
#
# Date: 2024-05-07
#
# Script Description: function to analyze the outputs of reciprocal scaling (lms)



# Get niche measures ------------------------------------------------------

#' Get mean and standard deviation
#'
#' Return the mean and standard deviation from a reciprocal scaling analysis
#'
#' @param recscal A dataframe (expected to be the output of [`ade4::reciprocal.coa`]).
#' The first columns must contain the reciprocal scaling scores and the last 3 columns
#' are `Row`, `Col` and `Weight`.
#' @param ax the axes for which to compute mean and variance
#' (the function will use the n-th first columns from the `recscal` dataframe).
#'
#' @return A list of 4 matrices, each containing the mean or standard deviation per axes. Each matrix
#' has one column per axis given in `ax` and as many rows as there are grouping levels.
#'
#' + `rowsd` is the standard deviation per row
#' + `rowmean` is the mean per row
#' + `colsd` is the standard deviation per column
#' + `colmean` is the mean per column
#'
#' Matrices rownames are the row or column groups (species names).
#' Column names are the same as in `recscal`.
#'
#' @export
get_mean_sd <- function(recscal, ax = 1:2) {

  if (length(ax) == 1) {
    warning("varfacwt does not work in one dimension with ade4 version 1.7-22")
  }

  # Get mean and variance for columns
  rowvar <- ade4::varfacwt(recscal[, ax], fac = recscal$Row,
                           wt = recscal$Weight)
  rowmean <- ade4::meanfacwt(recscal[, ax], fac = recscal$Row,
                             wt = recscal$Weight)

  # Get mean and variance for rows
  colvar <- ade4::varfacwt(recscal[, ax], fac = recscal$Col,
                           wt = recscal$Weight)
  colmean <- ade4::meanfacwt(recscal[, ax], fac = recscal$Col,
                             wt = recscal$Weight)


  res <- list(rowsd = sqrt(rowvar),
              rowmean = rowmean,
              colsd = sqrt(colvar),
              colmean = colmean)
  return(res)
}


# Linear models ------------------------------------------------------------


#' Get best model
#'
#' Compare 2 linear models and returns the best
#'
#' @param lmsimple First model (if using `method = "LRT"`,
#' must be nested in `lm2`)
#' @param lm2 Second model
#' @param method The method to use to compare models:
#' likelihood ratio test (`method = "LRT"`) for nested linear
#' models comparison or AIC (`method = "AIC"`)
#' @param alpha Significance threshold to consider when using
#' `method = "LRT"`
#'
#' @return The best model
#' @export
get_best_model <- function(lmsimple, lm2, method = c("LRT", "AIC"),
                           alpha = 0.05) {

  if (length(method) > 1) {
    method <- method[1]
  }

  if (method == "LRT") {
    # We assume lmsimple is nested in lm2
    lrt <- lmtest::lrtest(lm2, lmsimple)
    pval <- lrt$`Pr(>Chisq)`[2]

    if(pval <= alpha) {
      # Significant difference: the more complex model is better
      return(lm2)
    } else {
      # Non-significant difference: the more complex model is not better
      return(lmsimple)
    }
  } else if (method == "AIC") {
    # Compute models AIC as 2k - 2 log(likelikood)
    AICs <- AIC(lmsimple, lm2)

    AICsimple <- AICs["lmsimple", "AIC"]
    AIC2 <- AICs["lm2", "AIC"]

    # Then we choose the model with the smaller AIC
    if(AICsimple < AIC2) {
      return(lmsimple)
    } else {
      return(lm2)
    }
  }
}

#' Get the predictions of a linear model
#'
#' Returns the model prediction on given data range.
#'
#' @param dat_predict The data to predict
#' @param lmpred The model. It must have only one explanatory variable named `mean`.
#' @param by The step to use for the range of `dat_predict` values
#' @param level Confidence level used for the prediction. Defaults to 0.95.
#' @param interval The type of interval to use: the argument is used in
#' [`stats::predict.lm`] and has the same interpretation.
#' `confidence` gives the confidence interval around the mean of the observations whereas `prediction` gives
#' the confidence interval of predicted values.
#'
#' @return The model's prediction as a dataframe with columns:
#'
#' + `fit`: the predicted value
#' + `lwr`: lower bound of the confidence interval (see arguments `level` and `interval`)
#' + `upr`: upper bound of the confidence interval (see arguments `level` and `interval`)
#' + `x`: the explanatory variable
#'
#' @export
get_pred <- function(dat_predict, lmpred, by = 0.001, level = 0.95,
                     interval = c("confidence", "prediction")) {

  if (length(interval) > 1) {
    int <- interval[1]
  } else {
    int <- interval
  }

  # Get the range of explanatory variable on which to predict values
  newdat <- seq(min(dat_predict), max(dat_predict), by = by)
  newdat <- data.frame(mean = newdat)

  # Predict values over the range of newdat
  pred <- predict(lmpred, interval = int,
                  level = level,
                  newdata = newdat)

  # Return a dataframe
  pred <- as.data.frame(pred)
  pred$x <- newdat$mean

  return(pred)
}

#' Labels of a linear model
#'
#' Get the labels of a linear model to display on a plot.
#'
#' @param mod The linear model. It is expected to have 2 or 3
#' coefficients of the form \eqn{y = ax + b} or \eqn{y = ax + cx^2 + b}.
#' @param a The axis to consider. Used for the subscript of the variables.
#'
#' @return A dataframe with columns `formula` and `r2` Containing respectively the model equation and coefficient of determination.
#' The formulas are written to be parsed later in the plot.
#' The variables' names are \eqn{s_a} for the predicted value
#' and \eqn{m_a} for the predictor.
#'
#' @details Inspired from <https://r-graphics.org/recipe-annotate-facet>
#'
#'
#' @export
lm_labels <- function(mod, a) {

  coef <- stats::coef(mod)
  if (length(coef) == 3) { # quadratic term, slope and intercept
    formula <- sprintf("italic(s[%.0f]) ==  %.2f * italic(m[%.0f]) %+.2f * italic(m[%.0f]^2) %+.2f",
                       a, coef[2], a, coef[3], a, coef[1])
  } else if (length(coef) == 2) { # slope and intercept
    formula <- sprintf("italic(s[%.0f]) ==  %.2f * italic(m[%.0f]) %+.2f",
                       a, coef[2], a, coef[1])
  }

  r2 <- summary(mod)$adj.r.squared
  r2 <- sprintf("italic(R^2) == %.2f", r2)

  res <- data.frame(formula = formula, r2 = r2)
  return(res)
}
