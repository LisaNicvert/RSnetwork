% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/analysis.R
\name{get_pred}
\alias{get_pred}
\title{Get the predictions of a linear model}
\usage{
get_pred(
  dat_predict,
  lmpred,
  by = 0.001,
  level = 0.95,
  interval = c("confidence", "prediction")
)
}
\arguments{
\item{dat_predict}{The data to predict}

\item{lmpred}{The model. It must have only one explanatory variable named \code{mean}.}

\item{by}{The step to use for the range of \code{dat_predict} values}

\item{level}{Confidence level used for the prediction. Defaults to 0.95.}

\item{interval}{The type of interval to use: the argument is used in
\code{\link[stats:predict.lm]{stats::predict.lm}} and has the same interpretation.
\code{confidence} gives the confidence interval around the mean of the observations whereas \code{prediction} gives
the confidence interval of predicted values.}
}
\value{
The model's prediction as a dataframe with columns:
\itemize{
\item \code{fit}: the predicted value
\item \code{lwr}: lower bound of the confidence interval (see arguments \code{level} and \code{interval})
\item \code{upr}: upper bound of the confidence interval (see arguments \code{level} and \code{interval})
\item \code{x}: the explanatory variable
}
}
\description{
Returns the model prediction on given data range.
}
