% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plot.R
\name{plot_reciprocal}
\alias{plot_reciprocal}
\title{Plot reciprocal scaling}
\usage{
plot_reciprocal(
  dudi,
  recscal,
  xax = 1,
  yax = 2,
  labsize = 1,
  psize = 1,
  xlab = NULL,
  ylab = NULL,
  group = c("li", "co"),
  col = "cornflowerblue",
  alpha = 0.2,
  plot_arrows = TRUE,
  plot_points = TRUE,
  plot_labels = TRUE,
  ...
)
}
\arguments{
\item{dudi}{An object of class \code{\link[ade4:dudi]{ade4::dudi}} to pair with the analysis.
Used for the eigenvalues (always) and for the arrows of explanatory variables
if required.
It must be of class \code{coa}, \code{pcaiv} or \code{dpcaiv} (obtained with the new function \code{dpcaiv2}).}

\item{recscal}{A dataframe (expected to be the output of \code{reciprocal.coa} function).
The first columns must contain the reciprocal scaling scores and the last 3 columns
are \code{Row}, \code{Col} and \code{Weight}.}

\item{xax}{The index of the multivariate axis to plot on the x-axis (the column \code{xax} from
the \code{recscal} dataframe).}

\item{yax}{The index of the multivariate axis to plot on the y-axis  (the column \code{yax} from
the \code{recscal} dataframe).}

\item{labsize}{Size of the ellipses labels}

\item{psize}{Size of the points}

\item{xlab}{Custom x-label}

\item{ylab}{Custom y-label}

\item{group}{The group to use for the ellipses (\code{co} or \code{li})}

\item{col}{The ellipses colors}

\item{alpha}{The ellipses transparency}

\item{plot_arrows}{Whether to plot the arrows for the variables}

\item{plot_points}{Plot the data points?}

\item{plot_labels}{Plot ellipses labels?}

\item{...}{Additional parameters passed to \code{s.class}}
}
\value{
A plot with the points grouped by ellipses.
}
\description{
Plot reciprocal scaling
}
