# Header #############################################################
#
# Author: Lisa Nicvert
# Email:  lisa.nicvert@univ-lyon1.fr
#
# Date: 2024-05-07
#
# Script Description: functions to analyze simulations results



# Niche measures ----------------------------------------------------------

#' Get niches
#'
#' Returns the fundamental and realized niches and the niches returned
#' by reciprocal scaling.
#'
#' @param rec A dataframe (expected to be the output of [`reciprocal.coa()`]).
#' The first 2 columns must contain the reciprocal scaling scores (there can be more axes, but only
#' the first 2 will be used.) The last 3 columns are `Row`, `Col` and `Weight`.
#' @param comm The observed community matrix. It is used to get the
#' resource species niche width.
#' @param consumer_niche A nx2x2 array with n rows (n species),
#' 2 columns corresponding to species niche optima and standard deviation
#' (corresponding to their multivariate normal niche) and 2 traits stored
#' in the dimension.
#' There is 1 row per species.
#' @param resource_traits A mx2 matrix. Each row corresponds to a
#' species and each column to a trait value.
#' @param rowname How to name the first element of the list
#' (resources/rows) (e.g. could be "plants")
#' @param colname How to name the second element of the list
#' (resources/rows) (e.g. could be "birds")
#'
#' @return A list of lists of lists. The first element (named `rowname`) contains the niches for the
#' resources and the second element (named `colname`) contains the niches for the consumers.
#' Each element contains 2 lists:
#'
#' + `mean`: niche optima
#' + `sd`: niche breadths
#'
#' For each of these sublists, there are 3 matrices elements named:
#'
#' + `fundamental`: the fundamental niche based on species traits
#' + `realized`: the realized niche based on interacting species traits
#' + `mvar`: niche parameters corresponding to the mean and standard deviations
#' given by reciprocal scaling.
#'
#' These matrices describe the value for each species (in rows)
#' along each trait/axis (in column).
#'
#' @details
#' The fundamental niche optimum of a consumer is the mean of its niche
#' (from the `consumer_niche` table).
#' For a resource, it corresponds to the traits from the `resource_traits` table.
#' The fundamental niche breadths are computed for consumers only,
#' and corresponds to the breadths of their niche (from the `consumer_niche` table).
#'
#' The realized niche optima are computed as the weighted means of the traits
#' of species' interacting partners: therefore, the realized niche
#' of a consumer is the weighted mean of the traits of the resource species
#' it interacts with, and reciprocally.
#' The realized niche breadths are computed as the weighted standard deviations
#' of the traits of a species' interacting partners. Therefore, the realized
#' niche breadth of a consumer is the weighted standard deviation of the
#' traits of the resource species it interacts with, and reciprocally.
#'
#' Niche optima computed with multivariate methods correspond to the reciprocal
#' scaling mean and niche breadths correspond to the reciprocal scaling
#' standard deviation.
#'
#' @export
get_niches <- function(rec, comm,
                       consumer_niche, resource_traits,
                       rowname = "row", colname = "col") {

  # ---
  # Resources
  # ---

  # Mean ---
  # Fundamental niche
  mean_fu_row <- resource_traits
  mean_fu_row <- as.matrix(mean_fu_row)
  colnames(mean_fu_row) <- paste0("t", 1:ncol(mean_fu_row))

  # Realized niche
  mean_re_row <- apply(comm, 1,
                       function(x) meanfacwt(consumer_niche[, "mean", ], wt = x))
  mean_re_row <- as.matrix(t(mean_re_row))
  colnames(mean_re_row) <- paste0("t", 1:ncol(mean_re_row))

  # Reciprocal scaling
  mean_mvar_row <- meanfacwt(rec[, 1:2],
                             fac = rec$Row,
                             wt = rec$Weight)
  mean_mvar_row <- as.matrix(mean_mvar_row)
  colnames(mean_mvar_row) <- paste0("ax", 1:ncol(mean_mvar_row))

  # Standard deviation ---
  # No fundamental variance for resources

  # Realized niche
  var_re_row <- apply(comm, 1,
                      function(x) varfacwt(consumer_niche[, "mean", ], wt = x))
  sd_re_row <- sqrt(var_re_row)
  # It is the same as computing the divc below
  # var_true_row <- sapply(1:2,
  #                        function(i) {
  #                          divc(data.frame(t(comm)),
  #                               dis = dist(as.matrix(consumer_niche[, "mean", i])))[[1]]})
  sd_re_row <- as.matrix(t(sd_re_row))
  colnames(sd_re_row) <- paste0("t", 1:ncol(sd_re_row))

  # Reciprocal scaling
  var_mvar_row <- varfacwt(rec[, 1:2],
                           fac = rec$Row,
                           wt = rec$Weight)
  sd_mvar_row <- sqrt(var_mvar_row)
  sd_mvar_row <- as.matrix(sd_mvar_row)
  colnames(sd_mvar_row) <- paste0("ax", 1:ncol(sd_mvar_row))

  # ---
  # Consumers
  # ---

  # Mean ---
  # Fundamental niche
  mean_fu_col <- consumer_niche[, "mean", ]
  mean_fu_col <- as.matrix(mean_fu_col)
  colnames(mean_fu_col) <- paste0("t", 1:ncol(mean_fu_col))

  # Realized niche
  mean_re_col <- apply(comm, 2,
                       function(x) meanfacwt(resource_traits, wt = x))
  mean_re_col <- as.matrix(t(mean_re_col))
  colnames(mean_re_col) <- paste0("t", 1:ncol(mean_re_col))

  # Reciprocal scaling
  mean_mvar_col <- meanfacwt(rec[, 1:2],
                             fac = rec$Col,
                             wt = rec$Weight)
  mean_mvar_col <- as.matrix(mean_mvar_col)
  colnames(mean_mvar_col) <- paste0("ax", 1:ncol(mean_mvar_col))

  # Standard deviation ---
  # Fundamental niche
  sd_fu_col <- consumer_niche[, "sd", ]
  sd_fu_col <- as.matrix(sd_fu_col)
  colnames(sd_fu_col) <- paste0("t", 1:ncol(sd_fu_col))

  # Realized niche
  var_re_col <- apply(comm, 2,
                       function(x) varfacwt(resource_traits, wt = x))
  sd_re_col <- sqrt(var_re_col)
  sd_re_col <- as.matrix(t(sd_re_col))
  colnames(sd_re_col) <- paste0("t", 1:ncol(sd_re_col))

  # Reciprocal scaling
  var_mvar_col <- varfacwt(rec[, 1:2],
                           fac = rec$Col,
                           wt = rec$Weight)
  sd_mvar_col <- sqrt(var_mvar_col)
  sd_mvar_col <- as.matrix(sd_mvar_col)
  colnames(sd_mvar_col) <- paste0("ax", 1:ncol(sd_mvar_col))

  # ---
  # Result list
  # ---
  res <- list(list("mean" = list("fundamental" = mean_fu_row,
                                 "realized" = mean_re_row,
                                 "mvar" = mean_mvar_row),
                   "sd" = list("realized" = sd_re_row,
                                "mvar" = sd_mvar_row)),
              list("mean" = list("fundamental" = mean_fu_col,
                                 "realized" = mean_re_col,
                                 "mvar" = mean_mvar_col),
                   "sd" = list("fundamental" = sd_fu_col,
                                "realized" = sd_re_col,
                                "mvar" = sd_mvar_col)))
  names(res) <- c(rowname, colname)

  return(res)
}

#' Get correlations
#'
#' Measure correlations between niche obtained from CA or reciprocal scaling
#' and the fundamental and realized niches.
#'
#' @param niches An object containing the fundamental, realized and multivariate niches
#' (output of [get_niches()])
#'
#' @return A list of lists. The first element (named `fundamental`) contains the correlations with the
#' fundamental niche values and the second element (named `realized`)
#' contains the correlations with the realized niche values.
#' Each sublist is a list of 2 dataframes:
#'
#' + `mean_cor`: correlations of niches optima with the CA coordinates.
#' + `sd_cor`: correlations of niches breadths with reciprocal scaling standard deviations.
#'
#' These dataframes have the following columns:
#'
#' + `axis`: correlations measured on which multivariate axis?
#' + `type`: the type of individuals for which the measure is (rows or columns?).
#' Column or row codes are chosen with `rowname` and `colname`
#' from the parameter `niches`.
#' + `cor`: the correlation value
#'
#' @export
get_cor <- function(niches) {

  # Get names from the niches object
  rowname <- names(niches)[1]
  colname <- names(niches)[2]

  # ---
  # Means
  # ---
  mean_col <- niches[[colname]]$mean
  mean_row <- niches[[rowname]]$mean

  # Correlation matrix (fundamental) ---
  cormat_col_fu <- cor(mean_col$fundamental,
                       mean_col$mvar)
  cormat_row_fu <- cor(mean_row$fundamental,
                       mean_row$mvar)

  # Correlation matrix (realized) ---
  cormat_col_re <- cor(mean_col$realized,
                       mean_col$mvar)

  cormat_row_re <- cor(mean_row$realized,
                       mean_row$mvar)

  # Choose trait/axis pair ---
  # Since trait 1 explains more variance, it should be correlated more with axis 1. In that case,
  # we will measure the correlation between trait 1 and axis 1. But if it correlates with axis 2,
  # then we will measure the correlation between trait 1 and axis 2.

  # The test is based on the realized niche and the first trait of the column/consumer.
  # It will determine how all other correlations are measured.
  bigger_t1_ax1 <- abs(cormat_col_re[1, 1]) > abs(cormat_col_re[2, 1])

  # Get correlations ---
  if (bigger_t1_ax1) { # Correlation for (column) trait 1 larger on axis 1: don't invert axes
    # Fundamental niche
    cor_mean_col_fu <- abs(c(cormat_col_fu[1, 1], cormat_col_fu[2, 2]))
    cor_mean_row_fu <- abs(c(cormat_row_fu[1, 1], cormat_row_fu[2, 2])) # same axis for the rows

    # Realized niche
    cor_mean_col_re <- abs(c(cormat_col_re[1, 1], cormat_col_re[2, 2]))
    cor_mean_row_re <- abs(c(cormat_row_re[1, 1], cormat_row_re[2, 2]))

  } else { # Correlation for (column) trait 1 larger on axis 2
    # Fundamental niche
    cor_mean_col_fu <- abs(c(cormat_col_fu[2, 1], cormat_col_fu[1, 2])) # order: axis 1 then axis 2
    cor_mean_row_fu <- abs(c(cormat_row_fu[2, 1], cormat_row_fu[1, 2])) # same axis for the rows

    # Realized niche
    cor_mean_col_re <- abs(c(cormat_col_re[2, 1], cormat_col_re[1, 2]))
    cor_mean_row_re <- abs(c(cormat_row_re[2, 1], cormat_row_re[1, 2]))
  }

  # ---
  # Columns standard deviation
  # ---
  sd_col <- niches[[colname]]$sd

  sd_fu_col <- sd_col$fundamental
  sd_re_col <- sd_col$realized
  sd_col <- sd_col$mvar

  # Get correlations ---
  # Here again the test is based on the first trait of the column/consumer and
  # the CA mean position.
  if (bigger_t1_ax1) {
    cor_sd_col_fu <- abs(c(cor(sd_fu_col[, 1], sd_col[, 1]),
                           cor(sd_fu_col[, 2], sd_col[, 2])))
    cor_sd_col_re <- abs(c(cor(sd_re_col[, 1], sd_col[, 1]),
                           cor(sd_re_col[, 2], sd_col[, 2])))
  } else {
    cor_sd_col_fu <- abs(c(cor(sd_fu_col[, 2], sd_col[, 1]), # axis 1
                           cor(sd_fu_col[, 1], sd_col[, 2]))) # axis 2
    cor_sd_col_re <- abs(c(cor(sd_re_col[, 2], sd_col[, 1]), # axis 1
                           cor(sd_re_col[, 1], sd_col[, 2]))) # axis 2
  }

  # ---
  # Rows standard deviation
  # ---
  sd_row <- niches[[rowname]]$sd

  sd_re_row <- sd_row$realized
  sd_row <- sd_row$mvar

  # Get correlations ---
  if (bigger_t1_ax1) {
    cor_sd_row_re <- abs(c(cor(sd_re_row[, 1], sd_row[, 1]),
                           cor(sd_re_row[, 2], sd_row[, 2])))
  } else {
    cor_sd_row_re <- abs(c(cor(sd_re_row[, 2], sd_row[, 1]), # axis 1
                           cor(sd_re_row[, 1], sd_row[, 2]))) # axis 2
  }

  # ---
  # Format data
  # ---
  dfimean_fu <- data.frame(axis = rep(c(1, 2), 2),
                           type = rep(c(rowname, colname), each = 2), # Resource and consumer
                           cor = c(cor_mean_row_fu, cor_mean_col_fu))
  dfisd_fu <- data.frame(axis = c(1, 2),
                          type = rep(colname, 2), # Only consumer
                          cor = cor_sd_col_fu)

  dfimean_re <- data.frame(axis = rep(c(1, 2), 2),
                           type = rep(c(rowname, colname), each = 2), # Resource and consumer
                        cor = c(cor_mean_row_re, cor_mean_col_re))
  dfisd_re <- data.frame(axis = rep(c(1, 2), 2),
                          type = rep(c(rowname, colname), each = 2), # Resource and consumer
                          cor = c(cor_sd_row_re, cor_sd_col_re))

  res <- list("fundamental" = list("mean_cor" = dfimean_fu,
                                   "sd_cor" = dfisd_fu),
              "realized" = list("mean_cor" = dfimean_re,
                                "sd_cor" = dfisd_re))
  return(res)
}


# Helpers -----------------------------------------------------------------

#' Add parameters
#'
#' Adds parameters value to simulation results dataframe.
#'
#' @param corsim a list with the correlation values (output of [get_cor()]).
#' @param paramvec A vector of parameters (must be a 1-row dataframe)
#'
#' @return The list `corsim` with added columns in each dataframe of the list.
#' The added columns are named like the columns of `paramvec` and contain repeated values
#' of `paramvec`
#'
#' @export
add_params <- function(corsim, paramvec) {

  res <- corsim

  pcol <- paramvec[rep(1, 4), ]

  pcol_len2 <- paramvec[rep(1, 2), ]

  res$fundamental$mean_cor[, names(paramvec)] <- pcol
  res$fundamental$sd_cor[, names(paramvec)] <- pcol_len2

  res$realized$mean_cor[, names(paramvec)] <- pcol
  res$realized$sd_cor[, names(paramvec)] <- pcol

  return(res)
}

#' Combine correlations
#'
#' Combine lists formatted by [get_cor()].
#'
#' @param a First list formatted by [get_cor()]
#' @param b Second list formatted by [get_cor()]
#'
#' @return A list as formatted by [get_cor()] (see this function's
#' "Return" section).
#' @export
combine_cor <- function(a, b) {

  # Function to
  res <- list("fundamental" = list("mean_cor" = data.frame(),
                                   "sd_cor" = data.frame()),
              "realized" = list("mean_cor" = data.frame(),
                                "sd_cor" = data.frame()))

  res$fundamental$mean_cor <- rbind(a$fundamental$mean_cor,
                                    b$fundamental$mean_cor)
  res$fundamental$sd_cor <- rbind(a$fundamental$sd_cor,
                                  b$fundamental$sd_cor)
  res$realized$mean_cor <- rbind(a$realized$mean_cor,
                                 b$realized$mean_cor)
  res$realized$sd_cor <- rbind(a$realized$sd_cor,
                               b$realized$sd_cor)

  return(res)
}

