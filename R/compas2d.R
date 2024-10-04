# Header #############################################################
#
# Author: Lisa Nicvert
# Email:  lisa.nicvert@univ-lyon1.fr
#
# Date: 2023-12-15
#
# Script Description: code to generate resource-consumer species interactions.
# (can also be interpreted as species abundances on an environmental gradient).
# adapted from Stéphane Dray (1/12/2004 Montreal from an idea of J Oksanen)


#' Simulate species interactions
#'
#' Simulate species interactions based on the matching between 2 traits.
#' This is based on a code that originally simulated species abundances
#' in different sites based on their environmental niches.
#'
#'
#' @param nconsumer Number of consumer species
#' @param nresource Number of resource species
#' @param le_grad gradient length for traits. The values are comprised between
#' 0 and `le_grad`.
#' @param buffer buffer to allow consumer species optima to fall outside the
#' gradient by 0 - `buffer` ; `le_grad` + `buffer`.
#' @param consumer_ab Optional vector of abundances for consumers.
#' If missing, abundances are taken from a uniform law between `min_ab` and `max_ab`.
#' @param resource_ab Optional vector of abundances for resource species.
#' If missing, abundances are taken from a uniform law between `min_ab` and `max_ab`.
#' @param min_ab minimal abundance value for the uniform law
#' (if `consumer_ab` and/or `resource_ab` are absent)
#' @param max_ab maximal abundance value for the uniform law
#' (if `consumer_ab` and/or `resource_ab` are absent)
#' @param mean_tol Mean value for the distribution of the normal law
#' in which consumers niche standard deviation are drawn (for each trait).
#' @param sd_tol Standard deviation for the distribution of the normal law
#' in which consumers niche standard deviation are drawn (for each trait).
#' @param ratio_grad length of the second gradient as a ratio of the first one.
#' @param col_prefix Prefix for the column names of the matrix (which will be
#' prefix + index).
#' Defaults to "c" for "consumers".
#' @param row_prefix Prefix for the row names of the matrix (which will be
#' prefix + index).
#' Defaults to "r" for "resources".
#' @param rowvar_prefix Prefix for the column names of the row variable
#' (which will be prefix + 1 or 2).
#' Defaults to "tr" for "resource traits".
#' @param remove_zeroes If there are unobserved species,
#' should we keep them in the final matrix?
#' @param ninter Total number of observations in the matrix
#' @param return_intermediate Return intermediate results? If yes, will return the
#' probability of interaction based on matching only, on abundances only and
#' the predicted mix of matching and abundances (before sampling).
#' @param delta Exponent between 0 and 1 to give more (1) or less (0)
#' weight to trait matching relatively to abundance.
#'
#' @details This function models species interactions as arising from the matching between their traits.
#' Traits are drawn from a uniform law with parameters
#' determined by `le_grad` and `ratio_grad` (for the first and second trait).
#' The matching is based on a multivariate normal distribution with means
#' corresponding to the difference between trait values, and variance determined by
#' the variance of consumers' niches.
#' It is then mixed to the abundances (using the parameter `delta`) and sampled to
#' obtain `ninter`interactions.
#'
#' @return An object of class compas, which is a list with 3 components
#' (4 if `return_intermediate`):
#'
#' + `comm`: the interaction matrix of dimension `nresource` x `nconsumer`.
#' Rows and columns are named like species (see `row_prefix` and `col_prefix`).
#' + `rowvar`: resource traits values. It is a `nresource` x 2 matrix.
#' Rows are named like resource species (see `row_prefix`) and
#' column are named like resource traits (see `rowvar_prefix`).
#' + `colvar`: consumer traits values. It is a `nconsumer` x 2 x 2 array.
#' Rows are named like consumer species (see `col_prefix`).
#' Column names in each array dimension are `mean` and `sd`
#' (for the mean and standard deviation of the species niche, respectively).
#' + if `return_intermediate`, a fourth component named intermediate is returned.
#' The first element is named `p_matching` and contains the matrix of
#' probabilities of interactions based only on matching.
#' The second `ab_neutral` contains the matrix of count of interactions as
#' predicted by species abundances.
#' The third `p_mix` contains the matrix of probabilities of
#' interactions taking into account the mix of matching and abundances.
#'
#' @export
"compas2d" <-
function (nconsumer = 40,  nresource = 100,
          le_grad = 100, ratio_grad = 0.8,
          consumer_ab = NULL,
          resource_ab = NULL,
          min_ab = 1, max_ab = 100,
          ninter = 100,
          delta = 1,
          mean_tol = 2, sd_tol = 10,
          buffer = 1,
          col_prefix = "c", row_prefix = "r",
          rowvar_prefix = "tr",
          remove_zeroes = FALSE,
          return_intermediate = FALSE)
{

  # Gradient for axis 2 ---
  # The length of the second gradient is a fraction of the length of the
  # first gradient
  gradmin2 <- (1-ratio_grad) / 2 * le_grad
  gradmax2 <- le_grad - (1-ratio_grad) / 2 * le_grad
  # For the columns, we include a buffer that is also on axis 2
  # (if buffer = 0 then gradmin/max2_buffer and gradmin/max2 are equal)
  gradmin2_buffer <- -buffer + gradmin2
  gradmax2_buffer <- buffer + gradmax2

  # Environment/resource species trait gradient ---
  # Initialize empty trait/environment values matrix
  x <- matrix(0, nrow = nresource, ncol = 2)
  # For each site, generate an environmental gradient value at random (uniform) (and sort them)
  x[, 1] <- sort(runif(nresource, min = 0, max = le_grad))
  x[, 2] <- runif(nresource, min = gradmin2, max = gradmax2)

  # Abundances ---
  # Species/consumer species
  if (is.null(consumer_ab)) {
    # Generate random abundances for each species (uniform law)
    consumer_ab <- runif(nconsumer, min = min_ab, max = max_ab)
  }
  # Sites/resource species
  if (is.null(resource_ab)) {
    # Generate random abundances for each site (uniform law)
    resource_ab <- runif(nresource, min = min_ab, max = max_ab)
  }

  # Species niche/consumer niche ---
  # Initialize an array for species optima and tolerances
  # The array has last dimension 2 (one for each trait)
  p <- array(0, dim = c(nconsumer, 2, 2))

  # Fill array p

  # -> first dimension
  # Generate environmental optima for each species
  p[, 1, 1] <- runif(nconsumer, min = 0-buffer, max = le_grad+buffer)
  # Generate random niche widths (sd) for each species
  p[, 2, 1] <- abs(rnorm(nconsumer, mean = mean_tol, sd = sd_tol))

  # -> second dimension
  # Optima
  p[, 1, 2] <- runif(nconsumer,
                     min = gradmin2_buffer,
                     max = gradmax2_buffer)
  # Widths (sd)
  p[, 2, 2] <- abs(rnorm(nconsumer, mean = mean_tol, sd = sd_tol))

  # Probability matrix (only matching) ---
  # Initialize empty community matrix
  p_matching <- matrix(0, nrow = nrow(x), ncol = nconsumer)

  for (i in 1:nrow(x)) {
    for (j in 1:nconsumer) {
      # Fill each cell with a "presence probability" or an "interaction probability"
      # from a multivariate normal law depending on:
      # - the matching between resource species trait and consumer species trait
      # - the matching between the site environmental value and the species niche optimum on this gradient
      p_matching[i,j] <- mvtnorm::dmvnorm(x[i,] - p[j, 1, ], sigma = diag(p[j, 2, ]^2))
    }
  }
  # Transform negative values to zero probability
  p_matching <- ifelse(p_matching > 0, p_matching, 0)

  # Make columns a proba distribution
  # - each species (column) is distributed in sites following a probability of choosing this site
  # - each consumer (column) chooses the resource according to a given proba of presence
  p_matching <- sweep(p_matching, 2, STATS = apply(p_matching, 2, sum), FUN = "/")

  # Quick patch (in case there are species with zero obs that became NA at the division step)
  p_matching[is.na(p_matching)] <- 0

  # Theoretical interaction count (based on abundances) ---
  # This code makes sense only for interaction matrices because the abundance of
  # resource species is a limiting factor here.
  # It makes less sense in the context of species x sites association, unless we consider
  # resources in the different sites are limited and a site has a limited number of species.
  prop_row <- resource_ab/sum(resource_ab) # Get the proportion of each plant (its "availability")
  ab_neutral <- prop_row %*% t(consumer_ab) # Get the predicted abundance -> the birds
  # choose a plant only based on its availability
  # Each bird abundance is exactly consumer_ab and the plants abundances are proportional to
  # resource_ab

  # Probability matrix (matching and abundance) ---
  # these theoretical abundances of interactions must then be confronted to the probability of interactions
  ab_mix <- ab_neutral*p_matching^delta # We multiply expected interaction number by trait matching
  # Create a vector of probabilities that takes into account the matching
  p_mix <- ab_mix/sum(ab_mix)
  p_mix[is.na(p_mix)] <- 0 # Handle divisions per zero

  # Final community matrix (with counts) ---
  # Then we need to sample observations in this web
  # The method is inconsumerired by Fründ et al 2016
  # We sample ninter interactions from a multinomial distribution, where the probability to
  # draw each interaction depends on the probability taking into account abundance and matching
  ab_obs_vec <- stats::rmultinom(n = 1,
                                 size = ninter,
                                 prob = as.numeric(p_mix))
  ab_obs <- matrix(ab_obs_vec, nrow = nresource, byrow = FALSE) # Reformat to a matrix

  if (remove_zeroes) {
    # Remove species that haven't been observed
    colkeep <- which(colSums(ab_obs) != 0)
    rowkeep <- which(rowSums(ab_obs) != 0)

    ab_obs <- ab_obs[, colkeep]
    ab_obs <- ab_obs[rowkeep, ]

    # Select traits only species that have been observed
    p <- p[colkeep, , ]
    x <- x[rowkeep, ]
  }

  # Format result
  # Rows and columns names
  rnames <- paste0(row_prefix, 1:nrow(ab_obs))
  cnames <- paste0(col_prefix, 1:ncol(ab_obs))

  rownames(ab_obs) <- rnames
  colnames(ab_obs) <- cnames

  rownames(x) <- rnames
  colnames(x) <- paste0(rowvar_prefix, 1:2)

  rownames(p) <- cnames
  colnames(p) <- c("mean", "sd")

  if (return_intermediate) {
    rownames(p_matching) <- rnames
    colnames(p_matching) <- cnames

    rownames(ab_neutral) <- rnames
    colnames(ab_neutral) <- cnames

    rownames(p_mix) <- rnames
    colnames(p_mix) <- cnames

    intermediate <- list(p_matching = p_matching,
                         ab_neutral = ab_neutral,
                         p_mix = p_mix)

    sol <- list(comm = ab_obs,
                rowvar = x,
                colvar = p,
                intermediate = intermediate)
  } else {
    sol <- list(comm = ab_obs,
                rowvar = x,
                colvar = p)
  }


  class(sol) <- "compas"
  sol$call <- match.call()

  return(sol)
}
