# Header #############################################################
#
# Author: Lisa Nicvert
# Email:  lisa.nicvert@univ-lyon1.fr
#
# Date: 2024-05-07
#
# Script Description: reciprocal scaling function for CA. Function coded by Stéphane Dray in 2023.

#' Reciprocal scaling for CA
#'
#' Performs reciprocal scaling from corrrespondence analysis
#'
#' @param x output of `dudi.coa`
#'
#' @return Results of the reciprocal scaling.
#' Each row corresponds to a correspondence in the original table (a nonzero occurrence).
#'
#' The first columns (`Scorexx`) give coordinates in the different dimensions.
#' `Scorexx` gives the coordinates of the correspondence in the multivariate
#' space (as given by canonical correlation analysis).
#'
#' `Row` and `Col` give the row and column this correspondence belongs to.
#'
#' `Weight` gives the correspondence weight (the count of the original cell divided
#' by the sum of the table).
#'
#' @references See the original publication: Thioulouse, J., & Chessel, D. (1992).
#' A Method for Reciprocal Scaling of Species Tolerance and Sample Diversity.
#' Ecology, 73(2), 670–680. \doi{doi:10.2307/1940773}

#'
#' @export
"reciprocal.coa" <- function (x) {

  # Check type
    if (!inherits(x, "coa"))
        stop("Object of class 'coa' expected")

  # Retrieve original data table oritab
    if (inherits(x, "witwit")) {
        y <- eval.parent(as.list(x$call)[[2]])
        oritab <- eval.parent(as.list(y$call)[[2]])
    }
    else oritab <- eval.parent(as.list(x$call)[[2]])
    l.names <- row.names(oritab)
    c.names <- names(oritab)
    oritab <- as.matrix(oritab)

    f1 <- function(x,oritab,xax){
      # Function to compute double coa coordinates starting from
      # coa coordinates (eq 11 in Thioulouse & Chessel)
      # ### Inputs
      # x: dudi.coa object
      # oritab: original data table
      # xax: axis to compute coordinates from
      a <- x$co[col(oritab), xax]
      a <- a + x$li[row(oritab), xax]
      a <- a/sqrt(2 * x$eig[xax] * (1 + sqrt(x$eig[xax])))
      a <- a[oritab > 0]
    }

    # Compute new coordinates from coa table
    res <- sapply(1:x$nf,f1,x=x,oritab=oritab)
    # Get names of columns/rows where there is at least one individual
    aco <- c.names[col(oritab)[oritab > 0]]
    aco <- factor(aco, levels = c.names)
    ali <- l.names[row(oritab)[oritab > 0]]
    ali <- factor(ali, levels = l.names)

    # Get weights for each observation -> relative frequency
    # of each combination of site/species
    aw <- oritab[oritab > 0]/sum(oritab)
    res <- cbind.data.frame(res,Row=ali,Col=aco,Weight=aw)
    names(res)[1:x$nf] <- paste("Scor",1:x$nf,sep="")
    rownames(res) <- paste(ali,aco,sep="")
    return(res)
}
