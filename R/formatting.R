# Header #############################################################
#
# Author: Lisa Nicvert
# Email:  lisa.nicvert@univ-lyon1.fr
#
# Date: 2024-05-07
#
# Script Description: functions to format data

#' Interaction dataframe to matrix
#'
#' Transforms an interaction dataframe into an incidence matrix.
#'
#' @param df the dataframe to transform. Must have species names
#' (columns which names are given by the parameters `rows`and `columns`),
#' and a third column `frequency` with the interaction frequency/strength index.
#' @param rows Name of the column with the species names that will be the rows of the matrix.
#' @param columns Name of the column with the species names that will be the columns of the matrix.
#'
#' @return Returns an incidence matrix with named rows and columns
#' filled with the values in `frequency`.
#'
#' @export
df_to_matrix <- function(df,
                         rows = "plant_species_code",
                         columns = "animal_species_code"){

  mat_prep <- df |>
    dplyr::group_by(across(all_of(c(rows, columns)))) |>
    dplyr::summarise(frequency = sum(frequency), .groups = "drop")

  # create bird x fruit matrix (l x c)
  mat <- mat_prep |>
    tidyr::pivot_wider(names_from = !!sym(columns),
                       values_from = frequency,
                       values_fill = 0) |>
    tibble::column_to_rownames(rows)

  return(mat)
}


#' Interaction matrix to dataframe
#'
#' Transform interaction matrix (incidence matrix) to a dataframe
#'
#' @param mat Interaction matrix. The function assumes this matrix has row and
#' column names, that are then used to fill the output dataframe columns.
#' @param colnames The name to give to the column
#' containing the matrix column names (defaults to "animals").
#' @param rownames The name to give to the column
#' containing the matrix row names (defaults to "plants").
#' @param tofactor Should the final columns `colnames` and `rownames`
#' be of type factor?
#'
#' @return A dataframe with 3 columns:
#'
#' + `colnames` containing the initial matrix column names
#' + `rownames` containing the initial matrix row names
#' + `value` containing the matrix values
#'
#' @export
matrix_to_df <- function(mat,
                         colnames = "animals",
                         rownames = "plants",
                         tofactor = TRUE) {

  df <- mat |>
    tibble::rownames_to_column(rownames)

  df <- df |>
    tidyr::pivot_longer(cols = 2:ncol(df),
                        names_to = colnames)

  # Factor
  if (tofactor) {
    df[[colnames]] <- factor(df[[colnames]], levels = colnames(mat))
    df[[rownames]] <- factor(df[[rownames]], levels = rownames(mat))
  }

  return(df)
}


#' Filter matrix
#'
#' Filter incidence matrix to keep only species that interact
#' with more partners that a given threshold.
#'
#' @param mat The incidence matrix.
#' @param thr The minimal number of different interacting partners a species
#' must have (a species must have `>= thr` different interacting partners).
#' It is not weighted (for instance, if `thr = 2` then a species interacting 100 times
#' with a unique other partner will be removed.)
#'
#' @return The interaction matrix but where the rows and/or columns with species that do not interact
#' above the threshold removed.
#'
#' @export
filter_matrix <- function(mat, thr){
  # Filter the interaction matrix so as all rows and columns
  # have at least thr interactions with different partners.

  mat_presabs <- ifelse(mat > 0, 1, 0) # if abundance > 0 set to 1

  # Initialize cond
  cond_sumrow <- (apply(mat_presabs, 1, sum) >= thr)
  cond_sumcol <- (apply(mat_presabs, 2, sum) >= thr)

  res <- mat

  while(!(all(cond_sumcol) & all(cond_sumcol))) {
    # Filter out rows, then columns
    if (nrow(res) != 0){
      res <- res[cond_sumrow, ]
    } else {
      message("No data left")
      return(NULL)
    }
    if (ncol(res) != 0){
      res <- res[ , cond_sumcol]
    } else {
      message("No data left")
      return(NULL)
    }

    # Recompute cond
    mat_presabs <- ifelse(res > 0, 1, 0) # if abundance > 0 set to 1

    cond_sumrow <- (apply(mat_presabs, 1, sum) >= thr)
    cond_sumcol <- (apply(mat_presabs, 2, sum) >= thr)
  }

  return(res)
}


#' Create code
#'
#' Create a code for each (unique) value of a vector
#'
#' @param col vector to code
#' @param code_pattern first letter to add before code
#' @param name optional column names for the output dataframe
#'
#' @return dataframe with 2 columns:
#' one with original vector (named `col` by default),
#' the other with the corresponding codes (named `col_code`
#' by default)
#' @export
create_code <- function(col, code_pattern = NA, name = NA){

  ccol <- unique(sort(col))

  dfcode <- data.frame(ccol)

  col_code <- paste(col, "code", sep = "_")

  if(is.na(code_pattern)){
    # Defaults to first letter of wanted name
    if(!is.na(name)){
      code_pattern <- toupper(substr(name, 1, 1))
    }else{
      code_pattern <- ""
    }
  }

  code_number <- 1:nrow(dfcode)
  code_number <- stringr::str_pad(code_number, 2,
                                  pad = '0') # Add leading zero
  dfcode$col_code <- paste0(code_pattern, code_number)

  if(!is.na(name)){
    colnames(dfcode) <- c(name, paste(name, "code", sep = "_"))
  }else{
    colnames(dfcode) <- c("col", "col_code")
  }

  return(dfcode)
}
