#' @title GLstring_genes
#'
#' @description This function processes a specified column in a data frame
#' that contains GL strings. It separates the GL strings, identifies the HLA
#' loci, and transforms the data into a wider format with loci as column names.
#'
#' @param data A data frame
#' @param gl_string The name of the column in the data frame that contains
#' GL strings
#'
#' @return A data frame with GL strings separated, loci identified, and data
#' transformed to a wider format with loci as columns.
#'
#' @examples
#'
#' file <- HLA_typing_1[, -1]
#' GL_string <- data.frame("GL_string" = HLA_columns_to_GLstring(
#'   file,
#'   HLA_typing_columns = everything()
#' ))
#' GL_string <- GL_string[1, , drop = FALSE] # When considering first patient
#' result <- GLstring_genes(GL_string, "GL_string")
#' print(result)
#'
#' @export
#'
#' @importFrom dplyr select
#' @importFrom dplyr rename
#' @importFrom dplyr mutate
#' @importFrom dplyr %>%
#' @importFrom tidyr separate_longer_delim
#' @importFrom tidyr pivot_wider
#' @importFrom stringr str_extract


GLstring_genes <- function(data, gl_string) {
  # Identify the columns to modify
  col2mod <- names(select(data, {{ gl_string }}))
  data %>%
    # Separate the GL string column by the delimiter "^" into multiple rows
    separate_longer_delim({{ col2mod }}, delim = "^") %>%
    # Rename the separated column to "gl_string"
    rename(gl_string = {{ col2mod }}) %>%
    # Extract the locus information from the GL string
    mutate(locus = str_extract(gl_string, "[[:alnum:]-]+(?=\\*)")) %>%
    # Transform the data from long to wide format, using locus names as new column names
    pivot_wider(names_from = locus, values_from = gl_string) %>%
    # Apply the HLA_column_repair function to the transformed data
    HLA_column_repair(.)
}

globalVariables(c("."))
