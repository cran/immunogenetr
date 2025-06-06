#' @title HLA_mismatch_number
#'
#' @description Calculates the number of mismatched HLA alleles between a
#' recipient and a donor across specified loci. Supports mismatch calculations
#' for host-vs-graft (HvG), graft-vs-host (GvH), or bidirectional.
#'
#' @param GL_string_recip A GL string representing the recipient's HLA genotype.
#' @param GL_string_donor A GL string representing the donor's HLA genotype.
#' @param loci A character vector specifying the loci to be considered for
#' mismatch calculation. HLA-DRB3/4/5 (and their serologic equivalents DR51/52/53)
#' are considered once locus for this function, and should be called in this argument
#' as "HLA-DRB3/4/5" or "HLA-DR51/52/53", respectively.
#' @param direction A character string indicating the direction of mismatch.
#' Options are "HvG" (host vs. graft), "GvH" (graft vs. host), "bidirectional"
#' (the max value of "HvG" and "GvH"), or "SOT" (host vs. graft, as is used for
#' mismatching in solid organ transplantation).
#' @param homozygous_count An integer specifying how to count homozygous mismatches.
#' Defaults to 2, where homozygous mismatches are treated as two mismatches,
#' regardless if one or two alleles are supplied in the GL string (in cases
#' where one allele is supplied, it is duplicated by the function). If
#' specified as 1, homozygous mismatches are only counted once, regardless of
#' whether one or two alleles are supplied in the GL string (in cases where
#' two alleles are supplied, the second identical allele is deleted).
#'
#' @return An integer value or a character string:
#' - If `loci` includes only one locus, the function returns an integer
#' mismatch count for that locus.
#' - If `loci` includes multiple loci, the function returns a character
#' string in the format "Locus1=Count1, Locus2=Count2, ...".
#'
#' @examples
#'
#' file <- HLA_typing_1[, -1]
#' GL_string <- HLA_columns_to_GLstring(file, HLA_typing_columns = everything())
#'
#' GL_string_recip <- GL_string[1]
#' GL_string_donor <- GL_string[2]
#'
#' loci <- c("HLA-A", "HLA-DRB3/4/5", "HLA-DPB1")
#'
#' # Calculate mismatch numbers (Host vs. Graft)
#' HLA_mismatch_number(GL_string_recip, GL_string_donor, loci, direction = "HvG")
#'
#' # Calculate mismatch numbers (Graft vs. Host)
#' HLA_mismatch_number(GL_string_recip, GL_string_donor, loci, direction = "GvH")
#'
#' # Calculate mismatch numbers (Bidirectional)
#' HLA_mismatch_number(GL_string_recip, GL_string_donor,
#'   loci,
#'   direction = "bidirectional"
#' )
#'
#' @export
#'
#' @importFrom dplyr join_by
#' @importFrom dplyr mutate
#' @importFrom dplyr summarise
#' @importFrom dplyr left_join
#' @importFrom dplyr na_if
#' @importFrom stringr str_count
#' @importFrom stringr str_flatten
#' @importFrom tidyr separate_longer_delim
#' @importFrom tidyr separate_wider_delim
#' @importFrom tidyr replace_na
#' @importFrom tibble tibble
#' @importFrom tidyr unite


HLA_mismatch_number <- function(GL_string_recip, GL_string_donor, loci, direction, homozygous_count = 2) {
  direction <- match.arg(direction, c("HvG", "GvH", "bidirectional", "SOT"))
  # Code to determine mismatch numbers if a single locus was supplied.
  if (length(loci) == 1) {
    # Determine mismatches for both directions.
    HvG <- replace_na(str_count(HLA_mismatch_base(GL_string_recip, GL_string_donor, loci, "HvG", homozygous_count), "(\\+|$)"), 0) # The regex matches the end of the string or a "+".
    GvH <- replace_na(str_count(HLA_mismatch_base(GL_string_recip, GL_string_donor, loci, "GvH", homozygous_count), "(\\+|$)"), 0)
    # Make a tibble with the results and determine bidirectional mismatch.
    MM_table <- tibble(HvG, GvH) %>%
      mutate(bidirectional = pmax(HvG, GvH, na.rm = TRUE))
    # Return the result based on the direction argument.
    if (direction == "HvG" | direction == "SOT") {
      return(MM_table$HvG)
    } else if (direction == "GvH") {
      return(MM_table$GvH)
    } else if (direction == "bidirectional") {
      return(MM_table$bidirectional)
    }
  } else {
    # Code to determine mismatch numbers if multiple loci were supplied.
    # Determine mismatches for both directions.
    HvG_table <- tibble("HvG" = HLA_mismatch_base(GL_string_recip, GL_string_donor, loci, "HvG", homozygous_count)) %>%
      # Add a row number to combine data at the end.
      mutate(case = row_number()) %>%
      # Separate the loci.
      separate_longer_delim(HvG, delim = ", ") %>%
      separate_wider_delim(HvG, delim = "=", names = c("locus", "mismatches")) %>%
      # Recode NA values to ensure accurate matching.
      mutate(mismatches = na_if(mismatches, "NA")) %>%
      # Count number of mismatches.
      mutate(HvG_number = replace_na(str_count(mismatches, "(\\+|$)"), 0)) %>%
      # Clean up table.
      select(-mismatches)

    GvH_table <- tibble("GvH" = HLA_mismatch_base(GL_string_recip, GL_string_donor, loci, "GvH", homozygous_count)) %>%
      # Add a row number to combine data at the end.
      mutate(case = row_number()) %>%
      # Separate the loci.
      separate_longer_delim(GvH, delim = ", ") %>%
      separate_wider_delim(GvH, delim = "=", names = c("locus", "mismatches")) %>%
      # Recode NA values to ensure accurate matching.
      mutate(mismatches = na_if(mismatches, "NA")) %>%
      # Count number of mismatches.
      mutate(GvH_number = replace_na(str_count(mismatches, "(\\+|$)"), 0)) %>%
      # Clean up table.
      select(-mismatches)

    # Join the GvH and HvG tables
    MM_table <- HvG_table %>% left_join(GvH_table, join_by(locus, case)) %>%
      # Calculate bidirectional mismatch number.
      mutate(bidirectional = pmax(HvG_number, GvH_number, na.rm = TRUE))

    # Return appropriate direction.
    # HvG
    if (direction == "HvG" | direction == "SOT") {
      MM_table <- MM_table %>%
        select(locus, case, HvG_number) %>%
        unite(locus, HvG_number, col = "MM", sep = "=") %>%
        summarise(MM = str_flatten(MM, collapse = ", "), .by = case)
      # GvH
    } else if (direction == "GvH") {
      MM_table <- MM_table %>%
        select(locus, case, GvH_number) %>%
        unite(locus, GvH_number, col = "MM", sep = "=") %>%
        summarise(MM = str_flatten(MM, collapse = ", "), .by = case)
      # Bidirectional
    } else if (direction == "bidirectional") {
      MM_table <- MM_table %>%
        select(locus, case, bidirectional) %>%
        unite(locus, bidirectional, col = "MM", sep = "=") %>%
        summarise(MM = str_flatten(MM, collapse = ", "), .by = case)
    }
    return(MM_table$MM)
  }
}

globalVariables(c(
  "mismatches", "case", "HvG_number", "GvH_number", "MM", "bidirectional",
  "str_count", "left_join", "join_by"
))
