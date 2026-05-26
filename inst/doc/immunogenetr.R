## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(
  collapse = FALSE,
  comment = "#>"
)

# Render data frames and tibbles as HTML tables in chunk output, instead of
# the default monospace console-style print. The kable HTML is wrapped in a
# pandoc raw-HTML block (```{=html} ... ```) so pandoc passes the content
# through verbatim. Without this, pandoc's markdown_in_html_blocks extension
# re-parses the cell contents as markdown, turning "^...^" in GL strings
# into <sup>...</sup>.
knit_print.data.frame <- function(x, ...) {
  html <- as.character(knitr::kable(x, format = "html"))
  knitr::asis_output(paste0("\n```{=html}\n", html, "\n```\n"))
}
registerS3method("knit_print", "data.frame", knit_print.data.frame, envir = asNamespace("knitr"))
registerS3method("knit_print", "tbl_df",     knit_print.data.frame, envir = asNamespace("knitr"))
registerS3method("knit_print", "tbl",        knit_print.data.frame, envir = asNamespace("knitr"))

## ----message=FALSE------------------------------------------------------------
library(immunogenetr)
library(dplyr)

## -----------------------------------------------------------------------------
# HLA_typing_1 contains typing for 10 individuals across all classical HLA loci.
head(HLA_typing_1, 3)

## -----------------------------------------------------------------------------
HLA_typing_GL <- HLA_typing_1 %>%
  # Convert all typing columns (A1 through DPB1_2) into a GL string.
  mutate(
    GL_string = HLA_columns_to_GLstring(., HLA_typing_columns = A1:DPB1_2),
    .after = patient
  ) %>%
  # Keep only patient ID and the new GL string column.
  select(patient, GL_string)

# View the GL strings.
(HLA_typing_GL)

## -----------------------------------------------------------------------------
# Take the first patient's GL string and split it into locus columns.
# Note: GLstring_genes and GLstring_genes_expanded use pivot_longer on all
# columns, so only pass the GL string column (no other data types).
single_patient <- HLA_typing_GL[1, "GL_string", drop = FALSE]
GLstring_genes(single_patient, "GL_string")

## -----------------------------------------------------------------------------
GLstring_genes_expanded(single_patient, "GL_string")

## -----------------------------------------------------------------------------
# Patient 7 is the recipient, patient 9 is the donor.
recip_gl <- HLA_typing_GL %>% filter(patient == 7) %>% pull(GL_string)
donor_gl <- HLA_typing_GL %>% filter(patient == 9) %>% pull(GL_string)

## -----------------------------------------------------------------------------
# Check if there is an HLA-A mismatch in the graft-vs-host direction.
HLA_mismatch_logical(recip_gl, donor_gl, "HLA-A", direction = "GvH")

# Check host-vs-graft direction.
HLA_mismatch_logical(recip_gl, donor_gl, "HLA-A", direction = "HvG")

## -----------------------------------------------------------------------------
# Count bidirectional mismatches across several loci at once.
HLA_mismatch_number(
  recip_gl, donor_gl,
  c("HLA-A", "HLA-B", "HLA-C", "HLA-DRB1"),
  direction = "bidirectional"
)

## -----------------------------------------------------------------------------
# Identify the specific mismatched alleles in the HvG direction.
HLA_mismatched_alleles(recip_gl, donor_gl, "HLA-A", direction = "HvG")

## -----------------------------------------------------------------------------
# Count the number of matches (complement of mismatches).
HLA_match_number(
  recip_gl, donor_gl,
  c("HLA-A", "HLA-B", "HLA-C", "HLA-DRB1"),
  direction = "bidirectional"
)

## -----------------------------------------------------------------------------
# X-of-8 matching (A, B, C, DRB1 bidirectional).
HLA_match_summary_HCT(recip_gl, donor_gl,
  direction = "bidirectional",
  match_grade = "Xof8"
)

# X-of-10 matching (adds DQB1).
HLA_match_summary_HCT(recip_gl, donor_gl,
  direction = "bidirectional",
  match_grade = "Xof10"
)

## -----------------------------------------------------------------------------
# Patient 3 is the recipient; compare against all 10 donors.
recipient <- HLA_typing_GL %>%
  filter(patient == 3) %>%
  select(GL_string) %>%
  rename(GL_string_recip = GL_string)

donors <- HLA_typing_GL %>%
  rename(GL_string_donor = GL_string, donor = patient) %>%
  # Cross-join to pair recipient with each donor.
  cross_join(recipient) %>%
  # Calculate 8/8 match grade for each pair.
  mutate(
    match_8of8 = HLA_match_summary_HCT(
      GL_string_recip, GL_string_donor,
      direction = "bidirectional",
      match_grade = "Xof8"
    ),
    .after = donor
  ) %>%
  # Sort best matches first.
  arrange(desc(match_8of8))

donors %>% select(donor, match_8of8)

## -----------------------------------------------------------------------------
# Truncate a four-field allele to two fields.
HLA_truncate("HLA-A*02:01:01:01", fields = 2)

# Works on full GL strings too.
HLA_truncate("HLA-A*02:01:01:01+HLA-A*03:01:01:02^HLA-B*07:02:01:01+HLA-B*44:02:01:01",
  fields = 2
)

## -----------------------------------------------------------------------------
# Remove all prefixes to get just the allele fields.
HLA_prefix_remove("HLA-A*02:01")

# Keep the locus designation but remove "HLA-".
HLA_prefix_remove("HLA-A*02:01", keep_locus = TRUE)

# Add the full prefix back.
HLA_prefix_add("02:01", "HLA-A*")

# "HLA-" is added by default.
HLA_prefix_add("A*02:01")

## -----------------------------------------------------------------------------
gl <- "HLA-A*02:01:01+HLA-A*68:01^HLA-B*07:01+HLA-B*15:01"

# A two-field search correctly matches the three-field allele.
pattern <- GLstring_regex("HLA-A*02:01")
stringr::str_detect(gl, pattern)

# But won't falsely match a longer allele number.
stringr::str_detect("HLA-A*02:149:01", GLstring_regex("HLA-A*02:14"))

## -----------------------------------------------------------------------------
# GLstring_genes returns tidyverse-friendly names by default.
repaired <- GLstring_genes(single_patient, "GL_string")
names(repaired)

# Convert back to WHO format with asterisks.
who_names <- HLA_column_repair(repaired, format = "WHO", asterisk = TRUE)
names(who_names)

## -----------------------------------------------------------------------------
# immunogenetr ships with two example HML files.
hml_path <- system.file("extdata", "HML_1.hml", package = "immunogenetr")
hml_result <- read_HML(hml_path)
hml_result

