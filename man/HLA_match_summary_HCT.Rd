% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/HLA_match_summary_HCT.R
\name{HLA_match_summary_HCT}
\alias{HLA_match_summary_HCT}
\title{HLA_match_summary_HCT}
\usage{
HLA_match_summary_HCT(
  GL_string_recip,
  GL_string_donor,
  direction = "bidirectional",
  match_grade
)
}
\arguments{
\item{GL_string_recip}{A GL string representing the recipient's HLA genotype,
and minimally containing the HLA-A, B, C and DRB1 loci (for Xof8 matching)
or the HLA-A, B, C, DRB1 and DQB1 loci (for Xof10 matching).}

\item{GL_string_donor}{A GL string representing the donor's HLA genotype,
and minimally containing the HLA-A, B, C and DRB1 loci (for Xof8 matching)
or the HLA-A, B, C, DRB1 and DQB1 loci (for Xof10 matching).}

\item{direction}{"GvH", "HvG" or "bidirectional". Default is "bidirectional".}

\item{match_grade}{"Xof8" for HLA-A, B, C and DRB1 matching or "Xof10" for
HLA-A, B, C, DRB1 and DQB1 matching.}
}
\value{
An integer value of the match grade summary.
}
\description{
Calculates the match summary for either the HLA-A, B, C and DRB1
loci (out-of-8 matching) or the HLA-A, B, C, DRB1 and DQB1 loci (out-of-10 matching),
as is commonly used for hematopoietic cell transplantation (HCT). Homozygous
mismatches are counted twice. Bidirectional matching is the default, but can
be overridden with the "direction" argument.
}
\examples{
# Example recipient and donor GL strings
file <- HLA_typing_1[, -1]
GL_string <- HLA_columns_to_GLstring(file, HLA_typing_columns = everything())

GL_string_recip <- GL_string[1]
GL_string_donor <- GL_string[2]

# Calculate mismatch numbers
HLA_match_summary_HCT(GL_string_recip, GL_string_donor,
  direction = "bidirectional", match_grade = "Xof8"
)

}
