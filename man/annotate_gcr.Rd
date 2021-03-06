% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/annotate_gcr.R
\name{annotate_gcr}
\alias{annotate_gcr}
\title{Add positional annotation to MSstats group comparison result}
\usage{
annotate_gcr(name, MSstats.GC, label)
}
\arguments{
\item{MSstats.GC}{MSstats group comparison result as generated by groupComparison()}

\item{label}{Label (character) of the comparision}

\item{name}{file name (character) of a fasta file in working directory}
}
\value{
This function returns a table (data.frame) in long format. For none proteotypic peptides this function returns multiple rows, each listing the peptide in a
different positional context. This ensures that all possible peptide origins can be considered in
the downstream analysis.
}
\description{
This function takes a peptide-level group comparison result as output by
MSstats::groupComparison and adds positional annotation to each peptide feature.
}
\details{
This function maps all peptide features to a reference proteome using their stripped aa sequence and perfect string matching.
In addition, it determines the N-terminal modification status of each feature.

Important prerequisites:

This functions expects the results of a MSstats peptide-level significance analysis as input,
since positional proteomics is a completely peptide centric approach.

Peptide-level significance analysis can be enforced easily by copying the content of PeptideSequence
to ProteinName prior to calling MSstats::dataProcess().
}
\examples{
annotate_gcr("uniprot-proteome_UP000000589.fasta", sample.data, "control vs. GluC")
}
