#' Find paired neo C-terminal neo N-terminal peptides
#'
#' @param x (IRanges S4 object)
#'
#' @param united (logical)
#'
#' @return Function returns an IRanges object by default. If united is set to FALSE a pairs object
#' returned.
#'
#' @description This function looks for adjacent peptides (peptides being positioned head-to-tail at zero distance)
#' taking the specificity of the test protease and the working protease into account.
#' The function also checks for the presence of N-terminal labels on the peptides. These indicate that a peptide lies downstream
#' of a cleavage site. Pairs of adjacent peptides, passing all specificity checks, are unified
#' and returned as a single IRange or Pairs object depending on united.
#'
#'
#' @details ---XXXXX.XXXXX--- ==> ---XXXXXXXXXX---
#'
#' This function expects test.prot.spec and work.prot.spec to be present in the environment.
#'
#' @export
#'
#' @examples
neoCN <- function(x, united = TRUE) {

  stopifnot(class(x) == "IRanges", class(united) == "logical")

  m <- S4Vectors::elementMetadata(x)
  if (!any(is.na(test.prot.spec))) {
    a <- x[m$nterm.label == "label" & m$last.aa %in% work.prot.spec & m$up.aa %in% test.prot.spec]
    b <- x[m$nterm.label == "free" & m$last.aa %in% test.prot.spec & m$up.aa %in% work.prot.spec]
  } else {
    a <- x[m$nterm.label == "label" & m$last.aa %in% work.prot.spec]
    b <- x[m$nterm.label == "free" & m$up.aa %in% work.prot.spec]
  }
  p <- IRanges::findOverlapPairs(a, b, maxgap=1L)
  if (united) {
    IRanges::punion(BiocGenerics::subset(p, start(first) == end(second) + 1L))
  } else {
    BiocGenerics::subset(p, start(first) == end(second) + 1L)
  }
}

#' Find paired neo N- and neo C-terminal peptides over a list of proteins
#'
#' @param l is a list of IRanges representing peptides matched to proteins
#'
#' @return By default this function returns a list of IRanges representing all paired peptides found in the input data.
#' In case united = FALSE is provided as an argument, it returns a list of Pairs objects.
#'
#' @description
#'
#' @details ---XXXXX.XXXXX--- ==> ---XXXXXXXXXX---
#'
#' @export
#'
#' @examples
find_neoCN_pairs <- function(l, ...) {

  stopifnot(class(l) == "list",
            class(work.prot.spec) == "character",
            class(test.prot.spec) == "character"
  )
  message(paste("Working protease specificity:", paste(work.prot.spec, collapse = "")))
  message(paste("Test protease specificity:", paste(test.prot.spec, collapse = "")))

  plyr::llply(.data = l, .fun = plyr::failwith(default = NA, f = neoCN), ..., .progress = "text")

}
