.make.ranges <- function(x){

## helper function for split+apply+combine

  r <- IRanges::IRanges(start = x$offset, width = nchar(x$s.pep), names = x$ID)
  S4Vectors::elementMetadata(r) <- x[c("log2FC", "adj.pvalue", "s.pep", "m.pep", "nterm.label", "first.aa", "last.aa", "up.aa", "down.aa")]
  return(r)

}

#' Transform annotated peptide data to list of IRanges
#'
#' @param df A peptide-level group comparision result (data.frame) annotated by posprot::annotate_gcr()
#' @param n The name of RangesList object returned by this function
#'
#' @return This function returns a RangesList object
#' @export
#'
#' @examples
df2RangesList <- function(df, n = "name") {
  if (anyNA(df$prot)) {

    warning("Detected cases with missing protein annotation! These will be removed.")
    df <- df[!is.na(df$prot),]

  }
  l <- plyr::dlply(.data = df, .variables = "prot", .fun = plyr::failwith(NA, .make.ranges), .progress = "text")
  rl <- as(l, "RangesList")
  #rl <- as(l, "IRangesList") # not implemented yet
  IRanges::universe(rl) <- n
  return(rl)
}
