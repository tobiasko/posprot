#' Retrieve protease specificity by name
#'
#' @param protease Name (character) of the protease
#'
#' @return A character vector containing amino acids cleaved by the protease
#' @export
#'
#' @examples spec("GluC")
spec <- function (protease = "Trypsin") {

  ## spec.table is in R/sysData

  message(paste("Retrieving "), protease, " specificity.")
  unlist(strsplit(x = as.character(spec.table[match(protease, spec.table$name), "aa"]), split = ""))

}
