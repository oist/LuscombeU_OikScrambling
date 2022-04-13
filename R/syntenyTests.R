#' Flag pairs that are on homologous chromosomes
#'
#' This function assumes that homologous chromosomes have the same name, but
#' allows differences in character case.
#'
#' @param gb A [`GenomicBreaks::GBreaks`] object.
#'
#' @return A logical vector of same length as the `gb` object

isSyntenic <- function(gb) {
  tolower(as.character(seqnames(gb))) == tolower(as.character(seqnames(gb$query)))
}

#' Reports syntenic chromosome name or NA
#'
#' This function converts chromosome names to lower case.
#'
#' @param gb A [`GenomicBreaks::GBreaks`] object.
#'
#' @return A character vector of same length as the `gb` object

nameSyntenic <- function(gb) {
  ifelse(isSyntenic(gb), tolower(as.character(seqnames(gb))), NA)
}

#' Flag pairs that are on homologous chromosome arms
#'
#' This function assumes that homologous chromosomes have the same name, but
#' allows differences in character case.
#'
#' @param gb A [`GenomicBreaks::GBreaks`] object.
#'
#' @return A logical vector of same length as the `gb` object

isSynbrachial <- function(gb) {
  paste(tolower(as.character(seqnames(gb))), gb$Arm) == paste(tolower(as.character(seqnames(gb$query))), gb$query$Arm)
}

#' Reports syntenic chromosome arm name or NA
#'
#' This function converts chromosome names to lower case.
#'
#' @param gb A [`GenomicBreaks::GBreaks`] object.
#'
#' @return A character vector of same length as the `gb` object

nameSynbrachial <- function(gb) {
  ifelse(isSynbrachial(gb),  paste(tolower(as.character(seqnames(gb))), gb$Arm), NA)
}
