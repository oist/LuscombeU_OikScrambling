#' Load all transcript information as GRanges
#'
#' Accessory function to load all transcript annotations provided by the
#' `BreakpointsData` package
#'
#' @returns Returns a [`SimpleList`] of `GRanges` objects representing the
#' annotations of the genomes we study.

loadAllTranscriptsGR <- function() {
  genomes <- OikScrambling:::loadAllGenomes()
  annots <- OikScrambling:::loadAllAnnotations() |> suppressWarnings()
  transcripts <- sapply(annots, GenomicFeatures::transcripts) |> SimpleList()
  # Start modifying the object here.
  transcripts
}
