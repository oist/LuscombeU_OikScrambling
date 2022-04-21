#' Calculate operon ranges by distance and colinearity
#'
#' Takes a gene annotation and computes the coordinates of operons based on
#' colinearity (all genes on the same strand with no interruption by
#' opposite-strand gens, and distance (all genes within a maximal distance of
#' each other).
#'
#' @param genes A `GRanges` object representing gene annotations.
#' @param window The maximum distance separating two genes of the same operon.
#' It has to be an even number.
#'
#' @note The `window` parameter must be an even number because the algorithm
#' expands the `genes` ranges by hald the window size on each side.  If an
#' odd number is provided, no error message is produced, but the results will be
#' the same as if the window were 1 base greater.
#'
#' @return A `GRanges` object representing operons, with a metadata column `n`
#' reporting the number of genes in each operon.
#'
#' @examples
#' genes <- GRanges(c("chr1:100-199:+", "chr1:300-400:+",
#'                    "chr1:500-600:+", "chr1:700-800:-"))
#' calcOperons(genes, window = 100)

calcOperons <- function(genes, window = 100) {
  # Divide the window by two because we will expand both sides of the range.
  halfwin <- round(window / 2)
  # Sort the genes ignoring strand.
  g <- sort(genes, ignore.strand = TRUE)
  # Expand by `halfwin` nucleotides
  # We suppress warnings because coordinates at the edges of seqfeatures can
  # become transiently negative.
  g <- (g + halfwin) |> suppressWarnings()
  # Give a unique ID to all runs of successive genes on the same strand
  # https://stackoverflow.com/questions/21421047
  # op$idx <- data.table::rleid(as.character(strand(g)))
  g$idx <- Rle(rep(seq_along(runValue(strand(g))), runLength(strand(g))))
  # Split the genes by run ID
  opL <- split(g, g$idx)
  # Remove all runs of length one
  opL <- opL[runLength(strand(g)) > 1]
  # Merge into operons
  operons   <- unlist(endoapply(opL, reduce, min.gapwidth = 0L)) |> suppressWarnings()
  # Removed un-merged genes
  operons$n <- countOverlaps(operons, genes)
  operons   <- operons[operons$n > 1]
  # Remove the flanking sequences
  operons - halfwin
}
