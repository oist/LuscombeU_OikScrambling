#' Calculate operon ranges by distance and colinearity
#'
#' Takes a gene annotation and computes the coordinates of operons based on
#' colinearity (all genes on the same strand with no interruption by
#' opposite-strand genes, and distance (all genes within a maximal distance of
#' each other).  Nested genes on the opposite strand do not interrupt
#' colinearity as they do not interrupt transcription.
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
#' reporting the number of genes in each operon.  If the input object had a
#' `gene_id` column, then the gene identifiers will be reported as a
#' `IRanges::CharacterList()` in a `gene_id` column in the operons object.
#'
#' @examples
#' genes <- GRanges(c("chr1:100-199:+", "chr1:300-400:+",
#'                    "chr1:500-600:+", "chr1:700-800:-"))
#' genes$gene_id <- LETTERS[seq_along(genes)]
#' OikScrambling:::calcOperons(genes, window = 100)

calcOperons <- function(genes, window = 100) {
  # Calculating operons on each strand separately so that nested genes do not
  # interrupt collinearity
  strands <- split(genes, strand(genes))
  ops.plus  <- calcOperons.onestrand(strands[['+']], window = window)
  ops.minus <- calcOperons.onestrand(strands[['-']], window = window)
  sort(c(ops.plus, ops.minus), ignore.strand = TRUE)
}

calcOperons.onestrand <- function(genes, window = 100) {
  if (identical(genes, GRanges())) return(GRanges())
  if(nrun(strand(genes)) != 1)
    stop("Do not use this function on GRanges objects with more than one strand")
  # Divide the window by two because we will expand both sides of the range.
  halfwin <- round(window / 2)
  # Sort the genes ignoring strand.
  g <- sort(genes, ignore.strand = TRUE)
  # Expand by `halfwin` nucleotides
  # We suppress warnings because coordinates at the edges of seqfeatures can
  # become transiently negative.
  g <- (g + halfwin) |> suppressWarnings()
  # Merge into operons
  operons   <- reduce(g, min.gapwidth = 0L)
  # Remove the flanking sequences
  operons <- operons - halfwin
  # Count number of genes in operons
  operons$n <- countOverlaps(operons, genes)
  # Annotate operons with gene names.
  if(!is.null(genes$gene_id)){
    ov <- findOverlaps(operons, genes)
    operons$gene_id <- CharacterList(split(genes$gene_id[subjectHits(ov)], queryHits(ov))) |> unname()
  }
  # Remove non-merged genes and return
  operons[operons$n > 1]
}
