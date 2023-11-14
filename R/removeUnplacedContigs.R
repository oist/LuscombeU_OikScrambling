#' Remove unplaced contigs
#'
#' Discard all rages whose sequence names do not match _chr1_, _Chr1_, _chr2_,
#' _Chr2_, _PAR_, _XSR_ or _YSR_.
#'
#' @param gb A [`GBreaks`] object.
#' @param target Filter _target_ ranges.
#' @param query Filter _query_ ranges.
#'
#' @export

removeUplacedCongits <- function (gb, target = TRUE, query = TRUE) {
  chromosomeScaleScaffolds <- c('chr1', 'Chr1', 'chr2', 'Chr2', 'PAR', 'XSR', 'YSR')
  if (isTRUE(target))
    gb <- plyranges::filter(gb, seqnames %in% chromosomeScaleScaffolds)
  if (isTRUE(query))
    gb <- plyranges::filter(gb, seqnames(query) %in% chromosomeScaleScaffolds)
  gb
}
