#' Load all transcript information as GRanges
#'
#' Accessory function to load all transcript annotations provided by the
#' `BreakpointsData` package
#'
#' @include addMetadataToGRanges.R
#' @returns Returns a [`SimpleList`] of `GRanges` objects representing the
#' annotations of the genomes we study.

loadAllTranscriptsGR <- function() {
  genomes <- OikScrambling:::loadAllGenomes()
  annots <- OikScrambling:::loadAllAnnotations() |> suppressWarnings()
  transcripts <- sapply(annots, GenomicFeatures::transcripts) |> SimpleList()

  # Start modifying the object here.
  dnds <- SimpleList()
  dnds$dnds_guidance   <- loaddNdSTable("v2.0.0/genewise_oikopleura.extended.txt", alignment_type = 'GUIDANCE2')
  dnds$dnds_hmmcleaner <- loaddNdSTable("v2.0.0/genewise_oikopleura.extended.txt", alignment_type = 'HmmCleaner')
  dnds$dnds_prank      <- loaddNdSTable("v2.0.0/genewise_oikopleura.extended.txt", alignment_type = 'Alignment')
  transcripts <- sapply(names(transcripts), function(sp){
    sp_dnds = addMetadataToGRanges(gr=transcripts[[sp]], df=dnds$dnds_guidance[[sp]], gr_key='tx_name', df_key='transcript_id_simple', df_col_name = 'dNdS', rename_meta = 'dNdS_GUIDANCE2')
    sp_dnds = addMetadataToGRanges(gr=sp_dnds, df=dnds$dnds_hmmcleaner[[sp]], gr_key='tx_name', df_key='transcript_id_simple', df_col_name = 'dNdS', rename_meta = 'dNdS_HmmCleaner')
    sp_dnds = addMetadataToGRanges(gr=sp_dnds, df=dnds$dnds_prank[[sp]], gr_key='tx_name', df_key='transcript_id_simple', df_col_name = 'dNdS', rename_meta = 'dNdS_PRANK')
    sp_dnds
  })

  transcripts |> SimpleList()
}
