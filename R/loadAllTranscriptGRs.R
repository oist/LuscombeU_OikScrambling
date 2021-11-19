#' Load all transcript information as GRanges
#'
#' Accessory function to load all transcript annotations provided by the
#' `BreakpointsData` package
#'
#' @author Charles Plessy
#' @author Michael Mansfield
#'
#' @returns Returns a [`SimpleList`] of `GRanges` objects representing the
#' annotations of the genomes we study.

loadAllTranscriptsGR <- function() {
  genomes <- OikScrambling:::loadAllGenomes()
  # Load TxDB annotations
  annots <- OikScrambling:::loadAllAnnotations() |> suppressWarnings()
  # Convert to GRanges on transcript features.
  transcripts <- sapply(annots, GenomicFeatures::transcripts) |> SimpleList()
  # Load dNdS values
  dnds <- SimpleList()
  dnds$dnds_guidance   <- loaddNdSTable("v2.0.0/genewise_oikopleura.extended.txt", alignment_type = 'GUIDANCE2')
  dnds$dnds_hmmcleaner <- loaddNdSTable("v2.0.0/genewise_oikopleura.extended.txt", alignment_type = 'HmmCleaner')
  dnds$dnds_prank      <- loaddNdSTable("v2.0.0/genewise_oikopleura.extended.txt", alignment_type = 'Alignment')
  # Stitch them in the objects
  transcripts <- sapply(names(transcripts), function(sp){
    sp_dnds = addMetadataToGRanges(gr=transcripts[[sp]], df=dnds$dnds_guidance[[sp]], gr_key='tx_name', df_key='transcript_id_simple', df_col_name = 'dNdS', rename_meta = 'dNdS_GUIDANCE2')
    sp_dnds = addMetadataToGRanges(gr=sp_dnds, df=dnds$dnds_hmmcleaner[[sp]], gr_key='tx_name', df_key='transcript_id_simple', df_col_name = 'dNdS', rename_meta = 'dNdS_HmmCleaner')
    sp_dnds = addMetadataToGRanges(gr=sp_dnds, df=dnds$dnds_prank[[sp]], gr_key='tx_name', df_key='transcript_id_simple', df_col_name = 'dNdS', rename_meta = 'dNdS_PRANK')
    sp_dnds
  })
  # Give transcript names to the rows
  transcripts <- sapply(transcripts, \(x) {names(x) <- x$tx_name ; x}) |> SimpleList()
  transcripts
}
