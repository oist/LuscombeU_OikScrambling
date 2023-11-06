#' Load all transcript information as GRanges
#'
#' Accessory function to load all transcript annotations provided by the
#' `BreakpointsData` package
#'
#' @author Charles Plessy
#' @author Michael Mansfield
#'
#' @param compat Set to TRUE to keep the old names.  This option will be removed
#' when the migration will be complete.
#'
#' @returns Returns a [`SimpleList`] of `GRanges` objects representing the
#' annotations of the genomes we study.

loadAllTranscriptsGR <- function(compat = TRUE) {
  genomes <- OikScrambling:::loadAllGenomes()
  # Load TxDB annotations
  annots <- OikScrambling:::loadAllAnnotations() |> suppressWarnings()
  # Convert to GRanges on transcript features.
  transcripts <- sapply(annots, GenomicFeatures::transcripts) |> SimpleList()
  transcripts$Nor$tx_name <- transcripts$Nor$tx_name |> sub(pat = ".t1", rep = "")
  # Load dNdS values
  dnds <- SimpleList()
  dnds$dnds_guidance   <- OikScrambling:::loaddNdSTable("v3.0.0/genewise_oikopleura.extended.txt", alignment_type = 'GUIDANCE2')
  dnds$dnds_hmmcleaner <- OikScrambling:::loaddNdSTable("v3.0.0/genewise_oikopleura.extended.txt", alignment_type = 'HmmCleaner')
  dnds$dnds_prank      <- OikScrambling:::loaddNdSTable("v3.0.0/genewise_oikopleura.extended.txt", alignment_type = 'Alignment')
  # Stitch them in the objects
  transcripts[1:6] <- sapply(names(transcripts[1:6]), function(sp){
    sp_dnds = OikScrambling:::addMetadataToGRanges(gr=transcripts[[sp]], df=dnds$dnds_guidance[[sp]], gr_key='tx_name', df_key='transcript_id_simple', df_col_name = 'dNdS', rename_meta = 'dNdS_GUIDANCE2')
    sp_dnds = OikScrambling:::addMetadataToGRanges(gr=sp_dnds, df=dnds$dnds_hmmcleaner[[sp]], gr_key='tx_name', df_key='transcript_id_simple', df_col_name = 'dNdS', rename_meta = 'dNdS_HmmCleaner')
    sp_dnds = OikScrambling:::addMetadataToGRanges(gr=sp_dnds, df=dnds$dnds_prank[[sp]], gr_key='tx_name', df_key='transcript_id_simple', df_col_name = 'dNdS', rename_meta = 'dNdS_PRANK')
    sp_dnds
  })
  # Give transcript names to the rows
  transcripts <- sapply(transcripts, \(x) {names(x) <- x$tx_name ; x}) |> SimpleList()

  if(isFALSE(compat))
    names(transcripts) <- c("OKI2018.I69", "OSKA2016v1.9", "Bar2.p4", "KUM.M3.7f", "AOM.5.5f", "OdB3", "Ply", "Ros")

  transcripts
}
