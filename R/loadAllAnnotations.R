#' Load all annotations
#'
#' Accessory function to load all annotations provided by the `BreakpointsData`
#' package
#'
#' @returns Returns a [`SimpleList`] of `TxDb` objects representing the
#' annotations of the genomes we study.

loadAllAnnotations <- function() {
  `%in%` <- BiocGenerics::`%in%`
  annots <- SimpleList()

  gff2txdb <- function(file, genome) {
    file <- system.file(paste0("extdata/Annotations/", file), package = "BreakpointsData")
    tx <- rtracklayer::import.gff(file)
    tx <- tx[seqnames(tx) %in% seqnames(genome)]
    tx <- GRanges(tx, seqinfo = seqinfo(genome))
    tx <- GenomicFeatures::makeTxDbFromGRanges(tx)
  }

  gff2txdb_Norway <- function(file, genome) {
    file <- system.file(paste0("extdata/Annotations/", file), package = "BreakpointsData")
    tx <- rtracklayer::import.gff(file)
    tx <- tx[!is.na(tx$mRNA)]
    # Remove "scaffoldA" objects in the OdB3 annotation
    tx <- tx[seqnames(tx) %in% seqnames(genome)]
    # Prepare gene entries
    tx$ID <- tx$mRNA
    tx$Parent <- ""
    tx.genes <- tx
    tx.genes$type <- "gene"
    # Prepare transcript IDs
    tx$type <- "transcript"
    tx$Parent <- tx$ID
    tx$ID <- paste0(tx$ID, ".t1")
    tx <- GRanges(c(tx, tx.genes), seqinfo = seqinfo(genome))
    tx <- GenomicFeatures::makeTxDbFromGRanges(tx)
  }

  gff2txdb_other <- function(file) {
    file <- system.file(paste0("extdata/Annotations/", file), package = "BreakpointsData")
    tx <- rtracklayer::import.gff(file)
    tx <- tx[!tx$type %in% c("stop_codon", "start_codon")]  # These cause many entries to be discarded by makeTxDbFromGRanges
    tx <- GenomicFeatures::makeTxDbFromGRanges(tx)
  }

  annots$Oki <- gff2txdb("OKI2018_I69.v2/OKI2018_I69.v2.gm.gff.gz",  BSgenome.Oidioi.OIST.OKI2018.I69)
  annots$Osa <- gff2txdb("OSKA2016v1.9/OSKA2016v1.9.gm.gff.gz",      BSgenome.Oidioi.OIST.OSKA2016v1.9)
  annots$Bar <- gff2txdb("Bar2_p4.Flye/Bar2_p4.Flye.gm.gff.gz",      BSgenome.Oidioi.OIST.Bar2.p4)
  annots$Kum <- gff2txdb("KUM-M3-7f/KUM-M3-7f.gm.gff.gz",            BSgenome.Oidioi.OIST.KUM.M3.7f)
  annots$Aom <- gff2txdb("AOM-5-5f/AOM-5-5f.gm.gff.gz",              BSgenome.Oidioi.OIST.AOM.5.5f)
  annots$Nor <- gff2txdb_Norway("OdB3/Oikopleura_annot_v1.0.gff.gz", BSgenome.Oidioi.genoscope.OdB3)
  annots$Ply <- gff2txdb_other("C_int_P/C_int_P.gm.gff.gz")
  annots$Ros <- gff2txdb_other("C_int_R/C_int_R.gm.gff.gz")

  annots
}
