#' Load all genomes
#'
#' Utility function to load all BSgenome packages in the active session and
#' reduce boilerplate in the vignettes.
#'
#' @param compat Set to TRUE to keep the old names.  This option will be removed
#' when the migration will be complete.
#'
#' @returns Returns a [`SimpleList`] of all loaded genomes.

loadAllGenomes <- function(compat = TRUE) {
suppressPackageStartupMessages({
  library("BSgenome.Oidioi.OIST.OKI2018.I69")
  library("BSgenome.Oidioi.OIST.OSKA2016v1.9")
  library("BSgenome.Oidioi.OIST.Bar2.p4")
  library("BSgenome.Oidioi.OIST.KUM.M3.7f")
  library("BSgenome.Oidioi.OIST.AOM.5.5f")
  library("BSgenome.Oidioi.genoscope.OdB3")
})

genomes <-
  c("OKI2018.I69", "OSKA2016v1.9", "Bar2.p4", "KUM.M3.7f", "AOM.5.5f", "OdB3") |>
    sapply(BSgenome::getBSgenome) |> SimpleList()

if(isTRUE(compat))
  names(genomes) <- c("Oki", "Osa", "Bar", "Kum", "Aom", "Nor")

genomes
}
