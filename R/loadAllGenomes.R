#' Load all genomes
#'
#' Utility function to load all BSgenome packages in the active session and
#' reduce boilerplate in the vignetters.
#'
#' @returns Returns a [`SimpleList`] of all loaded genomes.

loadAllGenomes <- function() {
  library("BSgenome.Odioica.local.OKI2018.I69")
  library("BSgenome.Odioica.local.OSKA2016v1.9")
  library("BSgenome.Odioica.local.Bar2.p4")
  library("BSgenome.Odioica.local.KUM.M3")
  library("BSgenome.Odioica.local.AOM.5")
  library("BSgenome.Odioica.local.Odioica.reference.v3.0")
  S4Vectors::SimpleList(
    # Chromosome assemblies
    Oki = OKI2018_I69, Osa = OSKA2016v1.9, Bar = Bar2_p4,
    # Less contiguous assemblies
    Kum = KUM_M3, Aom = AOM_5, Nor = OdB3)
}
