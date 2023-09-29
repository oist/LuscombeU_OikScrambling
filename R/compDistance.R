#' Transform pair name into distance name.
#'
#' We have 3 main evolutionary distances: _same population_, _North Atlantic vs
#' North Pacific_, and _Okinawa vs the other populations_.  To ease plotting
#' this function transforms the names of pairwise comparisons in one of these
#' three classes.
#'
#' @param x A vector of pair names such as `Oki_Kum`, `Osa_Nor`, etc
#' @param short Further group all Ciona and Drosophila together.
#'
#' @return A vector of "`same_pop"`, `"Oki – North"` and `"North – North"`.

compDistance <- function(x, short = FALSE) {
  x[x %in% c("Oki_Kum", "Osa_Aom", "Bar_Nor")]                       <- "In same pop"
  x[x %in% c("Oki_Osa", "Oki_Bar", "Oki_Kum", "Oki_Aom", "Oki_Nor",
             "Osa_Oki", "Osa_Kum", "Bar_Oki", "Bar_Kum")]            <- "Oki – North"
  x[x %in% c("Osa_Bar", "Osa_Aom", "Osa_Nor",
             "Bar_Osa", "Bar_Aom", "Bar_Nor")]                       <- "North – North"
  x[x %in% c("Rob_Sav", "Ply_Sav")]                                  <- "Int/Rob – Sav"
  x[x %in% c("Rob_Oki", "Ply_Oki")]                                  <- "Int/Rob – Oki"
  x[x %in% c("Ply_Ros")]                                             <- "Int – Int"
  x[x %in% c("Ply_Rob", "Rob_Ros", "Rob_Ply")]                       <- "Int – Rob"
  if (isTRUE(short)) {
    x[x %in% c("Int/Rob – Sav", "Int/Rob – Oki", "Int – Int", "Int – Rob")] <- "Ciona"
    x[x %in% c("Dme_Dya", "Dme_Dma", "Dme_Dsu", "Dme_Dbu")] <-                 "Drosophila"
    x[x %in% c("Nig_Nig", "Nig_Bri", "Nig_Rem", "Nig_Ino",
               "Bri_Bri", "Bri_Nig", "Bri_Rem", "Bri_Ele")] <-                 "Caenorhabditis"
  }
  x
}

#' Transform pair name into genus name.
#'
#' We compare genomes within 3 genera: _Oikopleura_, _Ciona_, and _Drosophila_.
#' To ease plotting this function transforms the names of pairwise comparisons
#' in one of these three names, or `NA` for comparisons between genera.
#'
#' The levels are sorted so that _Oikopleura_ is at the top in default `ggplot`
#' figures, followed by _Ciona_ as both are tunicates.
#'
#' @param x A vector of pair names such as `Oki_Kum`, `Osa_Nor`, etc
#'
#' @return A factor of `Oikopleura`, `Ciona`, and `Drosophila`, or `NA`.

compGenus <- function(x) {
  x[x %in% c("Oki_Kum", "Osa_Aom", "Bar_Nor",
             "Oki_Osa", "Oki_Bar", "Oki_Kum", "Oki_Aom", "Oki_Nor",
             "Osa_Oki", "Osa_Kum", "Bar_Oki", "Bar_Kum",
             "Osa_Bar", "Osa_Aom", "Osa_Nor",
             "Bar_Osa", "Bar_Aom", "Bar_Nor")]                       <- "Oikopleura"
  x[x %in% c("Ply_Sav", "Ply_Ros", "Ply_Rob",
             "Rob_Sav", "Rob_Ros", "Rob_Ply")]                       <- "Ciona"
  x[x %in% c("Rob_Oki", "Ply_Oki")]                                  <- NA
  x[x %in% c("Dme_Dya", "Dme_Dma", "Dme_Dsu", "Dme_Dbu")]            <- "Drosophila"
  x[x %in% c("Nig_Nig", "Nig_Bri", "Nig_Rem", "Nig_Ino",
             "Bri_Bri", "Bri_Nig", "Bri_Rem", "Bri_Ele")]            <- "Caenorhabditis"

  factor(x, levels = c("Caenorhabditis", "Drosophila", "Ciona", "Oikopleura"))
}

#' Transform pair name into distance class.
#'
#' We have 4 main evolutionary distances: _same population_, _close_,
#' _intermediate_, _distant_.  To ease plotting this function transforms the
#' names of pairwise comparisons in one of these three classes.
#'
#' @param x A vector of pair names such as `Oki_Kum`, `Osa_Nor`, etc
#'
#' @return A vector of "`same_pop"`, `"close"`, `"intermediate"`, and `distant`.

compDistClass <- function(x) {
  x[x %in% c("Oki_Kum", "Osa_Aom", "Bar_Nor")]                       <- "same_or_sister"
  x[x %in% c("Oki_Osa", "Oki_Bar", "Oki_Kum", "Oki_Aom", "Oki_Nor",
             "Osa_Oki", "Osa_Kum", "Bar_Oki", "Bar_Kum")]            <- "distant"
  x[x %in% c("Osa_Bar", "Osa_Aom", "Osa_Nor",
             "Bar_Osa", "Bar_Aom", "Bar_Nor")]                       <- "intermediate"
  x[x %in% c("Rob_Sav", "Ply_Sav")]                                  <- "distant"
  x[x %in% c("Ply_Rob", "Rob_Ros", "Rob_Ply")]                       <- "close"
  x[x %in% c("Rob_Oki", "Ply_Oki")]                                  <- "different_genus"
  x[x %in% c("Ply_Ros")]                                             <- "same_or_sister"
  x[x %in% c("Dme_Dma", "Nig_Nig", "Bri_Bri")]                       <- "same_or_sister"
  x[x %in% c("Dme_Dya", "Nig_Bri", "Bri_Nig")]                       <- "close"
  x[x %in% c("Dme_Dsu", "Nig_Rem", "Bri_Rem")]                       <- "intermediate"
  x[x %in% c("Dme_Dbu", "Nig_Ino", "Bri_Ele")]                       <- "distant"
  factor(x, levels = c("same_or_sister", "close", "intermediate", "distant", "different_genus"))
}
