#' Transform pair name into distance class.
#'
#' We have 3 main evolutionary distances: _same population_, _North Atlantic vs
#' North Pacific_, and _Okinawa vs the other populations_.  To ease plotting
#' this function transforms the names of pairwise comparisons in one of these
#' three classes.
#'
#' @param x A vector of pair names such as `Oki_Kum`, `Osa_Nor`, etc
#'
#' @return A vector of "`same_pop"`, `"Oki – North"` and `"North – North"`.

compDistance <- function(x) {
  x[x %in% c("Oki_Kum", "Osa_Aom", "Bar_Nor")]                       <- "same_pop"
  x[x %in% c("Oki_Osa", "Oki_Bar", "Oki_Kum", "Oki_Aom", "Oki_Nor",
             "Osa_Oki", "Osa_Kum", "Bar_Oki", "Bar_Kum")]            <- "Oki – North"
  x[x %in% c("Osa_Bar", "Osa_Aom", "Osa_Nor",
             "Bar_Osa", "Bar_Aom", "Bar_Nor")]                       <- "North – North"
  x[x %in% c("Rob_Sav", "Ply_Sav")]                                  <- "Int/Rob – Sav"
  x[x %in% c("Rob_Oki", "Ply_Oki")]                                  <- "Int/Rob – Oki"
  x[x %in% c("Ply_Ros")]                                             <- "Int – Int"
  x[x %in% c("Ply_Rob", "Rob_Ros", "Rob_Ply")]                       <- "Int – Rob"
  x
}
