#' Load all repeat annotations
#'
#' Load all the repeat annotations distributed by the `BreakpointsData` data
#' package, and sort them in broad classes by parsing their `Target` attribute.
#'
#' @returns Returns a `SimpleList` of all loaded repeat annotations.

loadAllRepeats <- function() {
  reps <- SimpleList()
  # Suppressing warning message: ‘In is.na(genome) : is.na() applied to non-(list or vector) of type 'S4'’
  reps$Oki <- rtracklayer::import(system.file("extdata/Annotations/OKI2018_I69.v2/OKI2018_I69.repeats.gff", package = "BreakpointsData"),
                                  genome = seqinfo(BSgenome.Odioica.local.OKI2018.I69)) |> suppressWarnings()
  reps$Osa <- rtracklayer::import(system.file("extdata/Annotations/OSKA2016v1.9/OSKA2016v1.9.repeats.gff",  package = "BreakpointsData"),
                                  genome = seqinfo(BSgenome.Odioica.local.OSKA2016v1.9)) |> suppressWarnings()
  reps$Bar <- rtracklayer::import(system.file("extdata/Annotations/Bar2_p4.Flye/Bar2_p4.Flye.repeats.gff",  package = "BreakpointsData"),
                                  genome = seqinfo(BSgenome.Odioica.local.Bar2.p4)) |> suppressWarnings()
  reps$Kum <- rtracklayer::import(system.file("extdata/Annotations/KUM-M3-7f/KUM-M3-7f.repeats.gff",        package = "BreakpointsData"),
                                  genome = seqinfo(BSgenome.Odioica.local.KUM.M3)) |> suppressWarnings()
  reps$Aom <- rtracklayer::import(system.file("extdata/Annotations/AOM-5-5f/AOM-5-5f.repeats.gff",          package = "BreakpointsData"),
                                  genome = seqinfo(BSgenome.Odioica.local.AOM.5)) |> suppressWarnings()
  reps$Nor <- rtracklayer::import(system.file("extdata/Annotations/OdB3/OdB3.repeats.gff",                  package = "BreakpointsData"),
                                  genome = seqinfo(BSgenome.Odioica.local.Odioica.reference.v3.0)) |> suppressWarnings()
  reps <- endoapply(reps, \(x) {
    x$Class <- x$Target                                              |>
      sub(pat = '".*',                     rep = "")                 |>
      sub(pat = "Motif:",                  rep = "")                 |>
      sub(pat = "\\(.+)n",                 rep = "tandem")           |>
      sub(pat = "_family.*",               rep = "")                 |>
      sub(pat = "SINE.*",                  rep = "SINE")             |>
      sub(pat = "MITE.*",                  rep = "MITE")             |>
      sub(pat = "i69_juicer.*",            rep = "unknown")          |>
      sub(pat = "oska2016v1.9_.*",         rep = "unknown")          |>
      sub(pat = "bar2_p4_.*",              rep = "unknown")          |>
      sub(pat = "kum-m3-7f_.*",            rep = "unknown")          |>
      sub(pat = "aom-5-5f_.*",             rep = "unknown")          |>
      sub(pat = "Chr1",                    rep = "unknown")          |>
      sub(pat = "Chr2",                    rep = "unknown")          |>
      sub(pat = "PAR",                     rep = "unknown")          |>
      sub(pat = "S107",                    rep = "unknown")          |>
      sub(pat = "S144",                    rep = "unknown")          |>
      sub(pat = "S215",                    rep = "unknown")          |>
      sub(pat = "XSR",                     rep = "unknown")          |>
      sub(pat = "YSR",                     rep = "unknown")          |>
      sub(pat = "contig.*",                rep = "unknown")          |>
      sub(pat = "HiC_scaffold_.*",         rep = "unknown")          |>
      sub(pat = ".*-rich",                 rep = "LowComplexity")    |>
      sub(pat = "rnd.*",                   rep = "rnd")              |>
      factor()
    x
  })
  reps |> as("SimpleList")
}
