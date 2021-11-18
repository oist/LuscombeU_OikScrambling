#' Load dNdS data from OikScrambling data
#'
#' Accessory function to load the dNdS table provided by the `BreakpointsData` package.
#' Requires `BreakpointsData` version >=3.6.1.
#'
#' @author Michael Mansfield
#'
#' @returns Returns a [`SimpleList`] of data frames with dN/dS annotations.
loaddNdSTable <- function(file, alignment_type='GUIDANCE2') {
  # The first part of the function loads the dN/dS data frame from the BreakpointsData package's extdata.
  get_dnds_df <- function(file, alignment_type){
    file <- system.file(paste0('extdata/dNdS/', file), package='BreakpointsData')
    # note that the 10th column is just NA. So I remove it. :-)
    df <- read.delim(file, header=F)#[,-10]
    df <- Filter(function(x)!all(is.na(x)), df)
    colnames(df) <- c('HOG', 'dNdS', 'Alignment', 'OKI2018_I69.v2', 'KUM-M3-7f', 'OSKA2016v1.9', 'AOM-5-5f', 'Bar2_p4_Flye', 'OdB3')
    df <- subset(df, df$Alignment == alignment_type)
    df
  }
  df <- get_dnds_df("v2.0.0/genewise_oikopleura.extended.txt", alignment_type=alignment_type)
  # Remove NAs just in case.
  df <- df[!is.na(df$dNdS),]

  # This section creates a SimpleList for every species containing a dataframe with
  # some identifiers parsed out as well as dN/dS values.
  dnds <- SimpleList()
  # Do note that the gene_id and transcript_id parsing is specific to each sample (i.e.,
  # the indices for parsing out [gene_id].[transcript_id] change).
  # This parsing could be improved, but, ah, this works.
  dnds$Oki <- data.frame(
    gene_id=unname(sapply(df$"OKI2018_I69.v2", function(s) unlist(strsplit(s, '\\.'))[3])),
    transcript_id_full=df$"OKI2018_I69.v2",
    transcript_id_simple=unname(sapply(df$"OKI2018_I69.v2", function(s) paste(unlist(strsplit(s, '\\.'))[3:4], collapse='.'))),
    dNdS=df$dNdS)
  dnds$Osa <- data.frame(
    gene_id=unname(sapply(df$"OSKA2016v1.9", function(s) unlist(strsplit(unlist(strsplit(s, '\\.'))[3], '\\.t'))[1])),
    transcript_id_full=df$"OSKA2016v1.9",
    transcript_id_simple=unname(sapply( df$"OSKA2016v1.9", function(s) paste(unlist(strsplit(s, '\\.'))[3:4], collapse='.'))),
    dNdS=df$dNdS)
  dnds$Bar <- data.frame(
      gene_id=unname(sapply(df$"Bar2_p4_Flye", function(s) unlist(strsplit(unlist(strsplit(s, '\\.'))[2], '\\.t'))[1])),
      transcript_id_full=df$"Bar2_p4_Flye",
      transcript_id_simple=unname(sapply( df$"Bar2_p4_Flye", function(s) paste(unlist(strsplit(s, '\\.'))[2:3], collapse='.'))),
      dNdS=df$dNdS)
  dnds$Kum <- data.frame(
    gene_id=unname(sapply(df$"KUM-M3-7f", function(s) unlist(strsplit(unlist(strsplit(s, '\\.'))[2], '\\.t'))[1])),
    transcript_id_full=df$"KUM-M3-7f",
    transcript_id_simple=unname(sapply( df$"KUM-M3-7f", function(s) paste(unlist(strsplit(s, '\\.'))[2:3], collapse='.'))),
    dNdS=df$dNdS)
  dnds$Aom <- data.frame(
    gene_id=unname(sapply(df$"AOM-5-5f", function(s) unlist(strsplit(unlist(strsplit(s, '\\.'))[2], '\\.t'))[1])),
    transcript_id_full=df$"AOM-5-5f",
    transcript_id_simple=unname(sapply( df$"AOM-5-5f", function(s) paste(unlist(strsplit(s, '\\.'))[2:3], collapse='.'))),
    dNdS=df$dNdS)
  dnds$Nor <- data.frame(
    gene_id=unname(sapply(df$"OdB3", function(s) unlist(strsplit(unlist(strsplit(s, '\\.'))[2], '\\.t'))[1])),
    transcript_id_full=df$"OdB3",
    transcript_id_simple=unname(sapply(df$"OdB3", function(s) paste(unlist(strsplit(s, '\\.'))[2], collapse='.'))),
    dNdS=df$dNdS)
  dnds |> as("SimpleList")
}
