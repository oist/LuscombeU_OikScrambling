#!/bin/Rscript

read_dnds <- function(file, annotation_list, alignment_type='GUIDANCE2') {
  # note that the 10th column is just NA. So I remove it.
  df = read.delim(file, header=F)[,-10]
  colnames(df) = c('HOG', 'dNdS', 'Alignment', 'OKI2018_I69.v2', 'KUM-M3-7f', 'OSKA2016v1.9', 'AOM-5-5f', 'Bar2_p4_Flye', 'OdB3')
  df = subset(df, df$Alignment == alignment_type)
  # Remove NAs
  df = df[!is.na(df$dNdS),]
  dnds = SimpleList()
  dnds$Oki = data.frame(gene_id=unname(sapply(df$"OKI2018_I69.v2", function(s) unlist(strsplit(s, '\\.'))[3])),  transcript_id=df$"OKI2018_I69.v2", dNdS=df$dNdS)
  dnds$Osa = data.frame(gene_id=unname(sapply(df$"OSKA2016v1.9", function(s) unlist(strsplit(unlist(strsplit(s, '\\.'))[3], '\\.t'))[1])),  transcript_id=df$"OSKA2016v1.9", dNdS=df$dNdS)
  dnds$Bar = data.frame(gene_id=unname(sapply( df$"Bar2_p4_Flye", function(s) unlist(strsplit(unlist(strsplit(s, '\\.'))[2], '\\.t'))[1])),  transcript_id=df$"Bar2_p4_Flye", dNdS=df$dNdS) # setNames(df$dNdS, df$"Bar2_p4_Flye") 
  dnds$Kum = data.frame(gene_id=unname(sapply( df$"KUM-M3-7f", function(s) unlist(strsplit(unlist(strsplit(s, '\\.'))[2], '\\.t'))[1])),  transcript_id=df$"KUM-M3-7f", dNdS=df$dNdS) # setNames(df$dNdS, df$"KUM-M3-7f") 
  dnds$Aom =  data.frame(gene_id=unname(sapply( df$"AOM-5-5f", function(s) unlist(strsplit(unlist(strsplit(s, '\\.'))[2], '\\.t'))[1])),  transcript_id=df$"AOM-5-5f", dNdS=df$dNdS) # setNames(df$dNdS, df$"AOM-5-5f") 
  dnds$Nor = data.frame(gene_id=unname(sapply( df$"OdB3", function(s) unlist(strsplit(unlist(strsplit(s, '\\.'))[2], '\\.t'))[1])),  transcript_id=df$"OdB3", dNdS=df$dNdS) # setNames(df$dNdS, df$"OdB3") 
  dnds |> as("SimpleList")
}

dnds_tb <- read_dnds('../oik_dnds/genes.extended.txt')

# Note that Norway is missing its annotations and so I exclude it.
dnds_granges = lapply(names(dnds_tb)[1:5], function(name) {
  dnds = dnds_tb[[name]]
  annot = genes(annots[[name]])
  gr = annot
  gr$dnds = NA
  gr$dnds = dnds[match(gr$gene_id, dnds$gene_id),]$dNdS
  gr = gr[!is.na(gr$dnds)]
  return(gr)
})
dnds_granges <- SimpleList(dnds_granges)
names(dnds_granges) <- names(dnds_tb)[1:5]

annotateWithdNdS.GRanges <- function(gr, dnds_gene) {
  # Inspired from CAGEr:::ranges2names
  o <- findOverlaps(gr, dnds_gene)
  o <- as(o, "List")
  o <- extractList(dnds_gene$dnds |> as.numeric(), o)
  o <- endoapply(o, unique)
  o[sapply(o, length) == 0] <- NA
  gr$dnds <- o
  gr
}

annotateWithdNdS.GBreaks <- function(gb, repT, repQ) {
  gb$query <- annotateWithdNdS.GRanges(gb$query, repQ)
  gb       <- annotateWithdNdS.GRanges(gb,       repT)
}

annotateWithdNdS.GBreaksSimpleList <- function(gbl, dnds) {
  gbl$Oki_Osa <- gbl$Oki_Osa |> annotateWithdNdS.GBreaks(dnds$Oki, dnds$Osa)
  gbl$Oki_Bar <- gbl$Oki_Bar |> annotateWithdNdS.GBreaks(dnds$Oki, dnds$Bar)
  gbl$Oki_Kum <- gbl$Oki_Kum |> annotateWithdNdS.GBreaks(dnds$Oki, dnds$Kum)
  gbl$Oki_Aom <- gbl$Oki_Aom |> annotateWithdNdS.GBreaks(dnds$Oki, dnds$Aom)

  gbl$Osa_Oki <- gbl$Osa_Oki |> annotateWithdNdS.GBreaks(dnds$Osa, dnds$Oki)
  gbl$Osa_Bar <- gbl$Osa_Bar |> annotateWithdNdS.GBreaks(dnds$Osa, dnds$Bar)
  gbl$Osa_Kum <- gbl$Osa_Kum |> annotateWithdNdS.GBreaks(dnds$Osa, dnds$Kum)
  gbl$Osa_Aom <- gbl$Osa_Aom |> annotateWithdNdS.GBreaks(dnds$Osa, dnds$Aom)

  gbl$Bar_Oki <- gbl$Bar_Oki |> annotateWithdNdS.GBreaks(dnds$Bar, dnds$Oki)
  gbl$Bar_Osa <- gbl$Bar_Osa |> annotateWithdNdS.GBreaks(dnds$Bar, dnds$Osa)
  gbl$Bar_Kum <- gbl$Bar_Kum |> annotateWithdNdS.GBreaks(dnds$Bar, dnds$Kum)
  gbl$Bar_Aom <- gbl$Bar_Aom |> annotateWithdNdS.GBreaks(dnds$Bar, dnds$Aom)
  gbl
}

gbs_dnds = annotateWithdNdS.GBreaksSimpleList(gbs, dnds_granges)
