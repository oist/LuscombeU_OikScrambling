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

a <- read_dnds('../oik_dnds/genes.extended.txt')

annotateWithDNDS.GRanges <- function(gr, dnds) {
    
}
}