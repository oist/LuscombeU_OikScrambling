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

#a <- read_dnds('../oik_dnds/genes.extended.txt')

loadAnnotationsWithDNDS <- function(dnds_list) {
  `%in%` <- BiocGenerics::`%in%`
  annots <- SimpleList()
  
  gff2txdb <- function(file, genome, dnds_meta) {
    file <- system.file(paste0("extdata/Annotations/", file), package = "BreakpointsData")
    tx <- rtracklayer::import.gff(file)
    tx <- tx[seqnames(tx) %in% seqnames(genome)]
    tx <- GRanges(tx, seqinfo = seqinfo(genome))
    tx$dnds = NA
    tx$dnds = dnds_meta[match(tx$ID, dnds_meta$gene_id),]$dNdS
    tx <- GenomicFeatures::makeTxDbFromGRanges(tx)
  }
  
  gff2txdb_Norway <- function(file, genome, dnds_meta) {
    file <- system.file(paste0("extdata/Annotations/", file), package = "BreakpointsData")
    tx <- rtracklayer::import.gff(file)
    tx <- tx[!is.na(tx$mRNA)]
    # Remove "scaffoldA" objects in the OdB3 annotation
    tx <- tx[seqnames(tx) %in% seqnames(genome)]
    tx$ID <- tx$mRNA
    tx <- GRanges(tx, seqinfo = seqinfo(genome))
    tx$dnds = NA
    tx$dnds = dnds_meta[match(tx$ID, dnds_meta$gene_id),]$dNdS
    tx <- GenomicFeatures::makeTxDbFromGRanges(tx)
  }
  
  annots$Oki <- gff2txdb("OKI2018_I69.v2/OKI2018_I69.v2.gm.gff.gz",  OKI2018_I69,  dnds_list$Oki)
  annots$Osa <- gff2txdb("OSKA2016v1.9/OSKA2016v1.9.gm.gff.gz",      OSKA2016v1.9, dnds_list$Osa)
  annots$Bar <- gff2txdb("Bar2_p4.Flye/Bar2_p4.Flye.gm.gff.gz",      Bar2_p4,      dnds_list$Bar)
  annots$Kum <- gff2txdb("KUM-M3-7f/KUM-M3-7f.gm.gff.gz",            KUM_M3,       dnds_list$Kum)
  annots$Aom <- gff2txdb("AOM-5-5f/AOM-5-5f.gm.gff.gz",              AOM_5,        dnds_list$Aom)
  annots$Nor <- gff2txdb_Norway("OdB3/Oikopleura_annot_v1.0.gff.gz", OdB3,         dnds_list$Nor)
  
  annots
}
annots2 = loadAnnotationsWithDNDS(a)
