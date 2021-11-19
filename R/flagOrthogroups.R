#' Flag genes by duplication status in Orthofinder table
#'
#' Creates useful flags for orthologues in the `BreakpointsData` package.
#'
#' @author Michael Mansfield
#'
#' The gist of this function is that it reads an OrthoFinder hierarchical
#' orthogroup (HOG) table and flags genes based on their gene counts.
#'
#' More specifically, the function attempts to categorize the types of
#' gene families that emerge from the count table. For example, all species
#' having one copy means it is a 1-to-1 orthologue. The more complicated logic
#' is explained in detail in the function.
#'
#' @returns IDK yet

flagOrthofinderGenes <- function(file, clades=NULL, species_name_map=NULL) {
  # Read table. Note that empty cells become NA.
  pho <- read.delim(file, check.names=FALSE, na.strings=c("", NA))
  # Remove all-NA columns
  pho <- pho[,sapply(pho, \(col) ! all(is.na(col)))]

  # The first 3 columns of an OrthoFinder HOG table are always HOG/OG/Parent.
  # So the number of species is (number of columns - 3).
  num_species = ncol(pho)-3

  # Rename species names, if desired.
  # Matches old names to a set of new names.
  # The new names are later checked in the clades, not the old ones.
  if(!is.null(species_name_map)){
    old_names = colnames(pho)[4:ncol(pho)]
    if(!all.equal(sort(old_names), sort(species_name_map$old_names))){
      print(paste("HOG columns:      ", paste0(old_names, collapse=',')))
      print(paste("species_name_map: ", paste0(species_name_map$old_names, collapse=',')))
      stop("Error. The column names in the HOG table do not match the species_name_map?")
    } else {
      new_names = species_name_map[match(old_names, species_name_map$old_names),]$new_names
      colnames(pho)[4:ncol(pho)] = new_names
    }
  }

  # Make sure that all identifiers listed in the clades are found in the
  # HOG table. It uses the HOG table's new column names if they were renamed.
  if(!all(unique(unlist(clades)) %in% colnames(pho))){
    print(paste("HOG columns:       ", paste0(old_names, collapse=',')))
    print(paste("Clade identifiers: ", paste0(unique(unlist(clades)), collapse=',')))
    stop("Error. Not all names in the specified clades match column names in the HOG table?")
  }

  # The next lines just count how many entries (genes) there are in every cell across the
  # HOG table.
  count_flag_orthofinder <- function(row){
    hog <- row[["HOG"]]
    species <- names(row[4:length(row)])
    cells <- row[4:length(row)]
    counts <- sapply(cells, function(cell){
      if(is.na(cell)){
        ngenes <- 0
      } else {
        ngenes <- length(unlist(strsplit(cell, split=', ')))
      }
      ngenes
    })

    # Here we lay out the logic for how HOG counts relate to gene family evolution.
    # Some of the "flags" are intuitive - e.g. 1-to-1 - some less so.
    # The function uses a list of clades and evaluates whether there has been
    # clade-specific expansions or contractions of that HOG.
    #
    #   If none missing (all 6 species >= 1):
    #     If all 6 species = 1:                                     1-to-1 orthologue
    #     If only one species has >1:                               species-specific duplication(s)
    #     If more than one species has >1:
    #          If ALL in one clade more than ANY outside the clade: clade-specific expansion
    #          E.g. Oki (4) Kum (5) Nor (1) Bar (1) Osa (2) Aom (2)
    #            = expansion in Oki-Kum
    #          E.g. Oki (1) Kum (1) Nor (4) Bar (4) Osa (3) Aom (3)
    #            = expansion in Nor-Bar-Osa-Aom
    #          The inverse could be true (clade-specific contraction!)
    #     Else:                                                     some irreconcilable XtoX paralogous relationship
    #   If any missing:
    #        If missing in 1 species:                               probable misannotation
    #        If ALL in one group missing and ALL outside group >=1: clade-specific deletion
    #        If ALL missing except one species:                     species-specific HOG
    #        Other:                                                 some other explanation, misannotation, etc.

    # Count the counts (e.g. count1=5, count2=1, NA=0)
    # Make this a list (index = count, contents = count)
    count_freq = as.list(table(counts, useNA="always"))
    # If none missing...
    if(all(counts != 0)){
      # Pretty self-explanatory
      if(all(counts==1)){
        flag = "one_to_one"
      # Check whether there are 1-count entries in the HOG (not all HOGs have 1-counts).
      } else if("1" %in% names(count_freq)){
        # If only a single species has >1 gene, it is a species-specific HOG.
        if(count_freq[["1"]] == num_species-1){
          flag = paste("expansion_", names(which.max(counts)), sep="")
        }
      # If more than one species has >1 copy...
      } else {
        # Check whether the counts are consistent with a named clade.
        # Loop over every clade, and check whether all of that clade's
        # counts exceed all of the counts outside the clade.
        # Yields a logical vector for every clade.
        satisfied_clades = unlist(lapply(names(clades), function(clade){
          species = clades[[clade]]
          clade_counts = counts[species]
          nonclade_counts = counts[!names(counts) %in% species]
          if(all(clade_counts > nonclade_counts)){
            return(TRUE)
          } else {
            return(FALSE)
          }
        }))
        names(satisfied_clades) <- names(clades)

        if(any(satisfied_clades)){
          if(sum(satisfied_clades, na.rm=T) > 1){
            possible_clades = names(which(c(satisfied_clades)))
            biggest_clade   = which.max(unlist(lapply(clades[match(possible_clades, clades)], length)))
            flag = paste("expansion_", names(clades)[biggest_clade],  sep="")
          } else {
            flag = paste("expansion_", names(which(c(satisfied_clades))),  sep="")
          }
        } else {
          flag = "irreconcilable_paralogue"
        }
      }
    # TODO deal with missing genes
    # If any missing...
    } else {
      if(all(counts == 1)){
        flag = "blah"
      } else {
        flag = "nah"
      }
    }
    #print(c(counts, flag=flag))

    #transcript_names = SimpleList(as.list(cells))
    #hog_list = SimpleList()
    #hog_list$counts = counts
    #hog_list
  }
  res = apply(head(pho, n=10000), 1, count_flag_orthofinder)
  res
  # Flag interesting entries
  #pho$orthofinder_counts_per_species =
  #pho$orthofinder_transcripts_in_hog =
  #pho$orthofinder_hog_flag           =
}
clades=list("Oki_Kum"=c("Kum", "Oki"), "Osa_Aom"=c("Osa","Aom"), "Bar_Nor"=c("Bar","Nor"), "Northern"=c("Bar","Nor","Osa","Aom"))

species_name_map = data.frame(old_names=c("AOM-5-5f.prot.longest.fa_1", "Bar2_p4.Flye.prot.longest.fa_1", "KUM-M3-7f.prot.longest.fa_1", "OKI2018_I69.v2.prot.longest.fa_1", "OSKA2016v1.9.prot.longest.fa_1", "OdB3.v1.0.prot.fa_1.nohaplo"),
                              new_names=c("Aom", "Bar", "Kum", "Oki", "Osa", "Nor")
                              )

b = flagOrthofinderGenes(system.file("extdata/OrthoFinder/N19.tsv", package = "BreakpointsData"), clades=clades, species_name_map=species_name_map)

