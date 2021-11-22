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

  # Make sure that none of the clades consist of exactly the same species identifiers.
  if( length(clades) != length(unique(lapply(clades, sort))) ){
    problem_clade = lapply(names(clades), function(x) {
      taxa = paste(sort(clades[[x]]), collapse=',')
      informative_df = data.frame("clade"=x, "sorted_taxa"=taxa)
      #b = paste(x, ": ", taxa, sep="")
      return(informative_df)
    })
    problem_clades = do.call(rbind, problem_clade)
    print(problem_clades)
    stop("Error. Not all clades contain unique taxa?")
  }

  # The next lines just count how many entries (genes) there are in every cell across the
  # HOG table.
  count_flag_orthofinder <- function(row){
    # Gather some information from the non-count cells for later.
    hog <- row[["HOG"]]
    og <- row[["OG"]]
    gtpc <- row[["Gene Tree Parent Clade"]]
    species <- names(row[4:length(row)])
    cells <- row[4:length(row)]
    transcript_names <- paste0(na.omit(cells), collapse=',')

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
    count_freq = as.list(table(counts))
    # If none missing...
    if(all(counts != 0)){
      present_in_all_species = "present_in_all_spp"
      # If all species have 1 of this HOG
      if(all(counts==1)){
        flag = "one_to_one"
      # If more than one species has >1 copy...
      } else {
        # Check whether the counts are consistent with a named clade.
        # Loop over every clade, and check whether all of that clade's
        # counts exceed all of the counts outside the clade.
        # Yields a logical for every clade.
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
          # This block checks whether one or more clades are satisfied.
          # If more than one clade filter is satisfied, then it defaults to the largest-possible
          # clade. E.g. if Bar-Nor is satisfied and Bar-Nor-Osa-Aom is satisfied, the latter
          # gets the flag.
          if(sum(satisfied_clades, na.rm=T) > 1){
            possible_clades = names(which(c(satisfied_clades)))
            matching_clades = clades[match(possible_clades, names(clades))]
            biggest_clade   = names(which.max(unlist(lapply(matching_clades, length))))
            flag = paste("expansion_", biggest_clade,  sep="")
            warning_message = paste("HOG ", hog, " has counts that satisfy more than one clade: ",  paste0(names(counts), collapse=", "), ", counts: ", paste0(counts, collapse=", "), "\nSatisfied clades: ", paste0(possible_clades, collapse=", "), "\nDefaulting to the biggest clade: ", biggest_clade, "\nFlag: ", flag, "\n", sep="")
            warning(warning_message)
            # If only one clade is satisfied, then the flag is that clade.
          } else {
            flag = paste("expansion_", names(which(c(satisfied_clades))),  sep="")
          }
        } else {
          # If only a single species has >1 gene, it is a species-specific HOG.
          if(names(count_freq[which.max(count_freq)]) == 1){
            flag = paste("expansion_", names(which.max(counts)), sep="")
            #print(flag)
          } else {
            # If a HOG family is not specific to any clade, it is labeled as a difficult case.
            flag = "irreconcilable_paralogue"
          }
        }
      }
    # If any missing...
    } else {
      present_in_all_species = "absent_in_1_or_more_spp"
      # If 1 or more genes are found in only one taxon, it is a lineage-specific gene.
      if(count_freq[["0"]] == num_species-1){
        flag = paste("lineage_specific_", names(which.max(counts)), sep="")
      } else {
        # Check whether counts are consistent with a named clade.
        # Loop over every clade, and check whether all of that clade's
        # counts exceed all of the counts outside the clade.
        # Yields a logical for every clade.
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
          # This block checks whether one or more clades are satisfied.
          # If more than one clade filter is satisfied, then it defaults to the largest-possible
          # clade. E.g. if Bar-Nor is satisfied and Bar-Nor-Osa-Aom is satisfied, the latter
          # gets the flag.
          if(sum(satisfied_clades, na.rm=T) > 1){
            possible_clades = names(which(c(satisfied_clades)))
            matching_clades = clades[match(possible_clades, names(clades))]
            biggest_clade = names(which.max(unlist(lapply(matching_clades, length))))
            flag = paste("lineage_specific_", biggest_clade,  sep="")
            warning_message = paste("HOG ", hog, " has counts that satisfy more than one clade: ",  paste0(names(counts), collapse=", "), ", counts: ", paste0(counts, collapse=", "), "\nSatisfied clades: ", paste0(possible_clades, collapse=", "), "\nDefaulting to the biggest clade: ", biggest_clade, "\nFlag: ", flag, "\n", sep="")
            warning(warning_message)
          # If only one clade is satisfied, then the flag is that clade.
          } else {
            flag = paste("lineage_specific_", names(which(c(satisfied_clades))),  sep="")
          }
        } else {
          # If there is a non-clade-specific contraction - for any reason - it is labeled as a difficult case.
          flag = "irreconcilable_patchiness"
        }
      }
    }

    # Make some nicely-formatted strings for appending to the HOG table.
    # Taxa in which the HOG exists
    taxa_string <- paste0(names(counts[which(counts>0)]), collapse=',')
    # Taxon counts within HOG
    taxa_count_string <- paste0(lapply(names(counts), function(x) {
      count <- counts[[x]]
      paste0(x, ":", count)
    }), collapse=",")
    result <- data.frame(hog=hog, og=og, gene_tree_parent_clade=gtpc)
    # All transcript IDs for all taxa within HOG
    taxa_transcripts <- t(as.data.frame(cells))
    rownames(taxa_transcripts) <- NULL
    result <- cbind(result, taxa_transcripts)
    result$transcripts <- transcript_names
    #
    result$taxa <- taxa_string
    result$taxa_count <- taxa_count_string
    result$flag <- flag
    result$present_in_all_species <- present_in_all_species
    return(result)
  }
  res = setNames(apply(pho, 1, count_flag_orthofinder), pho$HOG)
  res |> SimpleList()
}
clades=list("Oki_Kum"=c("Kum", "Oki"), "Osa_Aom"=c("Osa","Aom"), "Bar_Nor"=c("Bar","Nor"), "Northern"=c("Bar","Nor","Osa","Aom"))

species_name_map = data.frame(old_names=c("AOM-5-5f.prot.longest.fa_1", "Bar2_p4.Flye.prot.longest.fa_1", "KUM-M3-7f.prot.longest.fa_1", "OKI2018_I69.v2.prot.longest.fa_1", "OSKA2016v1.9.prot.longest.fa_1", "OdB3.v1.0.prot.fa_1.nohaplo"),
                              new_names=c("Aom", "Bar", "Kum", "Oki", "Osa", "Nor")
                              )

b = do.call(rbind, flagOrthofinderGenes(system.file("extdata/OrthoFinder/N19.tsv", package = "BreakpointsData"), clades=clades, species_name_map=species_name_map))

# ggplot(b) + aes(x=reorder(flag, flag, function(x) length(x))) + geom_bar(stat="count", position="dodge") + facet_wrap(~ present_in_all_species, scales = 'free' ) + xlab(NULL) + coord_flip()
