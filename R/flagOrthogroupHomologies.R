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
#' @returns A data frame with phylogenetic hierarchical orthogroups and annotations.

flagOrthogroupHomologies <- function(file, clades=NULL, species_name_map=NULL, multiple_clade_resolution="largest", deletion_clade_resolution="largest") {
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

  if(!multiple_clade_resolution %in% c("largest", "mean")){
     stop("Error. The only valid options for this function are \"largest\" or \"mean\".\n\"largest\" breaks ties by assigning it to clade that includes the most species.\n\"mean\" breaks ties by assigning it to the clade with the highest mean count.")
  }

  if(!deletion_clade_resolution %in% c("largest", "mean")){
    stop("Error. The only valid options for this function are \"largest\" or \"mean\".\n\"largest\" breaks ties by assigning it to clade that includes the most species.\n\"mean\" breaks ties by assigning it to the clade with the lowest mean count.")
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

    # Count the counts (e.g. count1=5, count2=1, NA=0)
    # Make this a list (index = count, contents = count)
    count_freq = as.list(table(counts))

    # Check whether the counts are consistent with a named clade.
    # Loop over every clade, and check whether all of that clade's
    # counts exceed all of the counts outside the clade.
    # Yields a logical for every clade.
    check_clades = function(clade_list, counts, comparison="gt") {
      satisfied_clades = lapply(names(clade_list), function(clade){
        result = SimpleList()
        species = clade_list[[clade]]
        clade_counts = counts[species]
        nonclade_counts = counts[!names(counts) %in% species]
        if(comparison == "gt") {
          if(all(clade_counts > nonclade_counts)){
            result$clade_counts     = paste0(clade_counts, collapse=',')
            result$nonclade_counts  = paste0(nonclade_counts, collapse=',')
            result$clade_size       = length(clade_counts)
            result$mean_clade_count = mean(clade_counts)
            result$satisfied        = TRUE
          } else {
            result$clade_counts     = paste0(clade_counts, collapse=',')
            result$nonclade_counts  = paste0(nonclade_counts, collapse=',')
            result$clade_size       = length(clade_counts)
            result$mean_clade_count = mean(clade_counts)
            result$satisfied        = FALSE
          }
        } else if(comparison == "lt") {
          if(all(clade_counts < nonclade_counts)){
            result$clade_counts     = paste0(clade_counts, collapse=',')
            result$nonclade_counts  = paste0(nonclade_counts, collapse=',')
            result$clade_size       = length(clade_counts)
            result$mean_clade_count = mean(clade_counts)
            result$satisfied        = TRUE
          # The following block is kind of speculative. It is used later to call
          # clade-specific counts - but just because all counts in one clade are
          # 0, that does not make a deletion lineage-specific. The only way to
          # truly call it clade-specific is if all in clade are 0 and all outside
          # are >0. Delete this comment and the following block if this behaviour
          # is undesirable.
          } else if(all(clade_counts == 0)) {
            result$clade_counts     = paste0(clade_counts, collapse=',')
            result$nonclade_counts  = paste0(nonclade_counts, collapse=',')
            result$clade_size       = length(clade_counts)
            result$mean_clade_count = mean(clade_counts)
            result$satisfied        = TRUE
          } else {
            result$clade_counts     = paste0(clade_counts, collapse=',')
            result$nonclade_counts  = paste0(nonclade_counts, collapse=',')
            result$clade_size       = length(clade_counts)
            result$mean_clade_count = mean(clade_counts)
            result$satisfied        = FALSE
          }
        }
        result
      })
      names(satisfied_clades) <- names(clade_list)
      satisfied_clades = as.data.frame(do.call(rbind, satisfied_clades))
      return(satisfied_clades)
    }

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

    # If none missing...
    if(all(counts != 0)){
      present_in_all_species = "present_in_all_spp"
      # If all species have 1 of this HOG
      if(all(counts==1)){
        flag = "one_to_one"
      # If more than one species has >1 copy...
      } else {
        # Check if counts are consistent with a clade.
        clade_check = check_clades(clade_list=clades, counts=counts, comparison="gt")
        satisfied_clades = clade_check[unlist(clade_check$satisfied)==TRUE,]
        # If any clade is satisfied...
        if(any(unlist(clade_check$satisfied))){
          # If more than one clade filter is satisfied...
          if(nrow(satisfied_clades) > 1){
            # Then settle the tie, either by selecting the biggest clade ("largest")
            # or by the clade with the highest mean gene count ("mean").
            possible_clades = rownames(satisfied_clades)
            if(multiple_clade_resolution=="largest"){
              clade = rownames(head(satisfied_clades[order(unlist(satisfied_clades$clade_size), decreasing=T),], n=1))
              flag = paste("expansion_", clade,  sep="")
              warning_message = paste("HOG ", hog, " has counts that satisfy more than one clade:\n", paste(capture.output(print(satisfied_clades)), collapse='\n'), "\nBreaking tie in favour of largest clade: ", clade, "\nFlag: ", flag, "\n", sep="")
              warning(warning_message)
            } else if(multiple_clade_resolution=="mean") {
              clade = rownames(head(satisfied_clades[order(unlist(satisfied_clades$mean_clade_count), decreasing=T),], n=1))
              flag = paste("expansion_", clade,  sep="")
              warning_message = paste("HOG ", hog, " has counts that satisfy more than one clade:\n", paste(capture.output(print(satisfied_clades)), collapse='\n'), "\nBreaking tie in favour of clade with greatest mean count: ", clade, "\nFlag: ", flag, "\n", sep="")
              warning(warning_message)
            }
          # If only one clade is satisfied, then the flag is that clade.
          } else {
            flag = paste("expansion_", rownames(satisfied_clades),  sep="")
          }
        } else {
          # If there is a single species that has the maximum number of copies, it is a species-specific HOG.
          if(names(count_freq[which.max(count_freq)]) == 1){
            flag = paste("expansion_", names(which.max(counts)), sep="")
          } else {
            # If a HOG family is not specific to any clade, it is labeled as a difficult case.
            flag = "irreconcilable_paralogue"
          }
        }
      }

    # If any missing...
    } else {
      present_in_all_species = "absent_in_1_or_more_spp"
      # If a HOG is found in only one species, it is a species-specific HOG.
      if(count_freq[["0"]] == num_species-1){
        flag = paste("species_specific_", names(which.max(counts)), sep="")
      } else {
        # If a single species has 0 while every other species has >=1
        if(count_freq[["0"]] == 1 ){
          flag = paste("deletion_", names(which.min(counts)), sep="")
        } else {
          # Check for clade-specific HOG, allowing for missing genes in one or more species.
          clade_check = check_clades(clade_list=clades, counts=counts, comparison="gt")
          satisfied_clades = clade_check[clade_check$satisfied==TRUE,]
          # If gene counts consistent with any clade...
          if(any(unlist(satisfied_clades$satisfied))){
            # If more than one clade is satisfied....
            if(nrow(satisfied_clades)>1){
              # Then settle the tie, either by selecting the biggest clade ("largest")
              # or by the clade with the highest mean gene count ("mean").
              possible_clades = rownames(satisfied_clades)
              if(multiple_clade_resolution=="largest"){
                clade = rownames(head(satisfied_clades[order(unlist(satisfied_clades$clade_size), decreasing=T),], n=1))
                flag = paste("expansion_", clade,  sep="")
                warning_message = paste("HOG ", hog, " has counts that satisfy more than one clade:\n", paste(capture.output(print(satisfied_clades)), collapse='\n'), "\nBreaking tie in favour of largest clade: ", clade, "\nFlag: ", flag, "\n", sep="")
                warning(warning_message)
              } else if(multiple_clade_resolution=="mean") {
                clade = rownames(head(satisfied_clades[order(unlist(satisfied_clades$mean_clade_count), decreasing=T),], n=1))
                flag = paste("expansion_", clade,  sep="")
                warning_message = paste("HOG ", hog, " has counts that satisfy more than one clade:\n", paste(capture.output(print(satisfied_clades)), collapse='\n'), "\nBreaking tie in favour of clade with greatest mean count: ", clade, "\nFlag: ", flag, "\n", sep="")
                warning(warning_message)
              }
            # If only one clade is satisfied, then it is a clade-specific HOG.
            } else {
              flag = paste("clade_specific_", rownames(satisfied_clades), sep="")
            }
          # If no clades are satisfied...
          } else {
            # Check for clade-specific deletion.
            # Note that here we check if ALL species in the clade are lower than ALL other clades.
            deletion_check = check_clades(clade_list=clades, counts=counts, comparison="lt")
            satisfied_deletions = deletion_check[deletion_check$satisfied==TRUE,]
            #print(deletion_check)
            # If more than one clade has a deletion...
            if(any(unlist(satisfied_deletions$satisfied))){
              # Settle the tie.
              # The tie should in almost all situations be resolved such that deletion is inferred
              # to belong to the largest clade. But there is the option to use something else if for
              # some reason this is not desired.
              if(nrow(satisfied_deletions) > 1){
                possible_clades = rownames(satisfied_deletions)
                if(deletion_clade_resolution=="largest"){
                  clade = rownames(head(satisfied_deletions[order(unlist(satisfied_deletions$clade_size), decreasing=T),], n=1))
                  flag = paste("clade_specific_deletion_", clade,  sep="")
                  warning_message = paste("HOG ", hog, " seems to be deleted in more than one clade:\n", paste(capture.output(print(satisfied_deletions)), collapse='\n'), "\nBreaking tie in favour of largest clade: ", clade, "\nFlag: ", flag, "\n", sep="")
                  warning(warning_message)
                } else if(deletion_clade_resolution=="mean") {
                  clade = rownames(head(satisfied_deletions[order(unlist(satisfied_deletions$mean_clade_count), decreasing=T),], n=1))
                  flag = paste("clade_specific_deletion_", clade,  sep="")
                  warning_message = paste("HOG ", hog, " seems to be deleted in more than one clade:\n", paste(capture.output(print(satisfied_deletions)), collapse='\n'), "\nBreaking tie in favour of clade with greatest mean count: ", clade, "\nFlag: ", flag, "\n", sep="")
                  warning(warning_message)
                }
              # If only one clade seems to have it deleted, it is a clade-specific deletion.
              } else {
                flag = paste("clade_specific_deletion_", rownames(satisfied_deletions),  sep="")
              }
            } else {
              # If there is a non-clade-specific deletion - for any reason - it is labeled as a difficult case.
              flag = "irreconcilable_patchiness"
            }
          }
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
    result$taxa <- taxa_string
    result$taxa_count <- taxa_count_string
    result$flag <- flag
    result$present_in_all_species <- present_in_all_species
    return(result)
  }
  res = setNames(apply(pho, 1, count_flag_orthofinder), pho$HOG) |> SimpleList()
  res = do.call(rbind, res)
  res
}
