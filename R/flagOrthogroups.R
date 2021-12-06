#' Flag genes by duplication status in Orthofinder table
#'
#' A suite of functions to flag orthogroups by their inferred gene histories
#' which may be useful later in the `BreakpointsData` package.
#'
#' @author Michael Mansfield
#'
#' The gist of these functions is that they annotate an OrthoFinder
#' hierarchical orthogroup (HOG) table with flags.
#'
#' More specifically, the function attempts to categorize the patterns
#' of gene family evolution using these counts. In pseudo-code, the logic
#' follows like so:
#'
#' If all counts = 1:                   one-to-one orthologue
#' If 1 species > all(other species):   species-specific expansion
#' If >1 species with > 1:
#'  If species constitute a clade:      clade-specific expansion
#'  Else:                               irreconcilable paralogue
#'
#' If 1 species = 0 & all others > 0:   species-specific deletion
#' If >1 species = 0:
#'  If species constitute a clade:      clade-specific deletion
#'  Else:                               irreconcilable patchiness
#'
#' @returns A matrix with annotations.

# Iterates over rows in a HOG count matrix and returns whether all species
# have at least one gene in the HOG (TRUE) or not (FALSE).
flagPHOPresenceAbsence <- function(count_matrix) {
  flags <- apply(count_matrix, 1, function(row) {
    if(all(row > 0)){
      return(TRUE)
    } else {
      return(FALSE)
    }
  })
  flags
}

# Iterate over rows in a HOG count matrix and return whether all species
# have exactly 1 gene (TRUE) or not (FALSE).
flagPHOOneToOnes <- function(count_matrix) {
  flags <- apply(count_matrix, 1, function(row) {
    if(all(row == 1)){
      return(TRUE)
    } else {
      return(FALSE)
    }
  })
  flags
}

# Iterate over rows in a HOG count matrix and return whether all species
# at least one gene AND at least one species has more than 1 gene (TRUE),
# or not (FALSE).
flagPHODuplications <- function(count_matrix) {
  flags <- apply(count_matrix, 1, function(row) {
    if(all(row > 0) & any(row > 1)){
      return(TRUE)
    } else {
      return(FALSE)
    }
  })
  flags
}

# Iterate over rows in a HOG count matrix and return whether any species
# has more genes than every other species (expansion_speciesName)
# or not (NA).
flagPHOSpeciesSpecificExpansions <- function(count_matrix) {
  flags <- apply(count_matrix, 1, function(row) {
    max_count <- row[which.max(row)]
    other_counts <- row[-which.max(row)]
    max_taxon <- names(max_count)
    if(all(max_count > other_counts)) {
      flag <- max_taxon
    } else {
      flag <- NA
    }
    flag
  })
  flags
}

# This function uses the rows in a HOG count matrix and a list of clades.
# The list of clades is a named list where each element contains a
# vector of species names.
# For each of the clades, this function checks whether the HOG count
# for ALL species in each clade exceeds all other clades.
flagPHOCladeSpecificExpansions <- function(count_matrix, clade_name, clade_spp) {
  matching_cols <- match(clade_spp, colnames(count_matrix))
  flags <- apply(count_matrix, 1, function(row) {
    clade_counts <- row[matching_cols]
    other_counts <- row[-matching_cols]
    if(min(clade_counts) > max(other_counts)) {
      flag <- TRUE
    } else {
      flag <- FALSE
    }
    flag
  })
  flags
}

# Convenience wrapper for flagging a bunch of clade expansions at once
flagPHOAllCladeExpansions <- function(count_matrix, clade_list) {
  flag_matrix <- do.call(cbind, lapply(names(clade_list), function(clade_name) {
    clade_flags <- flagPHOCladeSpecificExpansions(count_matrix, clade_name = clade_name, clade_spp = clade_list[[clade_name]])
    clade_flags
  }))
  colnames(flag_matrix) <- paste("expansion_clade_", names(clade_list), sep='')
  flag_matrix
}

# This function checks clade-specific expansion flags for each clade in
# clade list, and uses their HOG counts from the count matrix to reconcile
# when two or more clades are TRUE.
# I would recommend using as input the output of `flagPHOAllCladeExpansions`.
# It reconciles these flags using a reconciliation type, either "largest_clade"
# or "mean_count".
#     If "largest_clade": Final expansion flag is the clade with
#                         the most constituent species.
#     If "mean_count":    Final expansion flag is the clade with
#                         the highest average count.
reconcilePHOExpansions <-  function(count_matrix, clade_list, reconciliation="largest_clade") {
  flag_matrix <- flagPHOAllCladeExpansions(count_matrix, clade_list=clade_list)
  expansion_reconciled  <- lapply(1:nrow(count_matrix), function(row) {
    counts <- count_matrix[row,]
    flags  <- flag_matrix[row,]
    # If more than one expansion flag is true...
    if(length(which(flags==TRUE))>1){
      matching_clades = clade_list[which(flags==TRUE)]
      matching_clades_lengths = unlist(lapply(matching_clades, length))
      matching_clade_counts = lapply(matching_clades, function(clade) count_matrix[row,clade])
      matching_clade_means = lapply(matching_clade_counts, function(x) mean(unlist(x)))
      if(reconciliation=="largest_clade") {
        # Need to make sure the "largest" clade is actually larger than every other clade.
        if(all(max(matching_clades_lengths) > matching_clades_lengths[-which.max(matching_clades_lengths)])) {
          reconciled_flag = paste("expansion_clade_", names(matching_clades[which.max(matching_clades_lengths)]), sep="")
          warning_message = paste("Counts satisfy more than one clade:\n", paste(capture.output(print(counts)), collapse="\n"), "\n", paste(capture.output(print(flags)), collapse='\n'), collapse="\n", "\nBreaking tie in favour of largest clade: ", names(matching_clades[which.max(matching_clades_lengths)]), "\nFlag: ", reconciled_flag, "\n", sep="")
          warning(warning_message)
        } else {
          reconciled_flag = paste0(names(matching_clades), collapse="_")
          warning_message = paste("Counts satisfy more than one clade:\n", paste(capture.output(print(counts)), collapse="\n"), "\n", paste(capture.output(print(flags)), collapse='\n'), collapse="\n", "\nMore than one equally-sized clade is supported, can't break tie. Returning concatenate: ", paste0(names(matching_clades), collapse="_"), "\nFlag: ", reconciled_flag, "\n", sep="")
          warning(warning_message)
        }
      } else if(reconciliation == "mean_count") {
        reconciled_flag = paste("expansion_clade_", names(matching_clades[which.max(matching_clade_means)]), sep="")
        warning_message = paste("Counts satisfy more than one clade:\n", paste(capture.output(print(counts)), collapse="\n"), "\n", paste(capture.output(print(flags)), collapse='\n'), collapse="\n", "\nBreaking tie in favour clade with highest mean: ", names(matching_clades[which.max(matching_clade_means)]), "\nFlag: ", reconciled_flag, "\n", sep="")
        warning(warning_message)
      }
    } else if(length(which(flags==TRUE)) == 1) {
      reconciled_flag = names(which(flags==TRUE))
    } else {
      reconciled_flag = NA
    }
    reconciled_flag
  })
  flag_matrix = cbind(flag_matrix, expansion_reconciled)
  flag_matrix
}

# Iterate over rows in a HOG count matrix and return whether any species
# has 0 genes while all others have >0 (TRUE) or not (FALSE).
flagSpeciesSpecificDeletions <- function(count_matrix) {
  flags <- apply(count_matrix, 1, function(row) {
    min_count <- row[which.min(row)]
    other_counts <- row[-which.min(row)]
    min_taxon <- names(min_count)
    if(min_count == 0 & all(other_counts > 0)) {
      flag <- min_taxon
    } else {
      flag <- NA
    }
    flag
  })
  flags
}

# This function uses the rows in a HOG count matrix and a list of clades.
# The list of clades is a named list where each element contains a
# vector of species names.
# For each of the clades, this function checks whether the HOG count
# for ALL species in each clade is exactly 0.
flagPHOCladeSpecificDeletions <- function(count_matrix, clade_name, clade_spp) {
  matching_cols <- match(clade_spp, colnames(count_matrix))
  flags <- apply(count_matrix, 1, function(row) {
    clade_counts <- row[matching_cols]
    other_counts <- row[-matching_cols]
    if(all(clade_counts==0)) {
      flag <- TRUE
    } else {
      flag <- FALSE
    }
    flag
  })
  flags
}

# Convenience wrapper for flagging a bunch of clade deletions at once
flagPHOAllCladeDeletions <- function(count_matrix, clade_list) {
  flag_matrix <- do.call(cbind, lapply(names(clade_list), function(clade_name) {
    clade_flags <- flagPHOCladeSpecificDeletions(count_matrix, clade_name = clade_name, clade_spp = clade_list[[clade_name]])
    clade_flags
  }))
  colnames(flag_matrix) <- paste("deletion_clade_", names(clade_list), sep='')
  flag_matrix
}

# This function is similar to reconcilePHOExpansions, but needs some extra logic.
reconcilePHODeletions <-  function(count_matrix, clade_list, reconciliation="largest_clade", verbose=F) {
  flag_matrix <- flagPHOAllCladeDeletions(count_matrix, clade_list=clade_list)
  deletion_reconciled  <- lapply(1:nrow(count_matrix), function(row) {
    counts <- count_matrix[row,]
    flags  <- flag_matrix[row,]
    # If more than one clade has a deletion flagged, we have to get fancy.
    if(length(which(flags==TRUE))>1){
      matching_clades = clade_list[which(flags==TRUE)]
      matching_clades_lengths = unlist(lapply(matching_clades, length))
      if(reconciliation=="largest_clade") {
        # Need to make sure the "largest" clade is actually larger than every other clade.
        if(all(max(matching_clades_lengths) > matching_clades_lengths[-which.max(matching_clades_lengths)])) {
          reconciled_flag = paste("deletion_clade_", names(matching_clades[which.max(matching_clades_lengths)]), sep="")
          warning_message = paste("Counts satisfy more than one clade:\n", paste(capture.output(print(counts)), collapse="\n"), "\n", paste(capture.output(print(flags)), collapse='\n'), collapse="\n", "\nBreaking tie in favour of largest clade: ", names(matching_clades[which.max(matching_clades_lengths)]), "\nFlag: ", reconciled_flag, "\n", sep="")
          if(verbose){
            warning(warning_message)
          }
        } else {
          reconciled_flag = paste("deletion_clade_", paste0(names(matching_clades), collapse="_"), sep="")
          warning_message = paste("Counts satisfy more than one clade:\n", paste(capture.output(print(counts)), collapse="\n"), "\n", paste(capture.output(print(flags)), collapse='\n'), collapse="\n", "\nMore than one equally-sized clade is supported, can't break tie. Returning concatenate: ", paste0(names(matching_clades), collapse="_"), "\nFlag: ", reconciled_flag, "\n", sep="")
          if(verbose){
            warning(warning_message)
          }
        }
      }
    # If only one flag is TRUE, then the reconciled flag is that one
    } else if(length(which(flags==TRUE)) == 1) {
      reconciled_flag = names(which(flags==TRUE))
    } else {
      reconciled_flag = NA
    }
    reconciled_flag
  })
  flag_matrix = cbind(flag_matrix, deletion_reconciled)
  flag_matrix
}

flagAllPHO <- function(count_matrix, clade_list, expansion_reconciliation="largest_clade", deletion_reconciliation="largest_clade", deletion_warnings=F){
  annotated_matrix                                <- as.data.frame(matrix(nrow=nrow(count_matrix)))
  colnames(annotated_matrix) <- "present_in_all_spp"
  annotated_matrix[,"present_in_all_spp"]         <- flagPHOPresenceAbsence(count_matrix)
  annotated_matrix[,"one_to_one"]                 <- flagPHOOneToOnes(count_matrix)
  annotated_matrix[,"common_gene_duplication"]    <- flagPHODuplications(count_matrix)
  annotated_matrix[,"expansion_species_specific"] <- flagPHOSpeciesSpecificExpansions(count_matrix)
  annotated_matrix                                <- cbind(annotated_matrix, reconcilePHOExpansions(count_matrix=count_matrix, clade_list = clade_list, reconciliation = expansion_reconciliation))
  annotated_matrix[,"deletion_species_specific"]  <- flagSpeciesSpecificDeletions(count_matrix)
  annotated_matrix                                <- cbind(annotated_matrix, reconcilePHODeletions(count_matrix=count_matrix, clade_list = clade_list, reconciliation=deletion_reconciliation, verbose=deletion_warnings))
  annotated_matrix
}
