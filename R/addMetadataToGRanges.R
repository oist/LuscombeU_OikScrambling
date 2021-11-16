#' Add metadata to a GenomicRanges object
#'
#' Function to add metadata columns to a GenomicRanges object. Potentially useful for `GenomicBreaks`.
#'
#' @returns An annotated GenomicRanges object with a new metadata column added.
#'

addMetadataToGRanges <- function(gr, df, gr_key="gene_id", df_key="gene_id", df_col_name=NA, rename_meta=NA) {
  # The input to this function is a GenomicRanges object and a dataframe.
  # The two must share some unique identifier (a key) in common; the dataframe
  # keys must be a subset of the GenomicRanges' keys.
  # Using these keys, a column from dataframe (the column with the name df_col_name)
  # is appended to the GenomicRanges' metadata columns.
  # Optionally, the appended column can be renamed to metadata_colname.
  gr_keys = mcols(gr)[[gr_key]]
  df_keys = as.list(df)[[df_key]]

  # Ensure all df keys are within gr_keys.
  if(!all(df_keys %in% gr_keys)){
    stop("Error. Looks like the keys in the dataframe don't match the keys in the GenomicRanges?")
  }

  # TODO add logic to check if column names are unique

  # Ensure the df column name is actually in the df.
  if(is.na(df_col_name)){
    stop("Error. A column name within the metadata data frame must be specified.")
  }
  if(! df_col_name %in% colnames(df)){
    stop(paste("Error. The column name ", df_col_name, " does not seem to be in the metadata data frame. Available columns:\n", paste(colnames(df), collapse=", ")))
  }

  # It's not easy to assign a column with a variable name to the granges metadata columns.
  # So, assign the metadata to a temporary column, then overwrite the column name later.
  gr$tmp_col_name = NA
  gr$tmp_col_name = df[match(gr_keys, df_keys),][[df_col_name]]
  df_values = as.list(df)[[df_col_name]]

  # df_col_name is the column name within the data frame to be added to the GenomicRanges.
  if(is.na(rename_meta)){
    names(mcols(gr)) = gsub("tmp_col_name", df_col_name, names(mcols(gr)))
  } else {
    names(mcols(gr)) = gsub("tmp_col_name", rename_meta, names(mcols(gr)))
  }
  gr
}
