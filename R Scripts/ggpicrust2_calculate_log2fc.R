# The ggpicrust2 pathway_daa function is currently broken - it uses incorrect code to calculate log2 fold changes.
# This has been confirmed for Maaslin2, but I recommend using this function instead for all methods for two reasons:
# Firstly, you know exactly how it's calculated (log2(mean relative abundance of group 1/mean of group 2)).
# Secondly, I haven't strictly confirmed that ALL of ggpicrust2's other differential abundance tools have accurate log2 fold changes.

# INSTRUCTIONS ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# The inputs to this function are the same as the pathway_daa command. 
# It will return a table with the log2 fold changes for your variable. 
# The 'reference' variable is required - i.e. which group is the reference group? log2(interesting_group/reference_group)
# If you have more than two groups, it'll calculate the log2 fold changes for each group compared to the reference group. 

# Run pathway_daa, then run this, and then OVERWRITE the pathway_daa log2 fold changes with these values.
# (I'll leave the details of the overwrite to you)
fix_lfc = function(abundance, metadata, group, reference){
  
  # The inputs are the same as the pathway_daa command.
  # Uncomment the following to easily troubleshoot fix_lfc().
  # abundance = path_filt
  # metadata = metadata
  # group = "Tert_VC"
  # reference = "VC Low"
  
  # Ensure that the samples in metadata and abundance match, 
  # and identify the column in metadata that corresponds to the sample names
  aligned = align_samples(abundance, metadata, verbose = FALSE)
  abundance = aligned$abundance
  metadata = aligned$metadata
  snames = aligned$sample_col
  
  # Convert the abundance to relative abundance
  abundance_relab = abundance %>% apply(2, function(x) x/sum(x)) %>% as.data.frame()
  
  # How many levels are there in the outcome of interest?
  if(!is.numeric(metadata[[group]])){
    levels = unique(metadata[[group]])
    num_levels = length(unique(metadata[[group]]))
    
    # Calculate average abundance for each group
    for(i in 1:length(levels)){
      # i=2
      lvl = levels[i]
      group_samples = metadata %>% filter(!!sym(group) == lvl) %>% pull(!!sym(snames))
      if(i==1){
        means_per_group = rowMeans(abundance_relab[, group_samples, drop = FALSE],na.rm = TRUE) %>% 
          as.data.frame %>% `colnames<-`(lvl)
      } else {
        # Join by rownames
        means_per_group = merge(means_per_group, 
                                rowMeans(abundance_relab[, group_samples, drop = FALSE], 
                                         na.rm = TRUE) %>% 
                                  as.data.frame %>% `colnames<-`(lvl),
                                by = "row.names", all = TRUE) %>% 
          rename(feature = `Row.names`)
      }
      
    }
    
    # Calculate the log 2 fold changes with the reference group as the denominator
    pseudocount = min(abundance_relab[abundance_relab > 0], na.rm = TRUE) / 2
    log2fc = means_per_group %>%
      pivot_longer(-all_of(c('feature',reference)), names_to = "group1", values_to = "abundance") %>%
      mutate(group2 = reference) %>%
      mutate(log2_fold_change = log2((abundance + pseudocount) / (.[[reference]] + pseudocount))) %>%
      select(feature, group1, group2, log2_fold_change)
    
  } else {
    print('!!!!!!ATTENTION!!!!!!')
    print('Group variable is numeric - calculating log2 fold change as slope of linear regression.')
    print('This means that the log2 fold change represents the change in abundance per unit increase in the group variable, rather than relative to a reference group.')
    print('For example, if your grouping variable is age, the log2 fold change will represent the change in abundance per year of age.')
    print('The interpretation of the value is therefore NOT the same is a regular log2 fold change!')
    print("I haven't changed the name to something more appropriate so that these results can easily be integrated into other things downstream without changing your code. Just be aware that the interpretation is different!")
    
    log2fc = abundance_relab %>%
      rownames_to_column("feature") %>%
      pivot_longer(-feature, names_to = "sample", values_to = "abundance") %>%
      left_join(metadata %>% select(!!sym(snames), !!sym(group)), by = c("sample" = snames)) %>%
      group_by(feature) %>%
      summarize(log2_fold_change = coef(lm(abundance ~ get(group)))[2]) %>%
      ungroup() %>% 
      mutate(variable = group, measurement = "change in abundance per unit variable") %>%
      select(feature, variable, measurement, log2_fold_change)
  }
  
  # result = log2fc %>% rownames_to_column("feature")
  return(log2fc)
}

# Helper functions ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# This is ggpicrust2's align_samples function, which for some reason is inaccessible for some students.
# Typically you can access it with ggpicrust2:::align_samples. 
# Two colons (::) lets you access common functions from a package; three colons (:::) allows you to access hidden functions.
align_samples = function (abundance, metadata, sample_col = NULL, verbose = TRUE) {
  abundance_samples <- colnames(abundance)
  if (length(abundance_samples) == 0) {
    stop("Abundance data has no column names (sample identifiers)")
  }
  if (is.null(sample_col)) {
    sample_col <- find_sample_column(metadata, abundance_samples)
    if (is.null(sample_col)) {
      stop("Cannot find matching sample identifiers between abundance and metadata.\n", 
           "  Abundance samples (first 5): ", paste(head(abundance_samples, 
                                                         5), collapse = ", "), if (length(abundance_samples) > 
                                                                                   5) 
                                                           "..."
           else "", "\n", "  Metadata columns: ", paste(colnames(metadata), 
                                                        collapse = ", "), "\n", "  Metadata rownames (first 5): ", 
           paste(head(rownames(metadata), 5), collapse = ", "), 
           if (nrow(metadata) > 5) 
             "..."
           else "")
    }
    if (verbose) {
      if (sample_col == ".rownames") {
        message("Using metadata rownames as sample identifiers")
      }
      else {
        message(sprintf("Using column '%s' as sample identifier", 
                        sample_col))
      }
    }
  }
  if (sample_col == ".rownames") {
    metadata <- as.data.frame(metadata)
    metadata$.sample_id <- rownames(metadata)
    sample_col <- ".sample_id"
  }
  if (!sample_col %in% colnames(metadata)) {
    stop(sprintf("Sample column '%s' not found in metadata", 
                 sample_col))
  }
  metadata_samples <- as.character(metadata[[sample_col]])
  common_samples <- intersect(abundance_samples, metadata_samples)
  if (length(common_samples) == 0) {
    stop("No common samples found between abundance and metadata.\n", 
         "  Abundance samples: ", paste(head(abundance_samples, 
                                             3), collapse = ", "), "...\n", "  Metadata samples: ", 
         paste(head(metadata_samples, 3), collapse = ", "), 
         "...")
  }
  missing_in_abundance <- setdiff(metadata_samples, abundance_samples)
  missing_in_metadata <- setdiff(abundance_samples, metadata_samples)
  if (verbose && (length(missing_in_abundance) > 0 || length(missing_in_metadata) > 
                  0)) {
    message(sprintf("Aligned %d samples (dropped %d from metadata, %d from abundance)", 
                    length(common_samples), length(missing_in_abundance), 
                    length(missing_in_metadata)))
  }
  abundance_aligned <- abundance[, common_samples, drop = FALSE]
  row_idx <- match(common_samples, metadata_samples)
  metadata_aligned <- as.data.frame(metadata[row_idx, , drop = FALSE])
  rownames(metadata_aligned) <- common_samples
  if (!identical(colnames(abundance_aligned), as.character(metadata_aligned[[sample_col]]))) {
    stop("Internal error: sample alignment failed")
  }
  
  return(list(abundance = abundance_aligned, 
              metadata = metadata_aligned, 
              sample_col = sample_col, 
              n_samples = length(common_samples)))
}

# Also a helper function from ggpicrust2
find_sample_column = function (metadata, abundance_samples) {
  abundance_samples <- as.character(abundance_samples)
  standard_names <- c("sample", "Sample", "SAMPLE", "sample_name", 
                      "Sample_Name", "SampleName", "samplename", "sample_id", 
                      "Sample_ID", "SampleID", "sampleid", "samples", "Samples")
  for (col in standard_names) {
    if (col %in% colnames(metadata)) {
      if (samples_match(metadata[[col]], abundance_samples)) {
        return(col)
      }
    }
  }
  for (col in colnames(metadata)) {
    if (samples_match(metadata[[col]], abundance_samples)) {
      return(col)
    }
  }
  if (samples_match(rownames(metadata), abundance_samples)) {
    return(".rownames")
  }
  return(NULL)
}

# Yet another helper function from ggpicrust2
samples_match = function (vec1, vec2, threshold = 0.5) {
  vec1 <- as.character(vec1)
  vec2 <- as.character(vec2)
  n_common <- length(intersect(vec1, vec2))
  min_length <- min(length(vec1), length(vec2))
  if (min_length == 0) 
    return(FALSE)
  n_common/min_length >= threshold
}
