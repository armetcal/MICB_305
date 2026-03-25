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
  aligned = ggpicrust2:::align_samples(abundance, metadata, verbose = FALSE)
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