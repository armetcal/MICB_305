pathway_errorbar_fixed = function (abundance, daa_results_df, Group, ko_to_kegg = FALSE, 
                            wrap_label = F, wraplength=20, fc_cutoff = 0, order_by_log = F,
                            p_values_threshold = 0.05, order = "group", select = NULL, 
                            p_value_bar = TRUE, colors = NULL, x_lab = NULL){
  
  abundance = ko_filt;
  daa_results_df = daa_annotated_results_df; order_by_log = T;
  Group = lpalmvstongue$body.site; wrap_label = T;wraplength=50; fc_cutoff = 0;
  p_values_threshold = 5e-11; order = "pathway_class"; select = NULL; ko_to_kegg = T; p_value_bar = TRUE;
  colors = NULL; x_lab = "pathway_name"
  
  require(stringr)
  
  ##################### Calculate log 2 fold changes first so we can filter them
  errorbar_abundance_mat <- as.matrix(abundance)
  pseudo = min(abundance[abundance>0])/2
  relative_abundance_mat <- apply(t(errorbar_abundance_mat), 
                                  1, function(x) (x+pseudo)/sum(x)) # Applying a pseudocount of 1
  
  sub_relative_abundance_mat <- relative_abundance_mat[rownames(relative_abundance_mat) %in% 
                                                         (daa_results_df %>% filter(!is.na(p_adjust),p_adjust<p_values_threshold) %>% pull(feature)), ,drop=F]
  
  error_bar_matrix <- cbind(sample = colnames(sub_relative_abundance_mat), 
                            group = Group, 
                            as.data.frame(t(sub_relative_abundance_mat)))
  error_bar_df <- as.data.frame(error_bar_matrix)
  error_bar_df$group <- factor(Group, levels = levels(as.factor(Group)))
  error_bar_pivot_longer_df <- tidyr::pivot_longer(error_bar_df, 
                                                   -c(sample, group)) %>% as_tibble()
  error_bar_pivot_longer_tibble_summarised <- error_bar_pivot_longer_df %>% 
    group_by(name, group) %>% summarise(mean = mean(value)) %>% 
    ungroup() %>% 
    group_by(name) %>% arrange(group) %>% 
    summarise(log2FC = log2(mean[1L]/mean[2L]))
  
  passes_fc = error_bar_pivot_longer_tibble_summarised %>% 
    filter(abs(log2FC)>fc_cutoff) %>% pull(name)
  
  log_order = error_bar_pivot_longer_tibble_summarised %>% 
    arrange(log2FC) %>% pull(name)
  
  #####################
  
  missing_pathways <- daa_results_df[is.na(daa_results_df$pathway_name), 
                                     "feature"]
  if (length(missing_pathways) > 0) {
    message("The following pathways are missing annotations and have been excluded: ", 
            paste(missing_pathways, collapse = ", "))
    message("You can use the 'pathway_annotation' function to add annotations for these pathways.")
  }
  column_names <- colnames(daa_results_df)
  group_columns <- grepl("^group", column_names)
  if (sum(group_columns) > 2) {
    p_value_bar <- FALSE
    message("There are more than two 'group' columns in the 'daa_results_df' data frame. As a result, it is not possible to compute the log2 fold values. The 'p_value_bar' has been automatically set to FALSE.")
  }
  daa_results_df <- daa_results_df[!is.na(daa_results_df[, x_lab]), ]
  if (is.null(x_lab)) {
    if (ko_to_kegg == TRUE) {
      x_lab <- "pathway_name"
    } else {
      x_lab <- "description"
    }
    if (is.null(daa_results_df$pathway_name) & is.null(daa_results_df$description)) {
      message("Please utilize the 'pathway_annotation' function to annotate the 'daa_results_df' data frame.")
    }
  }
  if (!(x_lab %in% colnames(daa_results_df))) {
    message("The 'x_lab' you defined does not exist as a column in the 'daa_results_df' data frame.")
  }
  if (nlevels(factor(daa_results_df$method)) != 1) {
    message("The 'method' column in the 'daa_results_df' data frame contains more than one method. Please filter it to contain only one method.")
  }
  if (nlevels(factor(daa_results_df$group1)) != 1 || nlevels(factor(daa_results_df$group2)) != 
      1) {
    message("The 'group1' or 'group2' column in the 'daa_results_df' data frame contains more than one group. Please filter each to contain only one group.")
  }
  if (is.null(colors)) {
    colors <- c("#d93c3e", "#3685bc", "#6faa3e", "#e8a825", 
                "#c973e6", "#ee6b3d", "#2db0a7", "#f25292")[1:nlevels(as.factor(Group))]
  }
  errorbar_abundance_mat <- as.matrix(abundance)
  daa_results_filtered_df <- daa_results_df[daa_results_df$p_adjust < p_values_threshold & 
                                              daa_results_df$feature %in% passes_fc, ]
  if (!is.null(select)) {
    daa_results_filtered_sub_df <- daa_results_filtered_df[daa_results_filtered_df$feature %in% 
                                                             select, ]
  } else {
    daa_results_filtered_sub_df <- daa_results_filtered_df
  }
  if (nrow(daa_results_filtered_sub_df) > 30) {
    message(paste0("The number of features with statistical significance exceeds 30, leading to suboptimal visualization. ", 
                   "Please use 'select' to reduce the number of features.\n", 
                   "Currently, you have these features: ", paste(paste0("\"", 
                                                                        daa_results_filtered_sub_df$feature, "\""), collapse = ", "), 
                   ".\n", "You can find the statistically significant features with the following command:\n", 
                   "daa_results_df %>% filter(p_adjust < 0.05) %>% select(c(\"feature\",\"p_adjust\"))"))
    stop()
  }
  if (nrow(daa_results_filtered_sub_df) == 0) {
    stop("Visualization with 'pathway_errorbar' cannot be performed because there are no features with statistical significance. ", 
         "For possible solutions, please check the FAQ section of the tutorial.")
  }
  relative_abundance_mat <- apply(t(errorbar_abundance_mat), 
                                  1, function(x) x/sum(x))
  
  sub_relative_abundance_mat <- relative_abundance_mat[rownames(relative_abundance_mat) %in% 
                                                         daa_results_filtered_sub_df$feature, ,drop=F]
 
  error_bar_matrix <- cbind(sample = colnames(sub_relative_abundance_mat), 
                            group = Group, t(sub_relative_abundance_mat))
  error_bar_df <- as.data.frame(error_bar_matrix)
  error_bar_df$group <- factor(Group, levels = levels(as.factor(Group)))
  error_bar_pivot_longer_df <- tidyr::pivot_longer(error_bar_df, 
                                                   -c(sample, group))
  error_bar_pivot_longer_tibble <- mutate(error_bar_pivot_longer_df, 
                                          group = as.factor(group))
  error_bar_pivot_longer_tibble$sample <- factor(error_bar_pivot_longer_tibble$sample)
  error_bar_pivot_longer_tibble$name <- factor(error_bar_pivot_longer_tibble$name)
  error_bar_pivot_longer_tibble$value <- as.numeric(error_bar_pivot_longer_tibble$value)
  error_bar_pivot_longer_tibble_summarised <- error_bar_pivot_longer_tibble %>% 
    group_by(name, group) %>% summarise(mean = mean(value), 
                                        sd = stats::sd(value))
  error_bar_pivot_longer_tibble_summarised <- error_bar_pivot_longer_tibble_summarised %>% 
    mutate(group2 = "nonsense")
  switch(order, p_values = {
    order <- order(daa_results_filtered_sub_df$p_adjust)
  }, name = {
    order <- order(daa_results_filtered_sub_df$feature)
  }, group = {
    daa_results_filtered_sub_df$pro <- 1
    for (i in levels(error_bar_pivot_longer_tibble_summarised$name)) {
      error_bar_pivot_longer_tibble_summarised_sub <- error_bar_pivot_longer_tibble_summarised[error_bar_pivot_longer_tibble_summarised$name == i, ]
      pro_group <- error_bar_pivot_longer_tibble_summarised_sub[error_bar_pivot_longer_tibble_summarised_sub$mean == 
                                                                  max(error_bar_pivot_longer_tibble_summarised_sub$mean),]$group
      pro_group <- as.vector(pro_group)
      daa_results_filtered_sub_df[daa_results_filtered_sub_df$feature == i, ]$pro <- pro_group
    }
    order <- order(daa_results_filtered_sub_df$pro, daa_results_filtered_sub_df$p_adjust)
  }, pathway_class = {
    if (!"pathway_class" %in% colnames(daa_results_filtered_sub_df)) {
      stop("The 'pathway_class' column is missing in the 'daa_results_filtered_sub_df' data frame. ", 
           "Please use the 'pathway_annotation' function to annotate the 'pathway_daa' results.")
    }
    order <- order(daa_results_filtered_sub_df$pathway_class, 
                   daa_results_filtered_sub_df$p_adjust)
  }, {
    order <- order
  })
  daa_results_filtered_sub_df <- daa_results_filtered_sub_df[order, ]
  error_bar_pivot_longer_tibble_summarised_ordered <- data.frame(name = NULL, 
                                                                 group = NULL, mean = NULL, sd = NULL)
  for (i in daa_results_filtered_sub_df$feature) {
    error_bar_pivot_longer_tibble_summarised_ordered <- rbind(error_bar_pivot_longer_tibble_summarised_ordered, 
                                                              error_bar_pivot_longer_tibble_summarised[error_bar_pivot_longer_tibble_summarised$name == i, ])
  }
  if (ko_to_kegg == FALSE) {
    error_bar_pivot_longer_tibble_summarised_ordered[, x_lab] <- rep(daa_results_filtered_sub_df[, x_lab], each = length(levels(factor(error_bar_pivot_longer_tibble_summarised_ordered$group))))
  }
  if (ko_to_kegg == TRUE) {
    error_bar_pivot_longer_tibble_summarised_ordered$pathway_class <- rep(daa_results_filtered_sub_df$pathway_class, 
                                                                          each = length(levels(factor(error_bar_pivot_longer_tibble_summarised_ordered$group))))
  }

  if(order_by_log == F){
    error_bar_pivot_longer_tibble_summarised_ordered$name <- factor(error_bar_pivot_longer_tibble_summarised_ordered$name, 
                                                                    levels = rev(daa_results_filtered_sub_df$feature))
  } else {
    error_bar_pivot_longer_tibble_summarised_ordered$name <- factor(error_bar_pivot_longer_tibble_summarised_ordered$name, 
                                                                    levels = rev(log_order))
  }
  
  
  daa_results_filtered_sub_df$pathway_name = sapply(daa_results_filtered_sub_df$pathway_name, function(x) str_split(x,' \\[')[[1]][1] %>% str_trim())
  if(wrap_label==T){
    daa_results_filtered_sub_df$pathway_name = stringr::str_wrap(daa_results_filtered_sub_df$pathway_name,wraplength)
  }
  
  bar_errorbar <- ggplot2::ggplot(error_bar_pivot_longer_tibble_summarised_ordered, 
                                  ggplot2::aes(mean, name, fill = group)) + 
    ggplot2::geom_errorbar(ggplot2::aes(xmax = mean + sd, xmin = 0), 
                           position = ggplot2::position_dodge(width = 0.8), 
                           width = 0.5, size = 0.5, color = "black") + 
    ggplot2::geom_bar(stat = "identity", 
                                                                                                                                                                 position = ggplot2::position_dodge(width = 0.8), width = 0.8) + 
    GGally::geom_stripped_cols(width = 10) + ggplot2::scale_fill_manual(values = colors) + 
    ggplot2::scale_color_manual(values = colors) + ggprism::theme_prism() + 
    ggplot2::scale_x_continuous(expand = c(0, 0))+ #, guide = "prism_offset_minor") + 
    ggplot2::scale_y_discrete(labels = rev(daa_results_filtered_sub_df[, x_lab])) + 
    ggplot2::labs(x = "Relative Abundance", y = NULL) + 
    ggplot2::theme(axis.ticks.y = ggplot2::element_blank(), 
                   axis.line.y = ggplot2::element_blank(), axis.line.x = ggplot2::element_line(size = 0.5), 
                   axis.ticks.x = ggplot2::element_line(size = 0.5), 
                   panel.grid.major.y = ggplot2::element_blank(), panel.grid.major.x = ggplot2::element_blank(), 
                   axis.text = ggplot2::element_text(size = 13, color = "black"), 
                   axis.text.x = ggplot2::element_text(margin = ggplot2::margin(r = 0)), 
                   axis.text.y = ggplot2::element_text(size = 13, color = "black", 
                                                       margin = ggplot2::margin(b = 6)), 
                   axis.title.x = ggplot2::element_text(size = 13, color = "black", hjust = 0.5), 
                   legend.position = "top", 
                   legend.key.size = ggplot2::unit(0.3, "cm"), 
                   legend.direction = "vertical", 
                   legend.justification = "left", 
                   legend.text = ggplot2::element_text(size = 13, face = "bold"), 
                   legend.box.just = "right", plot.margin = ggplot2::margin(0, 0.5, 0.5, 0, unit = "cm")) + 
    ggplot2::coord_cartesian(clip = "off")
  
  if (ko_to_kegg == TRUE) {
    pathway_class_group_mat <- daa_results_filtered_sub_df$pathway_class %>% 
      table() %>% data.frame() %>% column_to_rownames(".")
    pathway_class_group <- data.frame(. = unique(daa_results_filtered_sub_df$pathway_class), 
                                      Freq = pathway_class_group_mat[unique(daa_results_filtered_sub_df$pathway_class), 
                                      ])
    start <- c(1, rev(pathway_class_group$Freq)[1:(length(pathway_class_group$Freq) - 
                                                     1)]) %>% cumsum()
    end <- cumsum(rev(pathway_class_group$Freq))
    ymin <- start - 1/2
    ymax <- end + 1/2
    nPoints <- length(start)
    pCol <- c("#D51F26", "#272E6A", "#208A42", "#89288F", 
              "#F47D2B", "#FEE500", "#8A9FD1", "#C06CAB", "#E6C2DC", 
              "#90D5E4", "#89C75F", "#F37B7D", "#9983BD", "#D24B27", 
              "#3BBCA8", "#6E4B9E", "#0C727C", "#7E1416", "#D8A767", 
              "#3D3D3D")[1:nPoints]
    pFill <- pCol
    for (i in 1:nPoints) {
      bar_errorbar <- bar_errorbar + ggplot2::annotation_custom(
        grob = grid::rectGrob(gp = grid::gpar(col = pCol[i],fill = pFill[i], lty = NULL, lwd = NULL, alpha = 0.2)),
        xmin = ggplot2::unit(-2, "native"), xmax = ggplot2::unit(0, "native"), ymin = ggplot2::unit(ymin[i], "native"), 
                                                                ymax = ggplot2::unit(ymax[i], "native"))
    }
  }
  daa_results_filtered_sub_df <- cbind(daa_results_filtered_sub_df, 
                                       negative_log10_p = -log10(daa_results_filtered_sub_df$p_adjust), 
                                       group_nonsense = "nonsense", log_2_fold_change = NA)
  for (i in daa_results_filtered_sub_df$feature) {
    mean <- error_bar_pivot_longer_tibble_summarised_ordered[error_bar_pivot_longer_tibble_summarised_ordered$name %in% 
                                                               i, ]$mean
    daa_results_filtered_sub_df[daa_results_filtered_sub_df$feature == 
                                  i, ]$log_2_fold_change <- log2(mean[1]/mean[2])
  }
  if(order_by_log == F) {
    daa_results_filtered_sub_df$feature <- factor(daa_results_filtered_sub_df$feature, 
                                                  levels = rev(daa_results_filtered_sub_df$feature))
  } else {
    daa_results_filtered_sub_df$feature <- factor(daa_results_filtered_sub_df$feature, 
                                                  levels = rev(log_order))
  }

  p_values_bar <- daa_results_filtered_sub_df %>% 
    ggplot2::ggplot(ggplot2::aes(feature, log_2_fold_change, fill = group_nonsense)) + 
    ggplot2::geom_bar(stat = "identity", position = ggplot2::position_dodge(width = 0.8), width = 0.8) + 
    ggplot2::labs(y = "Log2 Fold Change", x = NULL) + GGally::geom_stripped_cols() + 
    ggplot2::scale_fill_manual(values = "#87ceeb") + ggplot2::scale_color_manual(values = "#87ceeb") + 
    ggplot2::geom_hline(ggplot2::aes(yintercept = 0), linetype = "dashed", 
                        color = "black") + 
    ggprism::theme_prism() +
    # ggplot2::scale_y_continuous(expand = c(0, 0), guide = "prism_offset_minor") + 
    ggplot2::theme(axis.ticks.y = ggplot2::element_blank(),
                   axis.line.y = ggplot2::element_blank(),
                   axis.line.x = ggplot2::element_line(size = 0.5), 
                   axis.ticks.x = ggplot2::element_line(size = 0.5), 
                   panel.grid.major.y = ggplot2::element_blank(), 
                   panel.grid.major.x = ggplot2::element_blank(),
                   axis.text = ggplot2::element_text(size = 13,color = "black"), 
                   axis.text.y = ggplot2::element_blank(), 
                   axis.text.x = ggplot2::element_text(size = 13, color = "black", margin = ggplot2::margin(b = 6)), 
                   axis.title.x = ggplot2::element_text(size = 14, color = "black", hjust = 0.5), legend.position = "non") + 
    ggplot2::coord_flip()
  if (ko_to_kegg == TRUE) {
    pathway_class_y <- (ymax + ymin)/2 - 0.5
    pathway_class_plot_df <- as.data.frame(cbind(nonsense = "nonsense", 
                                                 pathway_class_y = pathway_class_y, 
                                                 pathway_class = rev(unique(daa_results_filtered_sub_df$pathway_class))))
    pathway_class_plot_df$pathway_class = sapply(pathway_class_plot_df$pathway_class, function(x) str_split(x,';')[[1]] %>% .[length(.)] %>% str_trim())
    pathway_class_plot_df$pathway_class = sapply(pathway_class_plot_df$pathway_class, function(x) str_remove_all(x,'KEGG Orthology [(]KO[)] \\[BR:ko00001\\];  '))

    pathway_class_plot_df$pathway_class_y <- as.numeric(pathway_class_plot_df$pathway_class_y)
    if(nrow(pathway_class_plot_df)==2){ 
      pathway_class_plot_df$pathway_class_y = 1; pathway_class_plot_df = pathway_class_plot_df[1,,drop=F]
      }
    pathway_class_annotation <- pathway_class_plot_df %>% 
      ggplot2::ggplot(ggplot2::aes(nonsense, pathway_class_y)) + 
      ggplot2::geom_text(ggplot2::aes(nonsense, pathway_class_y, 
                                      label = pathway_class), size = 4.5, color = "black", 
                         fontface = "bold", family = "sans") + ggplot2::scale_y_discrete(position = "right") + 
      ggprism::theme_prism() + ggplot2::theme(axis.ticks = ggplot2::element_blank(), 
                                              axis.line = ggplot2::element_blank(), panel.grid.major.y = ggplot2::element_blank(), 
                                              panel.grid.major.x = ggplot2::element_blank(), panel.background = ggplot2::element_blank(), 
                                              axis.text = ggplot2::element_blank(), plot.margin = ggplot2::unit(c(0, 0.2, 0, 0), "cm"), 
                                              axis.title.y = ggplot2::element_blank(), 
                                              axis.title.x = ggplot2::element_blank(), legend.position = "non")
  }
  daa_results_filtered_sub_df$p_adjust <- as.character(signif(daa_results_filtered_sub_df$p_adjust,2))
  daa_results_filtered_sub_df$unique <- nrow(daa_results_filtered_sub_df) - 
    seq_len(nrow(daa_results_filtered_sub_df)) + 1

  p_annotation <- daa_results_filtered_sub_df %>% 
    ggplot2::ggplot(ggplot2::aes(group_nonsense,p_adjust)) + 
    ggplot2::geom_text(ggplot2::aes(group_nonsense, unique, label = p_adjust), size = 5, 
                       color = "black",fontface = "bold", family = "sans") +
    ggplot2::labs(y = "p-value (adjusted)") + 
    ggplot2::scale_y_discrete(position = "right") + ggprism::theme_prism() + 
    ggplot2::theme(axis.ticks = ggplot2::element_blank(), 
                   axis.line = ggplot2::element_blank(), panel.grid.major.y = ggplot2::element_blank(), 
                   panel.grid.major.x = ggplot2::element_blank(), panel.background = ggplot2::element_blank(), 
                   axis.text = ggplot2::element_blank(), 
                   plot.margin = ggplot2::unit(c(0,0.2, 0, 0), "cm"), 
                   axis.title.y = ggplot2::element_text(size = 14, color = "black", vjust = 0), 
                   axis.title.x = ggplot2::element_blank(), 
                   legend.position = "none")
  if (p_value_bar == TRUE) {
    if (ko_to_kegg == TRUE) {
      combination_bar_plot <- pathway_class_annotation + 
        bar_errorbar + p_values_bar + p_annotation + 
        patchwork::plot_layout(ncol = 4, widths = c(1.2, 
                                                    1.0, 0.5, 0.3))
    } else {
      combination_bar_plot <- bar_errorbar + p_values_bar + 
        p_annotation + patchwork::plot_layout(ncol = 3, 
                                              widths = c(2.3, 0.7, 0.5))
    }
  } else {
    if (ko_to_kegg == TRUE) {
      combination_bar_plot <- pathway_class_annotation + 
        bar_errorbar + p_annotation + patchwork::plot_layout(ncol = 3, 
                                                             widths = c(1, 1.2, 0.3))
    } else {
      combination_bar_plot <- bar_errorbar + p_annotation + 
        patchwork::plot_layout(ncol = 2, widths = c(2.5, 
                                                    0.4))
    }
  }
  return(combination_bar_plot)
}