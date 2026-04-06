# This script calculates how much each taxon contributes to each pathway.
# The idea is to figure out if there are any taxa that contribute a large 
#   proportion of the reads to a given pathway, which would suggest that they 
#   are the main drivers of that pathway in the community.
# If there are no clear high-contributing taxa, that would suggest that the 
#   pathway is more evenly distributed across the community, and there are no 
#   clear 'drivers' of that pathway.

# Credit: Izaak Yip generated the stratified PICRUSt2 data used here.
#         Sam Donato drafted the first version of this script. 
#         Dr. Avril Metcalfe-Roach completed the second version of this script. 
#         Sam Donato visualized the taxonomic bar plot. (March 2026)

# Stratified PICRUSt2 reference: https://github-wiki-see.page/m/picrust/picrust2/wiki/PICRUSt2-Tutorial-(v2.5.0)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# READ THE CODE CAREFULLY so you know what's going on! It can be a little confusing.

# Talk with your TA about the best way to show the data. 
#  A taxonomic bar plot-type figure might work well, where each column is one 
#    function and the colours show the proportion of reads from each taxon.
#  Or a table to show how many pathways each key taxa contribute significantly to.
#    (significance: use some cutoff, such as a taxon has to contribute 33% of 
#     the reads for a given pathway)
#  Or many other approaches.

# Load packages 
library(tidyverse)
library(readr)
library(dplyr)
library(ggplot2)
library(phyloseq)

# Define Pathways of interest
energy_pathways <- c(
  "ANAGLYCOLYSIS-PWY",
  "GLYCOLYSIS",
  "GLYCOLYSIS-E-D",
  "FERMENTATION-PWY",
  "METH-ACETATE-PWY",
  "FAO-PWY",
  "COA-PWY",
  "COA-PWY-1",
  "GLUCONEO-PWY",
  "PENTOSE-P-PWY",
  "NONOXIPENT-PWY",
  "GLYOXYLATE-BYPASS"
)

# PART 1 - Getting the data ready
# * Aggregate to the desired taxonomic level 
#   CURRENTLY FAMILY LEVEL - WILL NEED TO CHANGE THROUGHOUT IF YOU WANT A DIFFERENT LEVEL!!!
# * Convert to relative abundance to control for sequencing depth

# Load stratified pathway contribution data (MetaCyc used for this example)
strat_mc_data_1 <- data.table::fread("Datasets/STRATIFIED_PICRUST_EXAMPLE_path_abun_contrib.tsv")

# Convert reads per sample to relative abundance to control for sequencing depth
strat_rel = strat_mc_data_1 %>% 
  group_by(sample) %>%
  mutate(rel_abun = taxon_function_abun/sum(taxon_function_abun)) %>%
  ungroup()

head(strat_rel)

# Only keep the pathways of interest
strat_filt = strat_rel %>% 
  # Only include pathways of interest
  filter(`function` %in% energy_pathways)

head(strat_filt)

# Aggregate reads to the Family level
# Load tax data, keep just the Family info
tax <- read_tsv("Datasets/STRATIFIED_PICRUST_EXAMPLE_taxonomy.tsv") %>%  # Load Taxonomy Data Names
  separate(`Taxon`, into = c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus"), sep = ";") %>%
  select(taxon = `Feature ID`, Family) %>% 
  mutate(Family = ifelse(is.na(Family),taxon, Family)) # If no Family name, keep the feature ID for now

# Aggregate the ASVs to the Family level
strat_Family = strat_filt %>%
  left_join(tax) %>%
  group_by(Family, `function`, sample) %>% 
  # This is the 'amount' of the function contributed by each Family
  summarize(rel_abun = sum(rel_abun)) %>% 
  ungroup()

head(strat_Family)

# PART 2 - Calculating the contributions of each taxon to each function

# Add together the counts from each sample for each taxon and function
strat_avg = strat_Family %>% 
  group_by(Family, `function`) %>% 
  # combine all the relative abundances for a given taxon and pathway across all samples to get an average contribution of that taxon to that pathway
  summarize(total_contrib_per_Family = sum(rel_abun)) %>% 
  ungroup()

head(strat_avg)

# Now we just have three columns: the Family, the function, and the function's abundance that belongs to that taxon.
# Next, we figure out the proportion of each function that belongs to each taxon.
strat_prop = strat_avg %>% 
  group_by(`function`) %>%
  mutate(prop_tax_reads_per_function = total_contrib_per_Family/sum(total_contrib_per_Family)) %>% 
  ungroup() %>% 
  select(-total_contrib_per_Family) # We don't need the total contribution per Family anymore, just the proportion

# Sanity check. At this point, the abundances for each function should add up to 1, 
# and the individual values represent proportions.
strat_prop %>% group_by(`function`) %>% 
  summarize(sum = sum(prop_tax_reads_per_function)) %>% 
  ungroup() # Every pathway adds up to 1! Good!

# PART 3 - get a feel for the data

# Filter for taxa that contribute at least 33% of the reads to a given pathway
strat_cutoff = strat_prop %>% 
  filter(prop_tax_reads_per_function >= 0.33) # None!

# Top hits:
head(strat_prop %>% arrange(-prop_tax_reads_per_function))

# Histogram of the proportions, just to see the distributions. 
# No clear high-contributing taxa - a taxonomic bar plot-style figure might work best.
hist(log10(strat_prop$prop_tax_reads_per_function),breaks=50)

# Part 4 - Visualizing taxonomic bar plot 

# Create the 5% cutoff
threshold <- 0.05  

# Group the <5% reads by function and family, sum all the rare taxa 
strat_prop_grouped <- strat_prop %>%
  group_by(`function`) %>%
  mutate(Family = ifelse(prop_tax_reads_per_function < threshold,
                         "<5% of reads",
                         Family)) %>%
  ungroup() %>%
  group_by(`function`, Family) %>%
  summarize(prop_tax_reads_per_function = sum(prop_tax_reads_per_function),
            .groups = "drop")

# Plot without a legend for visibility
plot_avg <- strat_prop_grouped %>%
  mutate(`function` = str_to_title(`function`)) %>%
  ggplot(aes(`function`, prop_tax_reads_per_function, fill = Family)) +
  labs(y = "Proportion of Family taxa reads per function") +
  geom_col(position = "stack", color = NA) +
  theme_classic(base_size = 22) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
       legend.position = "none")+
  xlab(NULL)

plot_avg

# Plot with a legend 
plot_avg <- strat_prop_grouped %>%
  mutate(`function` = str_to_title(`function`)) %>%
  ggplot(aes(`function`, prop_tax_reads_per_function, fill = Family)) +
  labs(y = "Proportion of Family taxa reads per function") +
  geom_col(position = "stack", color = NA) +
  theme_classic(base_size = 22) +
  theme(strip.text = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1),
        legend.title = element_blank(),
        legend.key.size = unit(1, "cm"),
        legend.text = element_text(size = 16)) +
  xlab(NULL)

plot_avg
ggsave("plot_avg.jpg", plot_avg)
