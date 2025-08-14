# if (!require("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# 
# BiocManager::install("phyloseq")
getwd()
# load libraries
library(tidyverse)
library(phyloseq)

# Load datasets
taxonomy = read.delim('Datasets/taxonomy.tsv', row.names = 1)
tree = read_tree('Datasets/tree.nwk')

counts = read.delim('Datasets/feature-table.txt', skip=1, row.names=1) # First line is not data
metadata = read.delim('Datasets/sample-metadata.tsv', row.names = 1) # 1st col are names

# Wrangle Tables ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Taxonomy
taxonomy_formatted = taxonomy %>% 
  separate(col = Taxon, 
           into = c('Domain','Phylum','Class','Order',
                    'Family','Genus','Species'),
           sep=';', fill='right') %>% 
  select(-Confidence) %>% 
  as.matrix()

# Counts
counts_formatted = counts %>% as.matrix()

# Metadata ~~~~~~~
View(metadata)
# Extract just the column names
meta_names = read.delim('Datasets/sample-metadata.tsv',row.names = 1) %>% names()
# Skip first two lines, as they are not data
# The first line is the header, the second is the data type
meta_data = read.delim('Datasets/sample-metadata.tsv',row.names = 1,skip=2,header = F) 
table(metadata[2,]== meta_data[1,]) # All values in the first row of real data are identical

# Set column names for metadata
names(meta_data) = meta_names
# meta_data = meta_data %>% `names<-`(meta_names) # `<-` commands allow for piping

# Create the phyloseq object
ps = phyloseq(sample_data(meta_data),
              otu_table(counts_formatted, taxa_are_rows = T),
              tax_table(taxonomy_formatted),
              tree)

# Save as .rds or .Rdata object
saveRDS(ps,'Datasets/phyloseq_taxonomy.rds')
