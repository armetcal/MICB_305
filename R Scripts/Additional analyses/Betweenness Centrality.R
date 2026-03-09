library(tidyverse)
library(phyloseq)

# Betweenness is a measure of how many shortest paths between other nodes in the network pass through a given node.
# Analogy: if looking at the network of transit stops in Vancouver, the stop with the highest betweenness would probably be 
# Commercial-Broadway or Waterfront. Even if that's not the final stop for most people, those are the 
# stops that transmit the most people/information.
# If using this for your final project, please look the concept up yourself in addition to the above.

# install.packages("igraph")
library(igraph)

# First I'm going to transform the counts into relative abundance.
# Then I'm going to keep only the top 10 most abundant taxa.
# This is an arbitrary choice so that the network is easier to visualize, 
# but you should treat this like the core microbiome analysis and 
# choose abundance/prevalence thresholds.

# Load object
ps = readRDS('Datasets/phyloseq_taxonomy.rds')

# Extract the top 10 taxa (just for demonstration)
avg_abundance = taxa_sums(ps)/sum(taxa_sums(ps)) 
# Sort high to low
avg_abundance = sort(avg_abundance, decreasing = T)
# Take the top 10
top_10 = avg_abundance[1:10]
# Extract taxa names 
top_10 = names(top_10)

# Normalize microbiome data and select only the significant taxa
# CLR transformation
ps_rel = ps %>% microbiome::transform('compositional') 
# Filter taxa
ps_filt = prune_taxa(top_10,ps_rel)

# Next, we're going to extract the abundance data into a table.
# Then we can correlate the taxa with each other.
otu = data.frame(ps_filt@otu_table)
# The cor() command automatically runs Pearson correlations between all the columns in the dataset.
# We want to see how our TAXA correlate with each other, but the taxa are rows.
# We'll therefore transpose our data with t() before running cor().
otu_cor = cor(t(otu))
View(otu_cor) # 10x10 grid of the correlations between taxa

# BUILD THE NETWORK~~~~~~~~~~~~~~~~~

# Decide on a correlation cutoff.
# Pearson coefficients: 0 means no correlation, 1 means perfect positive correlation, -1 means perfect negative correlation.
cor_cutoff = 0.3 # this is a reasonably strong correlation

# VERSION 1: UNWEIGHTED MATRIX~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# In this version, all correlations that are stronger than
# cor_cutoff are treated equally. This is therefore the simpler approach.

# Create an adjacency matrix, which is a table of 0s and 1s that 
# indicates whether each pair of taxa is connected (1) or not (0).
# Connected means that the absolute value of the coefficient is higher than the cutoff.
adj_matrix = ifelse(abs(otu_cor) > cor_cutoff, 1, 0)
# We don't want taxa to be connected to themselves, so we set the diagonal to 0.
diag(adj_matrix) = 0

# Create igraph object (undirected)
g = graph_from_adjacency_matrix(adj_matrix, mode = "undirected", diag = FALSE)

# Calculate betweenness (unweighted)
btw = betweenness(g, normalized = TRUE)

# Result: named numeric vector
btw %>% sort
# 154709e160e8cada6bfb21115acc80f5 has the highest betweenness score.

# VERSION 2: WEIGHTED MATRIX ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# In this version, the strength of the correlation is taken into account.

# Use absolute correlations as weights, set diag zero
wmat = abs(otu_cor)
diag(wmat) = 0
# Remove small weights
wmat[wmat < cor_cutoff] = 0

# Create igraph with weights
g_w = graph_from_adjacency_matrix(wmat, mode = "undirected", weighted = TRUE, diag = FALSE)

# now we have a network similar to version 1 - all correlations above 0.3 are maintained.
# The only difference is that in this version, we can see how strong each connection is
# (think of it like the thickness of a line connecting two taxa).

# However, by default, the graph_from_adjacency_matrix function treats our correlation coefficients
# as distances, meaning that bigger numbers are assumed to be LESS correlated.
# To fix this, we'll simply take 1/X for each value of the correlation matrix.

# Convert to PROPER distances: distance = 1 / weight
# pmax(vector,B) takes the maximum of each value in the vector and B. This is a way to avoid dividing by zero, which would give us infinite distances.
E(g_w)$weight = 1 / pmax(E(g_w)$weight, .Machine$double.eps)

# Compute betweenness using weights as distances
btw_w <- betweenness(g_w, weights = E(g_w)$weight, normalized = TRUE)

btw_w %>% sort # now a different taxon is the most important!

# PLOTTING OPTIONS~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# These are just examples

# OPTION 1: HEATMAP ~~~~~~~~~~~~~~~~~~~~~~
# Packages
install.packages(c("pheatmap"))
library(pheatmap)

# Assume cor_mat is your Spearman correlation matrix with row/col names
# Optionally reorder with clustering, or set row/col order manually
# note that for this to be publication-ready, you'd probably have to change the ASV names (they're too long)
pheatmap(otu_cor,
         clustering_method = "complete",
         color = colorRampPalette(c("blue","white","red"))(50),
         breaks = seq(-1, 1, length.out = 51),
         main = "Pearson's Correlation Heatmap")

# OPTION 2: NETWORK PLOT ~~~~~~~~~~~~~~~~~~~~~~
install.packages(c("ggraph","tidygraph","viridis"))
library(tidygraph); library(ggraph); library(viridis)

# Convert adjacency/weight matrix to tbl_graph
# Reuse wmat from above
g_tbl = as_tbl_graph(g_w)

# Just for fun, add centrality degree (number of connections per taxon)
# Note that this is a bit different from betweenness centrality. 
# Don't feel like you need to do all of these for your final project, 
# I'm only including multiple types of measurements so you know what's available.
g_tbl = g_tbl %>% activate(nodes) %>% mutate(deg = centrality_degree())

# Plot with ggraph
set.seed(421)
ggraph(g_tbl, layout = "fr") +
  geom_edge_link(aes(width = weight), alpha = 0.7, colour = "grey50") +
  geom_node_point(aes(size = deg), color = "tomato") +
  geom_node_text(aes(label = name), repel = TRUE, size = 3) +
  scale_edge_width(range = c(0.2, 2.5)) +
  theme_void() +
  ggtitle("Correlation Network (weighted)")
