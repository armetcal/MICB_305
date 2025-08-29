library(devtools)
install_github('biobakery/MTX_model')

library(MTXmodel) #add the library to the working path

#read-in each file from the package then convert it to a data.frame
### RNA abundances
input_data <- system.file(
  'extdata','HMP2_pwyRNA.tsv', package="MTXmodel")

df_input_data = read.table(file             = input_data,
                           header           = TRUE,
                           sep              = "\t", 
                           row.names        = 1,
                           stringsAsFactors = FALSE)
df_input_data[1:5, 1:5]

# RNA/DNA ratio data 
#input_dataratio <- system.file(
#  'extdata','HMP2_pwy.RNA_DNA_ratio.tsv', package="MTXmodel")
df_input_dataratio = read.table(file             = "~/Practice Repo/HMP2_pwy.RNA_DNA_ratio.tsv",
                                header           = TRUE,
                                sep              = "\t", 
                                row.names        = 1,
                                stringsAsFactors = FALSE)
df_input_dataratio[1:5, 1:5]

# Metadata from the HMP2
input_metadata <-system.file(
  'extdata','HMP2_metadata.tsv', package="MTXmodel")
df_input_metadata = read.table(file             = input_metadata,
                               header           = TRUE,
                               sep              = "\t", 
                               row.names        = 1,
                               stringsAsFactors = FALSE)
df_input_metadata[1:5, 1:5]

# DNA data 
input_dnadata <- system.file(
  'extdata','HMP2_pwyDNA.tsv', package="MTXmodel")
df_input_dnadata = read.table(file             = input_dnadata,
                              header           = TRUE,
                              sep              = "\t", 
                              row.names        = 1,
                              stringsAsFactors = FALSE)
df_input_dnadata[1:5, 1:5]


# Raw RNA Abundances with MaAsLin 2
library(Maaslin2)
fit_rna <- Maaslin2(
  df_input_data, df_input_metadata, 'demo_output_rna', transform = "LOG",
  fixed_effects = c('diagnosis', 'dysbiosisnonIBD','dysbiosisUC','dysbiosisCD', 'antibiotics', 'age'),
  random_effects = c('site', 'subject'),
  reference = "diagnosis,nonIBD",
  normalization = 'NONE',
  standardize = FALSE
)

# RNA/DNA Ratios with MaAsLin 2
fit_rna_ratio <- Maaslin2(
  df_input_dataratio, df_input_metadata, 'demo_output_rna_ratio', transform = "LOG",
  fixed_effects = c('diagnosis', 'dysbiosisnonIBD','dysbiosisUC','dysbiosisCD', 'antibiotics', 'age'),
  random_effects = c('site', 'subject'),
  reference = "diagnosis,nonIBD",
  normalization = 'NONE',
  standardize = FALSE
)

# RNA abundance adjusted by DNA abundance with MTXmodel
library(MTXmodel)
fit_rna_dna <- MTXmodel(
  df_input_data, df_input_metadata, 'demo_output_rna_dna', transform = "LOG",
  fixed_effects = c('diagnosis', 'dysbiosisnonIBD','dysbiosisUC','dysbiosisCD', 'antibiotics', 'age'),
  random_effects = c('site', 'subject'),
  reference = "diagnosis,nonIBD",
  normalization = 'NONE',
  standardize = FALSE,
  input_dnadata = df_input_dnadata
)

# Compare output
#compare the raw counts of features associated with CD
results_rna = subset(fit_rna$results, metadata == "diagnosis")
table(results_rna$value)
#CD/UC = 186

results_rna_ratio = subset(fit_rna_ratio$results, metadata == "diagnosis")
table(results_rna_ratio$value)
#CD/UC = 37

results_rna_dna = subset(fit_rna_dna$results, metadata == "diagnosis")
table(results_rna_dna$value)
#CD/UC = 178

# features called by the RNA/DNA ratios compared to the Raw RNA abundances 
intersect(results_rna_ratio$feature, results_rna$feature) # 37 features overlapped 

# features called by the Raw RNA abundances compared to the RNA abundances adjusted by the DNA abundances 
intersect(results_rna$feature, results_rna_dna$feature) #178 features overlapped between the models

# features called by the RNA/DNA ratios compared to the RNA abundances adjusted by the DNA abundances 
intersect(results_rna_ratio$feature, results_rna_dna$feature) # 37 features overlapped 

#compare the raw counts of features associated with CD
head(results_rna$feature, 3)
results_rna$model = "RNA" # create a column that describes the data

#compare the raw counts of features associated with CD
head(results_rna_ratio$feature, 3)
results_rna_ratio$model = "RNA/DNA Ratio" # create a column that describes the data

#compare the raw counts of features associated with CD
head(results_rna_dna$feature, 3)
results_rna_dna$model = "RNA adjusted by DNA" # create a column that describes the data
list = unique(results_rna_dna$feature)[c(1:10)] # select the top pathway results for the next step

#plot top results
library(ggplot2)
results = rbind(results_rna, results_rna_dna, results_rna_ratio) # create one object from the results 
results = results[results$feature %in% list, ] # subset to just the top features in the MTXmodel
ggplot(results, aes(x = coef, y = feature, color = model)) + geom_point(aes(shape = value)) + theme_classic()
ggsave("model_comparison.png", height = 5, width = 8) # save the plot to the working directory 