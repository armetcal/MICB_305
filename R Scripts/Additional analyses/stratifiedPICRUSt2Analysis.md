# Stratified PICRUSt2 Analysis Script

Author: IY

Created: May 3, 2026

Last modified: May 3, 2026

---

The following Bash script is for a stratified PICRUSt2 analysis. This is relatively similar to the typical input used for regular (unstratified) PICRUSt2 analyses,
only with the addition of the ```--stratified``` flag to the command. The full process will be noted here for clarity.

## Data Preparation

Prepare the feature table as you would for unstratified PICRUSt2 analyses by filtering it for the desired number of features.
Use your desired minimum frequency in ```DESIRED_FREQUENCY```. This is a number. Talk to your TA about what frequency you should filter your dataset to.

```
# Import and filter your feature table
qiime feature-table filter-features \
  --i-table YOUR_FILTERED_TABLE.qza \
  --p-min-frequency DESIRED_FREQUENCY \
  --o-filtered-table OUTPUT_FREQ_FILTERED_FEATURE_TABLE.qza # Make sure to give it a memorable name!
```
As this is the same for unstratified PICRUSt2 preparation, you can alternatively reuse your filtered feature table for your stratified PICRUSt2 analysis.
```
# Make the stratified PICRUSt2 directory and export the data
mkdir strat_picrust

qiime tools export \
   --input-path OUTPUT_FREQ_FILTERED_FEATURE_TABLE.qza \ # This is your new file from the previous command
   --output-path strat_picrust

qiime tools export \
   --input-path YOUR_REP_SEQS_FILE.qza \ # Refer back to the QIIME2 modules for help generating this file
   --output-path strat_picrust

# Switch from the Qiime2 to the PICRUSt2 environment

conda deactivate
conda activate picrust2
```

## Stratified PICRUSt2 Analysis

Now we will run PICRUSt2 for a stratified analysis. Make sure you are in the PICRUSt2 environment before proceeding!

```
picrust2_pipeline.py \
-s strat_picrust/dna-sequences.fasta \
-i strat_picrust/feature-table.biom \
--stratified \ # The stratified flag is appended to the command for the stratified analysis to run
-o picrust_strat_output
```
Now you can use R to visualize your stratified PICRUSt2 analysis.

Happy coding!

> Made with ❤️ from team 7 (2025W2) 🐟🐠🐡
