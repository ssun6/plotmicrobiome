# plotmicrobiome
Plotmicrobiome is a statistical analysis and visualization pipeline for microbiome analysis. Plotmicrobiome intergrates novel approaches in analyses and visualization of microbiome data in a user-friendly way and generates publication ready statistical results and figures for both 16S rRNA gene amplicon sequencing and shotgun metagenome sequencing. We also developed a interaction application that is freely available at https://ssun6.shinyapps.io/plotmicrobiome/. 

## Installation
**Website**
Plotmicrobiome is freely available at https://ssun6.shinyapps.io/plotmicrobiome/ without installation. 
**R package**
install.packages("devtools")
devtools::install_github("ssun6/plotmicrobiome")

## Tutorials
### Data input
Plotmicrobiome supports the common output files from sequencing analysis pipelines.
For amplicon sequencing data:
1. A single file with all taxonomic levels and in .csv, .tsv and .biom formats.
2. Multiple files with one for each taxonomic level separately
3. DADA2 amplicon sequence variant (ASV) table with taxonomic classification in .csv, .tsv and .biom formats.

**App examples:**


**R package Example:**
``
taxa_table=format_asv(taxa_file ="biom-with-taxonomy.biom",biom=T,onefile = T,ASV=T)
``
For shotgun metagenome sequencing data, plotmicrobiome supports the output of the most commonly used
software MetaPhlAn and Kraken. 
``
taxa_table=format_wgs(taxa_file = "kraken2.txt")
taxa_table=format_wgs(taxa_file = "metaphlan2.txt")
``

### Metadata input
Metadata should include sample IDs that are consistent with Data input files and variables for testing.
**App examples:**




**R package Example:**
``
metadata=meta_format(metadata=metadata_dir,metadata_sep="\t",meta_sample_name_col=1)
``
### MDS plot
Plotmicrobiome can generate both PCoA and non-parametric NMDS plots with selected dissimilarity measures.
PERMANOVA test results of the association between microbiome and selected variable were shown as headers. 
**App examples:**


### Subset of samples
The data can be stratified by metadata to only include samples that are used for statistical analysis. Abundance
and Prevalence cutoffs can be used to exclude 
**Example:**
``
tab_s=table_subset(taxa_table = taxa_tab1,metadata=metadata1,stratify_by_metadata="Study",stratify_by_value="Sugar")

``


