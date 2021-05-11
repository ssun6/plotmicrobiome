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
```
taxa_table=format_asv(taxa_file ="biom-with-taxonomy.biom",biom=T,onefile = T,ASV=T)
```
For shotgun metagenome sequencing data, plotmicrobiome supports the output of the most commonly used
software MetaPhlAn and Kraken. 
```
taxa_table=format_wgs(taxa_file = "kraken2.txt")
taxa_table=format_wgs(taxa_file = "metaphlan2.txt")
```

### Metadata input
Metadata should include sample IDs that are consistent with Data input files and variables for testing.
**App examples:**




**R package Example:**
```
metadata=meta_format(metadata=metadata_dir,metadata_sep="\t",meta_sample_name_col=1)
```
### MDS plot
Plotmicrobiome can generate both PCoA and non-parametric NMDS plots with selected dissimilarity measures.
PERMANOVA test results of the association between microbiome and selected variable were shown as headers. 
**App examples:**

**R package Example:**
```
mds_plot(taxa_table = taxa_table, metadata=metadata,test_metadata="group",method_mds = "pcoa",palette_group=c("red","blue","orange","green"),distance_type="bray")

```

### Alpha-diversity plot
Plotmicrobiome provides test results and visualization of the alpha-diversity difference between groups at each taxonomic level.
The metrics analyzed include Shannon index, Simpson index, Inverse Simpson index and numbers of taxa. 

**App examples:**

**R package Example:**
```
alpha_plot(taxa_table = taxa_table, metadata=metadata,test_metadata="group",palette_group=c("red","blue","orange","green"))
```

### Data filter
The data can be stratified by metadata to only include samples that are used for statistical analysis. 
Abundance and Prevalence cutoffs can be used to exclude taxa of low abundance and/or low prevalence.
**App examples:**

**Example:**
``
tab_s=table_subset(taxa_table = taxa_table,metadata=metadata,stratify_by_metadata="Study",stratify_by_value="Sugar")

``
### Statistical tests
Statistical tests can be used to determine differential taxa between groups. 
Plotmicrobiome supports t-test, wilcoxon test, ANOVA and kruskal-wallis tests for now and more tests will be added in new updates.

**App examples:**

**Example:**
``
fdrs1=stat_test(taxa_table =tab_s,metadata=metadata,test_metadata="group",method="wilcoxon")

``
### Boxplots
Boxplots are used to visualize the differetial abundance of taxa identified in statistical tests.

**App examples:**

**Example:**
``
taxa_boxplot(taxa_table = tab_s, metadata=metadata1,test_metadata="group",fdrs=fdrs1,log_norm=T,cutoff=0.01,palette_group=c("red","blue","orange","green"))

``

### Tree plots
Significant taxa are highlighted in a taxonomic tree.

**App examples:**

**Example:**
``
plot1=tree_view(taxa_table =tab_s, metadata=metadata1,fdrs=fdrs1,test_metadata="group",prevalence_cutoff=0.1, abundance_cutoff=0)

``

### Correlation plots
The associations between taxa and continuous variables can be tested and visualized in this step.

**App examples:**

**Example:**
``
cor_plot1=meta_corplot(taxa_table =tab_s, metadata=metadata,test_metadata="test_score",col_metadata="group",fdr_cutoff=0.1)

``
