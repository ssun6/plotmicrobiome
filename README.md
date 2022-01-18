# plotmicrobiome
Plotmicrobiome is a user-friendly statistical analysis and visualization pipeline for microbiome analysis. Plotmicrobiome intergrates novel approaches in analyses and visualization of microbiome data in a user-interactive way and generates publication ready statistical results and figures for both 16S rRNA gene amplicon sequencing and shotgun metagenome sequencing. 

# Installation
There are three ways to use the pipeline: 1.GUI with R on your computer, 2. shiny apps website and 3. R package.

## 1. R shiny (recommended)

Run in Command line:
```
cd $HOME
git clone https://github.com/ssun6/plotmicrobiome.git
```

Run in R:
```
library(shiny)
runApp('$HOME/plotmicrobiome')
```
If you want to look at the code:
```
runApp('$HOME/plotmicrobiome', display.mode = "showcase")
```



## 2. Website

Plotmicrobiome is freely available at https://ssun6.shinyapps.io/plotmicrobiome/ without installation. But the website is limited to 1GB memory and 25 active hours per month. 




## 3. R package

```
install.packages("devtools")
devtools::install_github("ssun6/plotmicrobiome")
```



# Tutorials
# Load data
## Getting ready
1. Count table files in .csv, .tsv or .biom formats. For .csv and .tsv files, samples should be in columns and features should be in rows.
2. Metadata table with matching sample IDs as the count table file. 

## Data input
Plotmicrobiome supports the common output files from sequencing analysis pipelines. The example files can be found in data-raw.

### Amplicon sequencing data:
1. DADA2 amplicon sequence variant (ASV) table with taxonomic classification in .csv, .tsv and .biom formats. 
   The .csv and .tsv tables have the ASVs in rows, and the taxonomic classification of ASVs in the last column. If you used 'biom convert' to convert .biom file to .tsv file, please remove the '#' from the first line.

2. A single file with all taxonomic levels and in .csv, .tsv and .biom formats.
   The taxonomic levels are in rows and should include from phylum to species level.
3. Multiple files with one for each taxonomic level separately.
   All the files should be stored in one directory and there shouldn't be other files in the folder.

### Shotgun metagenome sequencing data
The output of the most commonly used software MetaPhlAn and Kraken. The examples can be found in data-raw.

### One level data 
Data of one level, for example, the pathway abundance table from HUMAnN. Functions such as the taxonomic tree cannot be applied to this type of data. 

**App example:**
![alt text](https://github.com/ssun6/plotmicrobiome/blob/main/pics/1.data.png)

**R package Example:**
```
taxa_table=format_asv(taxa_file ="./data-raw/16S_biom_taxonomy.biom",biom=T,onefile = T,ASV=T)
```

```
taxa_table=format_wgs(taxa_file = "./data-raw/wgs_kraken2.txt",method="kraken",sep="\t")
taxa_table=format_wgs(taxa_file = "./data-raw/wgs_metaphlan2.txt",method="metaphlan",sep="\t")
```
```
taxa_table=format_pathway(taxa_file = "./data-raw/humann2_pathway.txt",sep="\t")
```

## Metadata input
Metadata should include sample IDs that are consistent with Data input files and variables for testing. Please specify the column of sample IDs. Shared samples between count table and metadata will be used for downstream analysis. 

**App example:**
![alt text](https://github.com/ssun6/plotmicrobiome/blob/main/pics/2.meta.png)

**R package Example:**
```
metadata=meta_format(metadata="./data-raw/metadata_cafe.csv",metadata_sep=",",meta_sample_name_col=2)
```


### MDS plot
Plotmicrobiome can generate both PCoA and non-parametric NMDS plots with selected dissimilarity measures.
PERMANOVA test results of the association between microbiome and the selected variable are shown as headers. 

**App example:**
![alt text](https://github.com/ssun6/plotmicrobiome/blob/main/pics/3.mds.png)

**R package Example:**
```
mds_plot(taxa_table = taxa_table, metadata=metadata,test_metadata="Study",method_mds = "pcoa",palette_group=c("red","blue","orange","green"),distance_type="bray")
mds_plot(taxa_table = taxa_table, metadata=metadata,test_metadata="Study",method_mds = "nmds",distance_type="bray")

```

### Alpha-diversity plot
Plotmicrobiome provides test results and visualization of the alpha-diversity difference between groups at each taxonomic level.
The metrics analyzed include Shannon index, Simpson index, Inverse Simpson index and numbers of taxa. 

**App example:**
![alt text](https://github.com/ssun6/plotmicrobiome/blob/main/pics/4.alpha.png)

**R package Example:**
```
alpha_plot(taxa_table = taxa_table, metadata=metadata,test_metadata="Study",method = "wilcoxon",palette_group=c("red","blue","orange","green"))
```

### Taxa barplot
Visualize the average taxonomic composition of each group from phylum to genus level.

**App example:**
![alt text](https://github.com/ssun6/plotmicrobiome/blob/main/pics/5.barplot.png)

**R package Example:**
```
taxa_barplot(taxa_table = taxa_table, metadata=metadata,test_metadata="Study",num_taxa=10,taxa_level="phylum",xlab_direction=1,legend_size=1)
```


### Data filter
The data can be stratified by metadata to only include samples that are used for statistical analysis. 
Abundance and Prevalence cutoffs can be used to exclude taxa of low abundance and/or low prevalence.

**App example:**
![alt text](https://github.com/ssun6/plotmicrobiome/blob/main/pics/6.datafilter.png)

**R package Example:**
```
tab_s=table_subset(taxa_table = taxa_table,metadata=metadata,stratify_by_metadata="Study",stratify_by_value="Sugar")
```


### Statistical tests
Statistical tests can be used to determine differential taxa between groups. 
Plotmicrobiome supports the following analyses.
Categorical data: t-test, Wilcoxon test, ANOVA and kruskal-wallis.
Continuous data: Pearson, Spearman and Kendall correlation.
Both categorical and continuous data: generalized linear models (glm), logistic regression (lr), mixed effects linear models (lme).

**Statistical tests use the output of data filter step.**

**App example:**
**Wilcoxon test:**
![alt text](https://github.com/ssun6/plotmicrobiome/blob/main/pics/7.wilcoxon.png)
**Generalized linear models (glm):**
![alt text](https://github.com/ssun6/plotmicrobiome/blob/main/pics/7.glm.png)
**Logistic regression (lr):**
![alt text](https://github.com/ssun6/plotmicrobiome/blob/main/pics/7.lr.png)
**Mixed effects linear models (lme):**
![alt text](https://github.com/ssun6/plotmicrobiome/blob/main/pics/7.lme.png)

**R package Example:**
```
fdrs1=stat_test(taxa_table =tab_s,metadata=metadata,test_metadata="Treatment",method="wilcoxon")
fdrs2=stat_test(taxa_table =tab_s,metadata=metadata,test_metadata="Treatment",method="lr",test_metadata_continuous=F,glm_anova=F,outcome_meta=T)
fdrs3=stat_test(taxa_table =tab_s,metadata=metadata,test_metadata="test_score",method="glm",test_metadata_continuous=T,glm_anova=F,outcome_meta=F,model_glm="factor(Treatment)")
fdrs4=stat_test(taxa_table =tab_s,metadata=metadata,test_metadata="test_score",method="lme",test_metadata_continuous=F,glm_anova=F,model_glm="factor(Treatment)",random_effect_var="RatID")
```

### Tree plots
Significant taxa are highlighted in a taxonomic tree.

**App example:**
![alt text](https://github.com/ssun6/plotmicrobiome/blob/main/pics/8.tree.png)

**R package Example:**
```
plot1=tree_view(taxa_table =tab_s, metadata=metadata,fdrs=fdrs1,test_metadata="Treatment",prevalence_cutoff=0.1, abundance_cutoff=0)
```


### Boxplots
Boxplots are used to visualize the differetial abundance of taxa identified in statistical tests.

**App example:**
![alt text](https://github.com/ssun6/plotmicrobiome/blob/main/pics/9.boxplot.png)

**R package Example:**
```
taxa_boxplot(taxa_table = tab_s, metadata=metadata,test_metadata="Study",fdrs=fdrs1,log_norm=T,cutoff=0.1,palette_group=c("red","blue","orange","green"))
```


### Correlation plots
The associations between taxa and continuous variables can be tested and visualized in this step.

**App example:**
![alt text](https://github.com/ssun6/plotmicrobiome/blob/main/pics/10.correlation_plot.png)


**R package Example:**
```
cor_plot1=meta_corplot(taxa_table =tab_s, metadata=metadata,test_metadata="test_score",col_metadata="group",fdr_cutoff=0.1)
```
