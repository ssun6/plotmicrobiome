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
library(shinyjs)
runApp('~/plotmicrobiome')
```
If you want to look at the code:
```
runApp('~/plotmicrobiome', display.mode = "showcase")
```



## 2. Website

Plotmicrobiome is freely available at https://ssun6.shinyapps.io/plotmicrobiome/ without installation. But the website is limited to 1GB memory and 25 active hours per month. 




## 3. R package

```
install.packages("devtools")
devtools::install_github("ssun6/plotmicrobiome")
```

# Load data
## Getting ready
1. Count table files in .csv, .tsv or .biom formats. For .csv and .tsv files, samples should be in columns and features should be in rows.
2. Metadata table with matching sample IDs as the count table file. 

## Data input
Plotmicrobiome supports the common output files from sequencing analysis pipelines. The example files can be found in data-raw.

### Amplicon sequencing data:
1. DADA2 amplicon sequence variant (ASV) table with taxonomic classification in .csv, .tsv and .biom formats. 
   The .csv and .tsv tables have the ASVs in rows, and the taxonomic classification of ASVs in the last column. If you used 'biom convert' to convert .biom file to .tsv file, please remove the '#' from the first line. Examples: 16S_biom_taxonomy.biom, table_taxonomy.txt.

2. A single file with all taxonomic levels and in .csv, .tsv and .biom formats.
   The taxonomic levels are in rows and should include from phylum to species level. Examples: taxa_all.txt, taxa_all.csv.
3. Multiple files with one for each taxonomic level separately.
   All the files should be stored in one directory and there shouldn't be other files in the folder. Examples: multiple_biom, multiple_tsv.

### Shotgun metagenome sequencing data
The output of the most commonly used software MetaPhlAn and Kraken. The examples can be found in data-raw: wgs_metaphlan2.txt, wgs_kraken2.txt

### One level data 
Data of one level, for example, the pathway abundance table from HUMAnN. Functions such as the taxonomic tree cannot be applied to this type of data.  Examples: humann2_pathway.txt.

### App example:

![alt text](https://github.com/ssun6/plotmicrobiome/blob/main/pics/1data.png)

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

### App example:
![alt text](https://github.com/ssun6/plotmicrobiome/blob/main/pics/2metadata.png)

**R package Example:**
```
metadata=meta_format(metadata="./data-raw/metadata_cafe.csv",metadata_sep=",",meta_sample_name_col=2)
```

# Analysis
## MDS plot
Plotmicrobiome can generate both PCoA and non-parametric NMDS plots with selected dissimilarity measures.
PERMANOVA test results of the association between microbiome and the selected variable are shown as headers. 

### App example:
![alt text](https://github.com/ssun6/plotmicrobiome/blob/main/pics/3mds.png)

**R package Example:**
```
mds_plot(taxa_table = taxa_table, metadata=metadata,test_metadata="Study",method_mds = "pcoa",palette_group=c("red","blue","orange","green"),distance_type="bray")
mds_plot(taxa_table = taxa_table, metadata=metadata,test_metadata="Study",method_mds = "nmds",distance_type="bray")

```

## Alpha-diversity plot
Plotmicrobiome provides test results and visualization of the alpha-diversity difference between groups at each taxonomic level.
The metrics analyzed include Shannon index, Simpson index, Inverse Simpson index and numbers of taxa. 

### App example:
![alt text](https://github.com/ssun6/plotmicrobiome/blob/main/pics/4alpha.png)

**R package Example:**
```
alpha_plot(taxa_table = taxa_table, metadata=metadata,test_metadata="Study",method = "wilcoxon",palette_group=c("red","blue","orange","green"))
```

## Taxa barplot
Visualize the average taxonomic composition of each group from phylum to genus level.

### App example:
![alt text](https://github.com/ssun6/plotmicrobiome/blob/main/pics/5barplot.png)

**R package Example:**
```
taxa_barplot(taxa_table = taxa_table, metadata=metadata,test_metadata="Study",num_taxa=10,taxa_level="phylum",xlab_direction=1,legend_size=1)
```


## Data filter
The data can be stratified by metadata to only include samples that are used for statistical analysis. 
Abundance and Prevalence cutoffs can be used to exclude taxa of low abundance and/or low prevalence.

### App example:
![alt text](https://github.com/ssun6/plotmicrobiome/blob/main/pics/6filter.png)

**R package Example:**
```
tab_s=table_subset(taxa_table = taxa_table,metadata=metadata,stratify_by_metadata="Study",stratify_by_value="Sugar")
```


## Statistical tests
**Plotmicrobiome supports the following analyses.** <br />
* Categorical data: t-test, Wilcoxon test, ANOVA and kruskal-wallis.<br />
* Continuous data: Pearson, Spearman and Kendall correlation.<br />
* Both categorical and continuous data: generalized linear models (glm), logistic regression (lr), mixed effects linear models (lme).<br />
<br />

Tips:
* The default family variable is binomial for categorical metadata and gaussian for continuous metadata. <br />
* For logistic regression (lr), the testing metadata should be a two level variable (e.g., case and control). And use the following parameters: Is the metadata outcome? True. Is the metadata for testing continuous? False.<br />
* Mixed effects linear models (lme) adjust for random effects. For example, when multiple samples were taken from one rat, the Rat ID can be included as random effect variable.<br />
* ANOVA can be run on linear models. It is helpful when the testing metadata is a multiple level variable.<br />
<br />

**Statistical tests use the output from the data filter step. Please run data filter first before running statistical tests**

### App example:
**Wilcoxon test:**
![alt text](https://github.com/ssun6/plotmicrobiome/blob/main/pics/7wilcoxon.png)
**Generalized linear models (glm):**
![alt text](https://github.com/ssun6/plotmicrobiome/blob/main/pics/7glm.png)
**Logistic regression (lr):**
![alt text](https://github.com/ssun6/plotmicrobiome/blob/main/pics/7lr.png)
**Mixed effects linear models (lme):**
![alt text](https://github.com/ssun6/plotmicrobiome/blob/main/pics/7lme.png)

**R package Example:**
```
fdrs1=stat_test(taxa_table =tab_s,metadata=metadata,test_metadata="Treatment",method="wilcoxon")
fdrs2=stat_test(taxa_table =tab_s,metadata=metadata,test_metadata="Treatment",method="lr",test_metadata_continuous=F,glm_anova=F,outcome_meta=T)
fdrs3=stat_test(taxa_table =tab_s,metadata=metadata,test_metadata="test_score",method="glm",test_metadata_continuous=T,glm_anova=F,outcome_meta=F,model_glm="factor(Treatment)")
fdrs4=stat_test(taxa_table =tab_s,metadata=metadata,test_metadata="test_score",method="lme",test_metadata_continuous=F,glm_anova=F,model_glm="factor(Treatment)",random_effect_var="RatID")
```

## Tree plots
Significant taxa are highlighted in a taxonomic tree. The taxonomic tree is pruned to remove some less prevalent and less abundant branches for visualization. These parameters can be adjusted in the tool to include less or more branches. 


### App example:
![alt text](https://github.com/ssun6/plotmicrobiome/blob/main/pics/8tree.png)

**R package Example:**
```
plot1=tree_view(taxa_table =tab_s, metadata=metadata,fdrs=fdrs1,test_metadata="Treatment",prevalence_cutoff=0.1, abundance_cutoff=0)
```


## Boxplots
Boxplots are used to visualize the differetial abundance of taxa identified in statistical tests. If there are not figures shown, please try adjusting the FDR cutoff. 

### App example:
![alt text](https://github.com/ssun6/plotmicrobiome/blob/main/pics/9box.png)

**R package Example:**
```
taxa_boxplot(taxa_table = tab_s, metadata=metadata,test_metadata="Study",fdrs=fdrs1,log_norm=T,cutoff=0.1,palette_group=c("red","blue","orange","green"))
```


## Correlation plots
The associations between taxa and continuous variables can be tested and visualized in this step. If there are not figures shown, please try adjusting the FDR cutoff. 

### App example:
![alt text](https://github.com/ssun6/plotmicrobiome/blob/main/pics/10cor.png)


**R package Example:**
```
cor_plot1=meta_corplot(taxa_table =tab_s, metadata=metadata,test_metadata="test_score",col_metadata="group",fdr_cutoff=0.1)
```

## P values vs P values plots
* P values vs P values plots can be used to compare the testing results from two groups. <br />
* The log10(P1) from dataset1 and log10(P2) from dataset2 are plotted against each other, and the correlations are calculated. <br />
* If there is a direction of changes, the log10(P) will be multiplied by +1/-1 to include the direction. ANOVA and Kruskal-wallis do not have direction of changes. <br />
* The P value files can be downloaded from the Statistical tests step or calculated at this step. <br />

### App example:
![alt text](https://github.com/ssun6/plotmicrobiome/blob/main/pics/11pvp.png)

**R package Example:**
```
p_compare(fdrs1,fdrs2,p_col1=2,p_col2=2,indicator1=4,indicator2=4,point_color="black",lab_cutoff=0.01,cor_method="spearman",x.reverse=T,exclude_unclassified=T,one_level=F,direction=T)
```
