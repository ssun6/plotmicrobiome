library(shiny)
library(shinyjs)
library(shinythemes)
#library(dplyr)
#library(stringr)
#library(rvest)
#library(httr)
#library(XML)


# Shiny UI -------
#increase upload size
options(shiny.maxRequestSize = 100*1024^2)
ui <- fluidPage(
  theme = shinytheme("sandstone"),
  tags$head(tags$style(
    HTML(".shiny-notification {
               position: fixed;
               top: 50%;
               left: 20%;
               padding: 20px 30px 40px 30px;
               text-align: center;
               font-weight: bold;
               font-size: 200%;
               color: darkred;
               background-color: white;
               z-index: 1000;
             }")
  )),
  navbarPage(
    title="PlotMicrobiome",
    id = "main_navbar",

    tabPanel(
      "Data input",
      sidebarLayout(
        sidebarPanel(
          useShinyjs(),
          width = 3,
          br(),
          h4("Data input"),
          selectInput("data_type", "Only select 16S or metagenomics for data with multiple taxonomic levels. See instructions for examples.", c("Amplicon with taxonomic structure","Metagenomics with taxonomic structure","Table without taxonomic structure")),
          br(),
          br(),
          div(id = "Amplicon with taxonomic structure",
            h4("Amplicon sequencing taxonomic data"),
            h5("Samples should be in column names and taxa should be in row names for .csv and .tsv files. If not, please transpose before uploading.", style="color:tomato"),
            fileInput("file_16s", "Choose taxa abundance file (accept .csv, .tsv and .biom files.)"),
            selectInput("biom", "Is the data in biom format?", c("True","False")),
            selectInput("ASV", "Is the data a DADA2 ASV table?", c("True","False")),
            selectInput("sep_16s", "What is the delimiter?", c(",","tab")),
            numericInput("n_reads_16s", "Exclude samples with the number of reads lower than", value = 1000),
            selectInput("norm_16s", "Normalize the data? (Supports proportion scaled by average sequencing depth and rarefaction)", c("True","False")),
            selectInput("rarefy_16s", "Use rarefaction for normalization? Default is proportion scaled by average sequencing depth.", c("False","True")),
            numericInput("rarefy_reads_16s", "Rarefy samples to how many reads?", value = 1000),
            numericInput("n_raw_16s", "Preview rows", value = 5, min = 1, step = 1),
            numericInput("n_raw_16s_col", "Preview columns", value = 10, min = 1, step = 1),
            br(),
            br(),
            br(),
            uiOutput("code_link_16s"),
            br()
          ),
          div(id = "Metagenomics with taxonomic structure",
            h4("Metagenomics taxonomic data"),
            h5("Samples should be in column names and taxa should be in row names for .csv and .tsv files. If not, please transpose before uploading.", style="color:tomato"),
            fileInput("file_wgs", "Choose metagenomics taxa file"),
            selectInput("method_wgs", "Which tool was used for taxonomic classification?", c("kraken","metaphlan")),
            selectInput("sep_wgs", "What is the delimiter?", c("tab",",")),
            numericInput("n_reads_wgs", "Exclude samples with the number of reads lower than", value = 1000),
            selectInput("norm_wgs", "Normalize the data? (Supports proportion scaled by average sequencing depth and rarefaction)", c("True","False")),
            selectInput("rarefy_wgs", "Use rarefaction for normalization? Default is proportion scaled by average sequencing depth.", c("False","True")),
            numericInput("rarefy_reads_wgs", "Rarefy samples to how many reads?", value = 1000),
            numericInput("n_raw_wgs", "Preview rows", value = 5, min = 1, step = 1),
            numericInput("n_raw_wgs_col", "Preview columns", value = 10, min = 1, step = 1),
            br(),
            br(),
            br(),
            uiOutput("code_link_wgs"),
            br()
          ),
          div(id = "Table without taxonomic structure",
            h4("Table without taxonomic structure"),
            h5("Samples should be in column names and taxa should be in row names for .csv and .tsv files. If not, please use the transpose option below.", style="color:tomato"),
            fileInput("file_tab", "Choose file"),
            selectInput("sep_tab", "What is the delimiter?", c(",","tab")),
            selectInput("transpose_tab", "Transpose the table? Samples should be in column names and taxa should be in row names.", c("False","True")),
            numericInput("n_reads_tab", "Exclude samples with the number of reads lower than", value = 0),
            selectInput("norm_tab", "Normalize the data? (Supports proportion scaled by average sequencing depth and rarefaction)", c("False","True")),
            selectInput("rarefy_tab", "Use rarefaction for normalization? Default is proportion scaled by average sequencing depth.", c("False","True")),
            numericInput("rarefy_reads_tab", "Rarefy samples to how many reads?", value = 1000),
            numericInput("n_raw_tab", "Preview rows", value = 5, min = 1, step = 1),
            numericInput("n_raw_tab_col", "Preview columns", value = 10, min = 1, step = 1),
            br(),
            br(),
            br(),
            uiOutput("code_link_tab"),
            br()
          ),
          actionButton("button_raw", "Run"),
          textOutput("data_dim"),
          br(),
          br(),
          h5("Download code:"),
          downloadButton("code_datainput_Download", label ="Download"),
          br(),
          br(),
          br(),
          h5("Download table:"),
          downloadButton("downloadInput", "Download")
        ),
        mainPanel(tableOutput("head_raw"))
      )
    ),


    tabPanel(
      "Metadata input",
      sidebarLayout(
        sidebarPanel(
          width = 3,
          br(),
          h4("Metadata input"),
          h5("Samples should be in row names and metadata should be in column names."),
          fileInput("file_meta", "Choose metadata file"),
          #textInput("dir_meta", "Metadata directory",value="./data-raw/metadata_cafe.csv"),
          selectInput("sep_meta", "What is the delimiter?", c(",","tab")),
          numericInput("meta_sample_name_col", "Which column are the sample names in?", value = 0,step = 1),
          numericInput("n_meta", "Preview rows", value = 5, min = 1, step = 1),
          actionButton("button_meta", "Run"),
          br(),
          br(),
          h5("Download code:"),
          downloadButton("code_metadata_Download", label ="Download"),
          br(),
          br(),
          br(),
          uiOutput("code_link_meta"),
          br()
        ),
        mainPanel(tableOutput("head_meta"))
      )
    ),

    tabPanel(
      "MDS plot",
      sidebarLayout(
        sidebarPanel(
          width = 3,
          br(),
          h4("MDS plot"),
          selectInput("stratify_by_metadata_mds", "If a subset of the data will be used, select the variable for subsetting (select no for no subsetting)",c("")),
          h5("Variable preview:"),
          textOutput("head_stratify_mds"),
          br(),
          selectInput("stratify_by_value_mds", "Select the values for subsetting (keep blank for no subsetting)",multiple = TRUE,""),
          #checkboxGroupInput("stratify_by_value_mds", "Select the values for subsetting",""),
          #textInput("stratify_by_value_mds", "Select the values for stratification (separated with comma)",value=NULL),
          selectInput("test_metadata_mds", "Select the metadata for testing",c("")),
          h5("Variable preview:"),
          textOutput("head_test_mds"),
          br(),
          uiOutput(outputId = 'taxa_level_mdsUI'),
          selectInput("method_mds", "Which method should be used for ordination?", c("pcoa","nmds")),
          selectInput("distance_type", "Distance type",c("bray","euclidean","manhattan","jaccard")),
          selectInput("log_normalization_mds", "Should the data be log10 normalized?", c("False","True")),
          textAreaInput("palette_group_mds", "Colors for plot", value = "#3498DB,#E74C3C,#2ECC71,#9B59B6,#F39C12,#1ABC9C,#D35400,hotpink,darkgrey"),
          sliderInput("dot_size_mds", label ="Size of the points", min = 0.5, max = 5, value = 1.5),
          sliderInput("dot_transparency_mds", label ="Transparency of the points", min = 0, max = 1, value = 0.6),
          sliderInput("label_size_mds", label ="Size of figure labels", min = 0.5, max = 5, value = 1.3),
          sliderInput("ellipse_label_size_mds", label ="Size of ellipse labels", min = 0.5, max = 5, value = 1.3),
          sliderInput("legend_label_size_mds", label ="Size of legend labels", min = 0.5, max = 5, value = 1.3),
          selectInput("show_sample_name_mds", "Show sample names?", c("False","True")),
          actionButton("button_mds", "Run"),
          br(),
          br(),
          br(),
          sliderInput("mds_h", label ="Figure height", min = 0.5, max = 5, value = 1),
          sliderInput("mds_w", label ="Figure width", min = 0.5, max = 5, value = 1),
          br(),
          h5("Download figure:"),
          downloadButton("plotMDSDownload", label ="Download"),
          br(),
          br(),
          h5("Download code:"),
          downloadButton("codeMDSDownload", label ="Download"),
          # br(),
          # htmlOutput("function_mds"),# to use new line \n
          # br(),
          # uiOutput("code_link_mds"),
          br(),
          br()
        ),
        mainPanel(uiOutput("plotMDS.ui"))
      )
    ),

    tabPanel(
      "Alpha diversity plot",
      sidebarLayout(
        sidebarPanel(
          width = 3,
          br(),
          h4("Alpha diversity plot"),
          selectInput("stratify_by_metadata_alpha", "Select the variable for stratification",c("")),
          h5("Variable preview:"),
          textOutput("head_stratify_alpha"),
          br(),
          selectInput("stratify_by_value_alpha", "Select the values for stratification","",multiple = TRUE),
          #textInput("stratify_by_value_alpha", "Select the values for stratification (separated with comma)",value=NULL),
          selectInput("test_metadata_alpha", "Select the metadata for testing",c("")),
          h5("Variable preview:"),
          textOutput("head_test_alpha"),
          br(),
          checkboxInput("test_metadata_order_alpha_checkbox","Use default order of metadata", TRUE),
          selectInput("test_metadata_order_alpha", "Select the order of metadata to change those in the figure ","",multiple = TRUE),
          textInput("xlab_alpha", "X axis label",value="default"),
          selectInput("method_alpha", "Select statistical test method",c("wilcoxon","t.test","kruskal-wallis","anova")),
          textAreaInput("palette_group_alpha", "Colors for plot", value = "#3498DB,#E74C3C,#2ECC71,#9B59B6,#F39C12,#1ABC9C,#D35400"),
          selectInput("x_dir_alpha", "Direction of X axis labels (1 is horizontal, 2 is vertical)", c(1,2)),
          actionButton("button_alpha", "Run"),
          br(),
          br(),
          sliderInput("alpha_h", label ="Figure height", min = 0.5, max = 5, value = 1.65),
          sliderInput("alpha_w", label ="Figure width", min = 0.5, max = 5, value = 1.5),
          br(),
          br(),
          h5("Download figure:"),
          downloadButton("plotAlphaDownload", "Figure Download"),
          br(),
          h5("Download alpha diversity table:"),
          downloadButton("downloadAlpha", "Table Download"),
          br(),
          br(),
          br(),
          br(),
          h5("Download code:"),
          downloadButton("code_alpha_Download", label ="Download"),
          br(),
          br(),
          #htmlOutput("function_alpha"),# to use new line \n
          br(),
          uiOutput("code_link_alpha"),
          br(),
          br()
        ),
        mainPanel(uiOutput("plotAlpha.ui"))
      )
    ),

    tabPanel(
      "Taxa barplot",
      sidebarLayout(
        sidebarPanel(
          width = 3,
          br(),
          h4("Taxa barplot"),
          selectInput("stratify_by_metadata_bar", "Select the variable for stratification",c("")),
          h5("Variable preview:"),
          textOutput("head_stratify_bar"),
          br(),
          selectInput("stratify_by_value_bar", "Select the values for stratification","",multiple = TRUE),
          #textInput("stratify_by_value_bar", "Select the values for stratification (separated with comma)",value=NULL),
          selectInput("test_metadata_bar", "Select the metadata for testing",c("")),
          h5("Variable preview:"),
          textOutput("head_test_bar"),
          br(),
          checkboxInput("test_metadata_order_bar_checkbox","Use default order of metadata", TRUE),
          selectInput("test_metadata_order_bar", "Select the order of metadata to change those in the figure ","",multiple = TRUE),
          #textInput("test_metadata_order_bar", "Type in the order of metadata separated with comma to change those in the figure ",value="default"),
          textInput("num_taxa_bar", "Select the number of taxa shown",value=8),
          selectInput("taxa_level_bar", "Select the taxonomic level shown",c("phylum","class","order","family","genus")),
          textAreaInput("palette_group_bar", "Colors for plot", value = "default"),
          selectInput("x_dir_bar", "Direction of X axis labels (1 is horizontal, 2 is vertical)",c(1,2)),
          textInput("legend_size_bar", "Select the legend size", value = 1.5),
          actionButton("button_bar", "Run"),
          br(),
          br(),
          sliderInput("bar_h", label ="Figure height", min = 0.5, max = 5, value = 1),
          sliderInput("bar_w", label ="Figure width", min = 0.5, max = 10, value = 1),
          br(),
          br(),
          h5("Download figure:"),
          downloadButton("plotBarDownload", "Download"),
          br(),
          h5("Download barplot composition table:"),
          downloadButton("downloadBar", "Table Download"),
          br(),
          br(),
          br(),
          h5("Download code:"),
          downloadButton("code_bar_Download", label ="Download"),
          br(),
          br(),
          br(),
          uiOutput("code_link_bar"),
          br(),
          br()
        ),
        mainPanel(uiOutput("plotBar.ui"))
      )
    ),

    tabPanel(
      "Data filter",
      sidebarLayout(
        sidebarPanel(
          width = 3,
          br(),
          h4("Data filter"),
          selectInput("stratify_by_metadata_filter", "Select the variable for stratification",c("")),
          h5("Variable preview:"),
          textOutput("head_stratify_metadata_filter"),
          br(),
          selectInput("stratify_by_value_filter", "Select the values for stratification","",multiple = TRUE),
          selectInput("exclude_ASV_filter", "Should ASV or strain level be excluded from the analysis (this changes the P-value distribution and FDR)?", c("True","False")),
          numericInput("prevalence_cutoff", "Prevalence cutoff", value = 0.25),
          numericInput("abundance_cutoff", "Abundance cutoff", value = 0),
          numericInput("n_filtered", "Preview rows", value = 5, min = 1, step = 1),
          numericInput("n_filtered_col", "Preview columns", value = 10, min = 1, step = 1),
          actionButton("button_filter", "Run"),
          textOutput("filtered_dim"),
          br(),
          br(),
          br(),
          h5("Download code:"),
          downloadButton("code_filter_Download", label ="Download"),
          br(),
          br(),
          br(),
          uiOutput("code_link_subset"),
          br(),
          br()
        ),
        mainPanel(tableOutput("head_filtered"))
      )
    ),

    tabPanel(
      "Statistical test",
      sidebarLayout(
        sidebarPanel(
          useShinyjs(),
          width = 3,
          br(),
          h4("Statistical test"),
          h5("Statistical test uses the output of data filter. Please run data filter first."),
          selectInput(inputId = "test_metadata_stat", "Select the metadata for testing",c("")),
          h5("Variable preview:"),
          textOutput("head_test_stat"),
          selectInput("test_metadata_continuous_stat", "Is the metadata for testing continuous?",c("False","True"),selected ="False"),
          uiOutput(outputId = 'method_statUI'),
          selectInput(inputId = "log_norm_stat", "Should the data be log10 normalization?",c("True","False"),selected="True"),
          #selectInput("outcome_meta", "Is the metadata used as the outcome (Should be True for lr test, can be True or False for glm and lme ?",c("False","True"),selected="False"),
          uiOutput(outputId = 'outcome_metaUI'),
          uiOutput(outputId = 'model_glmUI'),
          uiOutput(outputId = 'glm_anovaUI'),
          uiOutput(outputId = 'glm_distUI'),
          uiOutput(outputId = 'random_effect_varUI'),
          #selectInput(inputId = "glm_anova", "Run ANOVA on linear models?",c("False","True"),selected="False"),
          #textAreaInput(inputId = "model_glm", "Covariates for adjusting (e.g., age+factor(sex)+factor(batch))", value = ""),
          #textAreaInput(inputId = "glm_dist", "Family variable of glm function (e.g.,binomial, gaussian, poisson). Default is gaussian for continuous outcomes, and binomial for two-category outcomes.", value = "default"),
          #selectInput(inputId = "random_effect_var", "Random effect variable for mixed effects models",c("")),
          br(),
          selectInput("sort_fdr", "Sort by FDR",c("True","False")),
          numericInput("n_fdrs", "Preview rows", value = 5, min = 1, step = 1),
          actionButton("button_fdrs", "Run"),
          br(),
          br(),
          h5("Download statistical test results:"),
          downloadButton("downloadStat", "Download"),
          br(),
          br(),
          br(),
          h5("Download code:"),
          downloadButton("code_stat_Download", label ="Download"),
          br(),
          br(),
          br(),
          uiOutput("code_link_stat"),
          br()
        ),
        mainPanel(
          tabsetPanel(
            tabPanel("Results", tableOutput("head_fdrs")),
            tabPanel("Histogram of P values", uiOutput("plotPvalsHist.ui"))
          )
        )
      )
    ),

    tabPanel(
      "Tree plot",
      sidebarLayout(
        sidebarPanel(
          width = 3,
          br(),
          h4("Tree plot"),
          h5("Parameters for pruning tree:"),
          numericInput("prevalence_cutoff_tree", "Prevalence cutoff", value = 0.1,min=0),
          numericInput("abundance_cutoff_tree", "Abundance cutoff", value = 0,min=0),
          selectInput("branch_tree", "Choose branch for the tree.",c("Bacteria","Archaea","Fungi")),
          h5("Please note that Archaea and Fungi trees are only recommended when they have high diversity and abundance. When their abundance are much lower than bacteria, the differential abundance results can be largely impacted by data compositionality."),
          br(),
          br(),
          selectInput("test_metadata_continuous_tree", "Is the test metadata continuous?",c("False","True")),
          numericInput("fdr_cutoff_tree", "FDR cutoff", value = 0.1, min = 0, step = 0.01),
          textAreaInput("node_size_breaks", "Breaks for node size", value = "0,0.01,0.05,0.5,5"),
          textAreaInput("palette_group_tree", "Colors for plot", value = "#3498DB,#E74C3C,#2ECC71,#9B59B6,#F39C12,#1ABC9C,#D35400"),
          textInput("taxa_removal_tree", "Remove taxa from plot?",value=NULL),
          selectInput("single_parent_branch_removal","Remove the parent branch if there is only one child branch within the parent branch?",c("False","True")),
          selectInput("single_child_branch_removal","Remove the child branch if thgvere is only one child branch within the parent branch?",c("False","True")),
          actionButton("button_tree", "Run"),
          br(),
          br(),
          sliderInput("tree_h", label ="Figure height", min = 0.5, max = 5, value = 1),
          sliderInput("tree_w", label ="Figure width", min = 0.5, max = 5, value = 1),
          br(),
          br(),
          h5("Download figure:"),
          downloadButton("plotTreeDownload", "Download"),
          br(),
          br(),
          h5("Download code:"),
          downloadButton("code_tree_Download", label ="Download"),
          br(),
          br(),
          uiOutput("code_link_tree"),
          br(),
          br()
        ),
        mainPanel(uiOutput("plotTree.ui"))
      )
    ),

    tabPanel(
      "Boxplot plot",
      sidebarLayout(
        sidebarPanel(
          width = 3,
          br(),
          h4("Boxplot plot"),
          selectInput("log_norm_box", "Log normalization?",c("True","False")),
          numericInput("fdr_cutoff_box", "FDR cutoff", value = 0.1, min = 0, step = 0.01),
          checkboxInput("test_metadata_order_box_checkbox","Use default order of metadata", TRUE),
          selectInput("test_metadata_order_box", "Select the order of metadata to change those in the figure","",multiple = TRUE),
          textInput("xlab_box", "X axis label",value="default"),
          textInput("ylab_box", "Y axis label",value="default"),
          #textInput("test_metadata_order_box", "Type in the order of metadata separated with comma to change those in the figure ",value="default"),
          textInput("taxa_shown_box", "Select specific taxa", value = ""),
          textAreaInput("palette_group_box", "Colors for plot", value = "red,blue,orange,green"),
          selectInput("x_dir_box", "Direction of X axis labels (1 is horizontal, 2 is vertical)", c(1,2)),
          br(),
          br(),
          numericInput("page_box", "Page number", value = 1, min = 1, step = 1),
          actionButton("button_box", "Run"),
          br(),
          br(),
          sliderInput("box_h", label ="Figure height", min = 0.5, max = 5, value = 1),
          sliderInput("box_w", label ="Figure width", min = 0.5, max = 5, value = 1),
          br(),
          br(),
          h5("Download figure:"),
          downloadButton("plotBoxDownload", "Download"),
          br(),
          br(),
          br(),
          h5("Download code:"),
          downloadButton("code_box_Download", label ="Download"),
          br(),
          br(),
          br(),
          uiOutput("code_link_box"),
          uiOutput("code_link_box_download"),
          br(),
          br(),
          br()
        ),
        mainPanel(
          mainPanel(uiOutput("plotBox.ui"))
        )
      )
    ),

    tabPanel(
      "Correlation plot",
      sidebarLayout(
        sidebarPanel(
          width = 3,
          br(),
          h4("Correlation plot"),
          selectInput("test_metadata_cor", "Select the metadata for testing",c("")),
          h5("Variable preview:"),
          textOutput("head_metadata_cor_test"),
          br(),
          selectInput("col_metadata_cor", "Select the metadata for color",c("")),
          h5("Variable preview:"),
          textOutput("head_metadata_cor_col"),
          br(),
          textInput("ylab_cor", "Y axis label",value="default"),
          selectInput("cor_method", "Correlation method",c("spearman","pearson","kendall")),
          selectInput("log_norm_cor", "Log normalization?",c("True","False")),
          numericInput("fdr_cutoff_cor", "FDR cutoff \n(Try increasing the cutoff if there is no taxa shown)", value = 0.1, min = 0, step = 0.01),
          textInput("taxa_shown_cor", "Select taxa", value = ""),
          textAreaInput("palette_group_cor", "Colors for plot", value = "red,blue,orange,green"),
          br(),
          br(),
          numericInput("page_cor", "Page number", value = 1, min = 1, step = 1),
          actionButton("button_cor", "Run"),
          br(),
          br(),
          sliderInput("cor_h", label ="Figure height", min = 0.5, max = 5, value = 1),
          sliderInput("cor_w", label ="Figure width", min = 0.5, max = 5, value = 1),
          br(),
          br(),
          h5("Download figure:"),
          downloadButton("plotCorDownload", "Download"),
          br(),
          br(),
          br(),
          h5("Download code:"),
          downloadButton("code_cor_Download", label ="Download"),
          br(),
          br(),
          br(),
          uiOutput("code_link_cor"),
          uiOutput("code_link_cor_download"),
          br(),
          br(),
          br()
        ),
        mainPanel(
          mainPanel(uiOutput("plotCor.ui"))
        )
      )
    ),

    tabPanel(
      "P vs P plot",
      sidebarLayout(
        sidebarPanel(
          width = 3,
          br(),
          h4("Pvalue vs Pvalue plot"),
          br(),
          h4("Data1:"),
          selectInput("upload_p1", "Upload existing P value files for data1?",c("Yes","No")),
          br(),
          br(),
          h5("Upload existing P value files for data1"),
          br(),
          fileInput("file_p1", "Choose P value file (accept .csv and .tsv.)"),
          selectInput("sep_p1", "What is the delimiter?", c(",","tab")),
          br(),
          br(),
          h5("Generate P values for data1"),
          selectInput("stratify_by_metadata_p1", "Select the variable for stratification",c("")),
          h5("Variable preview:"),
          textOutput("head_stratify_p1"),
          br(),
          selectInput("stratify_by_value_p1", "Select the values for stratification","",multiple = TRUE),
          #textInput("stratify_by_value_p1", "Select the values for stratification (separated with comma)",value=NULL),
          selectInput("exclude_ASV_p1", "Should ASV or strain be excluded from the analysis (this changes the P-value distribution and FDR)?", c("True","False")),
          numericInput("prevalence_cutoff_p1", "Prevalence cutoff", value = 0.25),
          numericInput("abundance_cutoff_p1", "Abundance cutoff", value = 0),
          selectInput("test_metadata_p1", "Select the metadata for testing",c("")),
          selectInput("test_metadata_continuous_p1", "Is the metadata for testing continuous?",c("False","True")),
          h5("Variable preview:"),
          textOutput("head_test_p1"),
          selectInput("method_stat_p1", "Select method",c("wilcoxon","t.test","kruskal-wallis","anova","pearson","spearman","kendall")),
          br(),
          br(),
          actionButton("button_p1_2", "Run"),
          br(),
          br(),
          br(),
          br(),
          br(),
          h4("Data2:"),
          selectInput("upload_p2", "Upload existing P value files for data2?",c("Yes","No")),
          br(),
          br(),
          h5("Upload existing P value files for data2"),
          br(),
          fileInput("file_p2", "Choose P value file (accept .csv and .tsv.)"),
          selectInput("sep_p2", "What is the delimiter?", c(",","tab")),
          br(),
          br(),
          h5("Generating P value files for data2"),
          selectInput("stratify_by_metadata_p2", "Select the variable for stratification",c("")),
          h5("Variable preview:"),
          textOutput("head_stratify_p2"),
          br(),
          selectInput("stratify_by_value_p2", "Select the values for stratification","",multiple = TRUE),
          #textInput("stratify_by_value_p2", "Select the values for stratification (separated with comma)",value=NULL),
          selectInput("exclude_ASV_p2", "Should ASV or strain be excluded from the analysis (this changes the P-value distribution and FDR)?", c("True","False")),
          numericInput("prevalence_cutoff_p2", "Prevalence cutoff", value = 0.25),
          numericInput("abundance_cutoff_p2", "Abundance cutoff", value = 0),
          selectInput("test_metadata_p2", "Select the metadata for testing",c("")),
          selectInput("test_metadata_continuous_p2", "Is the metadata for testing continuous?",c("False","True")),
          h5("Variable preview:"),
          textOutput("head_test_p2"),
          selectInput("method_stat_p2", "Select method",c("wilcoxon","t.test","kruskal-wallis","anova","pearson","spearman","kendall")),
          br(),
          br(),
          actionButton("button_p2_2", "Run"),
          br(),
          br(),
          br(),
          br(),
          br(),
          h4("Correlation plot parameters:"),
          selectInput("direction_pvp", "Is there a direction of changes? Please select False for ANOVA and Kruskal-wallis tests.",c("True","False")),
          numericInput("p1_col", "Select the column showing P values in data 1 file", value = 2, min = 1, step = 1),
          numericInput("ind1_col", "Select the column showing the group indicators in data 1 file", value = 4, min = 1, step = 1),
          numericInput("p2_col", "Select the column showing P values in data 2 file", value = 2, min = 1, step = 1),
          numericInput("ind2_col", "Select the column showing the group indicators in data 2 file", value = 4, min = 1, step = 1),
          textAreaInput("point_color", "Colors for points", value = "red"),
          numericInput("lab_cutoff", "P value cutoff for labeling points", value = 0.05, min = 0, step = 0.01),
          selectInput("cor_method_p", "Select the correlation methods comparing P values",c("pearson","spearman","kendall")),
          selectInput("x_reverse", "Reverse the x axis?", c("False","True")),
          selectInput("y_reverse", "Reverse the y axis?", c("False","True")),
          selectInput("exclude_unclassified", "Exclude labels of unclassified taxa from plot?", c("True","False")),
          br(),
          br(),
          actionButton("button_cor_p", "Run"),
          br(),
          br(),
          sliderInput("pvp_h", label ="Figure height", min = 0.5, max = 5, value = 1),
          sliderInput("pvp_w", label ="Figure width", min = 0.5, max = 5, value = 1),
          br(),
          br(),
          br(),
          br(),
          br(),
          h5("Download figure:"),
          downloadButton("plotPvalsDownload", "Download"),
          br(),
          br(),
          br(),
          h5("Download code:"),
          downloadButton("code_pvp_Download", label ="Download"),
          br(),
          br(),
          br(),
          uiOutput("code_link_pvp"),
          br(),
          br()
        ),
        mainPanel(
          tabsetPanel(
            tabPanel("Data1", tableOutput("head_p1")),
            tabPanel("Data2", tableOutput("head_p2")),
            tabPanel("Correlation plot", uiOutput("plotPvals.ui"))
          )
        )
      )
    ),
    tags$head(tags$style(type="text/css", "
             #loadmessage {
               position: fixed;
               top: 50%;
               left: 50%;
               padding: 20px 30px 40px 30px;
               text-align: center;
               font-weight: bold;
               font-size: 200%;
               color: darkgrey;
               background-color: white;
               z-index: 1000;
             }
          ")),
    conditionalPanel(condition="$('html').hasClass('shiny-busy')",
                     tags$div("Loading...",id="loadmessage"))

  )
)

server <- function(input, output, session) {
  #taxonomic table input
  observeEvent(input$data_type, {

    if(input$data_type == "Amplicon with taxonomic structure"){
      shinyjs::show(id = "Amplicon with taxonomic structure")
      shinyjs::hide(id = "Metagenomics with taxonomic structure")
      shinyjs::hide(id = "Table without taxonomic structure")
    }else if(input$data_type == "Metagenomics with taxonomic structure"){
      shinyjs::show(id = "Metagenomics with taxonomic structure")
      shinyjs::hide(id = "Amplicon with taxonomic structure")
      shinyjs::hide(id = "Table without taxonomic structure")
    }else{
      shinyjs::show(id = "Table without taxonomic structure")
      shinyjs::hide(id = "Amplicon with taxonomic structure")
      shinyjs::hide(id = "Metagenomics with taxonomic structure")
    }
  })

  data_raw <- eventReactive(input$button_raw,{
    if(input$data_type=="Amplicon with taxonomic structure"){
      if(input$sep_16s=="tab"){
        sep_char_16s="\t"
      }else{
        sep_char_16s=input$sep_16s
      }
      try(format_asv(taxa_file =input$file_16s$datapath,biom=as.logical(input$biom),ASV=as.logical(input$ASV),sep=sep_char_16s,reads_cutoff=as.numeric(input$n_reads_16s),normalization=as.logical(input$norm_16s),rarefy=as.logical(input$rarefy_16s),rarefy_num=as.numeric(input$rarefy_reads_16s)))
    }else if (input$data_type=="Metagenomics with taxonomic structure"){
      if(input$sep_wgs=="tab"){
        sep_char_wgs="\t"
      }else{
        sep_char_wgs=input$sep_wgs
      }
      try(format_wgs(taxa_file =input$file_wgs$datapath,sep=sep_char_wgs,method=input$method_wgs,reads_cutoff=as.numeric(input$n_reads_wgs),normalization=as.logical(input$norm_wgs),rarefy=as.logical(input$rarefy_wgs),rarefy_num=as.numeric(input$rarefy_reads_wgs)))
    }else if (input$data_type=="Table without taxonomic structure"){
      if(input$sep_tab=="tab"){
        sep_char_tab="\t"
      }else{
        sep_char_tab=input$sep_tab
      }
      try(format_tabs(taxa_file =input$file_tab$datapath,sep=sep_char_tab,transpose=input$transpose_tab,reads_cutoff=as.numeric(input$n_reads_tab),normalization=as.logical(input$norm_tab),rarefy=as.logical(input$rarefy_tab),rarefy_num=as.numeric(input$rarefy_reads_tab)))
   }})

  #notifications when sample type is not correct
  observeEvent(input$button_raw, {
    # Show a modal when the button is pressed
    if(!is.null(class(data_raw()))){
      if(class(data_raw())[1]=="try-error"){
        showNotification(paste("Please double check that you selected the correct data type. Consider use 'table without taxonomic structure' as data type or check formats at https://github.com/ssun6/plotmicrobiome"),duration = 10, type = "error")
      }
    }
  })


  output$head_raw <- renderTable({
    if(input$data_type=="Amplicon with taxonomic structure"){
      if(ncol(data_raw())<as.numeric(input$n_raw_16s_col)){
        n_col=ncol(data_raw())
      }else{
        n_col=as.numeric(input$n_raw_16s_col)
      }
      head(data_raw()[,1:n_col], input$n_raw_16s)
    }else if (input$data_type=="Metagenomics with taxonomic structure"){
      if(ncol(data_raw())<as.numeric(input$n_raw_wgs_col)){
        n_col=ncol(data_raw())
      }else{
        n_col=as.numeric(input$n_raw_wgs_col)
      }
      head(data_raw()[,1:n_col], input$n_raw_wgs)
      }else if (input$data_type=="Table without taxonomic structure"){
        if(ncol(data_raw())<as.numeric(input$n_raw_tab_col)){
          n_col=ncol(data_raw())
        }else{
          n_col=as.numeric(input$n_raw_tab_col)
        }
        head(data_raw()[,1:n_col], input$n_raw_tab)
      }},rownames = TRUE)

  output$downloadInput <- downloadHandler(
    filename = function() {
      "input.csv"
    },
    content = function(file) {
      write.csv(data_raw(), file, row.names = TRUE)
    }
  )

  one_level_all = eventReactive(input$data_type,{
    if(input$data_type=="Amplicon with taxonomic structure" | input$data_type=="Metagenomics with taxonomic structure"){
      "FALSE"
    }else{
      "TRUE"
    }
  })

  # output data dimensions
  output$data_dim <- renderText({
    paste("Number of rows :",nrow(data_raw()),"\nNumber of columns :",ncol(data_raw()))
  })




  #format_asv(taxa_file =input$file_16s$datapath,biom=as.logical(input$biom),ASV=as.logical(input$ASV),sep=sep_char_16s,reads_cutoff=as.numeric(input$n_reads_16s),normalization=as.logical(input$norm_16s),rarefy=as.logical(input$rarefy_16s),rarefy_num=as.numeric(input$rarefy_reads_16s))
  # output$function_16s <- renderText({
  #   if(input$sep_16s=="tab"){
  #     sep_char_16s="\\t"
  #   }else{
  #     sep_char_16s=input$sep_16s
  #   }
  #   paste0("tab1 = format_asv(taxa_file = \"",input$file_16s$name,"\", biom = ",input$biom,", ASV = ",input$ASV,", sep = \"",sep_char_16s,"\", reads_cutoff = ",input$n_reads_16s,", normalization = ",input$norm_16s,", rarefy = ",input$rarefy_16s,", rarefy_num = ",input$rarefy_reads_16s,")")
  # })
  #
  # output$function_wgs <- renderText({
  #   if(input$sep_wgs=="tab"){
  #     sep_char_wgs="\\t"
  #   }else{
  #     sep_char_wgs=input$sep_wgs
  #   }
  #   #format_wgs(taxa_file =input$file_wgs$datapath,sep=sep_char_wgs,method=input$method_wgs,reads_cutoff=as.numeric(input$n_reads_wgs),normalization=as.logical(input$norm_wgs),rarefy=as.logical(input$rarefy_wgs),rarefy_num=as.numeric(input$rarefy_reads_wgs))
  #   paste0("tab1 = format_wgs(taxa_file = \"",input$file_wgs$name,"\", sep = \"",sep_char_wgs,"\", reads_cutoff = ",input$n_reads_wgs,", method = \"",input$method_wgs,"\", normalization = ",input$norm_wgs,", rarefy = ",input$rarefy_wgs,", rarefy_num = ",input$rarefy_reads_wgs,")")
  # })
  #
  # output$function_tab <- renderText({
  #   if(input$sep_tab=="tab"){
  #     sep_char_tab="\\t"
  #   }else{
  #     sep_char_tab=input$sep_tab
  #   }
  #   #format_tabs(taxa_file =input$file_tab$datapath,sep=sep_char_tab,reads_cutoff=as.numeric(input$n_reads_tab),normalization=as.logical(input$norm_tab),rarefy=as.logical(input$rarefy_tab),rarefy_num=as.numeric(input$rarefy_reads_tab))
  #   paste0("tab1 = format_tabs(taxa_file = \"",input$file_tab$name,"\", sep = \"",sep_char_tab,"\", reads_cutoff = ",input$n_reads_tab,", normalization = ",input$norm_tab,", rarefy = ",input$rarefy_tab,", rarefy_num = ",input$rarefy_reads_tab,")")
  # })
  #
  # output$function_meta <- renderText({
  #   if(input$sep_meta=="tab"){
  #     sep_char_meta="\\t"
  #   }else{
  #     sep_char_meta=input$sep_meta
  #   }
  #   #meta_format(metadata=input$file_meta$datapath,metadata_sep=sep_char_meta,meta_sample_name_col=input$meta_sample_name_col)
  #   paste0("metadata1 = meta_format(metadata = \"",input$file_meta$name,"\", metadata_sep = \"",sep_char_meta,"\", meta_sample_name_col = ",input$meta_sample_name_col,")")
  # })

  #MDS show code run
  # output$function_mds <- renderUI({
  #   #dataf1=table_subset(taxa_table = data_raw(),metadata=data_meta(),stratify_by_metadata=input$stratify_by_metadata_mds,stratify_by_value=input$stratify_by_value_mds,one_level=as.logical(one_level_all()))
  #   #mds_plot(taxa_table =dataf1,metadata=data_meta(),test_metadata=input$test_metadata_mds,taxa_level=input$taxa_level_mds,method_mds=input$method_mds,one_level=as.logical(one_level_all()),log_norm=as.logical(input$log_normalization_mds),palette_group=strsplit(input$palette_group_mds, ",\\s*")[[1]],distance_type=input$distance_type,dot_transparency=as.numeric(input$dot_transparency_mds),dot_size=as.numeric(input$dot_size_mds),show_sample_name=as.logical(input$show_sample_name_mds),ellipse_label_size=as.numeric(input$ellipse_label_size_mds),label_size=as.numeric(input$label_size_mds),legend_label_size=as.numeric(input$legend_label_size_mds))
  #
  #   HTML(paste(paste0("data_mds = table_subset(taxa_table = tab1, metadata = metadata1, stratify_by_metadata = \"",input$stratify_by_metadata_mds,"\", stratify_by_value = \"",input$stratify_by_value_mds,"\", one_level = ",as.logical(one_level_all()),")"),
  #         paste0("mds_plot(taxa_table = data_mds, metadata = metadata1, test_metadata = \"",input$test_metadata_mds,"\", taxa_level = \"",input$taxa_level_mds,"\", method_mds = ",input$method_mds,"\", one_level = ",as.logical(one_level_all()),", log_norm = ",input$log_normalization_mds,
  #                 ", palette_group=c(\"",paste(strsplit(input$palette_group_mds, ",\\s*")[[1]],collapse = "\", \""),"\"), distance_type = \"",input$distance_type,"\", dot_transparency = ",as.numeric(input$dot_transparency_mds),", dot_size = ",as.numeric(input$dot_size_mds),
  #                 ", show_sample_name = ",as.logical(input$show_sample_name_mds),", ellipse_label_size = ",as.numeric(input$ellipse_label_size_mds),", label_size = ",as.numeric(input$label_size_mds),", legend_label_size = ",as.numeric(input$legend_label_size_mds),")"),
  #         sep='<br/> <br/>'))
  # })
  #
  # #alpha-div show code run
  # output$function_alpha <- renderUI({
  #   if(input$test_metadata_order_alpha[1]=="default"){
  #     test_metadata_order_alpha1="\"default\""
  #   }else{
  #     test_metadata_order_alpha1=paste0("c(\"",paste(input$test_metadata_order_alpha,collapse = "\", \""),"\")")
  #   }
  #   #dataf2=table_subset(taxa_table = data_raw(),metadata=data_meta(),stratify_by_metadata=input$stratify_by_metadata_alpha,stratify_by_value=input$stratify_by_value_alpha,one_level=as.logical(one_level_all()))
  #   #alpha_plot(taxa_table =dataf2,metadata=data_meta(),test_metadata=input$test_metadata_alpha,one_level=as.logical(one_level_all()),test_metadata_order=input$test_metadata_order_alpha,method=input$method_alpha,xlab_direction=as.integer(input$x_dir_alpha),palette_group=strsplit(input$palette_group_alpha, ",\\s*")[[1]],xlab=input$xlab_alpha)
  #   HTML(paste(paste0("data_alpha = table_subset(taxa_table = tab1, metadata = metadata1, stratify_by_metadata = \"",input$stratify_by_metadata_alpha,"\", stratify_by_value = \"",input$stratify_by_value_alpha,"\", one_level = ",as.logical(one_level_all()),")"),
  #              paste0("alpha_plot(taxa_table = data_alpha, metadata = metadata1, test_metadata = \"",input$test_metadata_alpha,"\", one_level = ",as.logical(one_level_all()),", test_metadata_order = ",test_metadata_order_alpha1,", method = \"",input$method_alpha,"\", xlab_direction = ",as.integer(input$x_dir_alpha),
  #                     ", palette_group=c(\"",paste(strsplit(input$palette_group_alpha, ",\\s*")[[1]],collapse = "\", \""),"\"), xlab = \"",input$xlab_alpha,"\")"),
  #              sep='<br/> <br/>'))
  # })
  #



  #metadata input
  data_meta <- eventReactive(input$button_meta,{
    if(input$sep_meta=="tab"){
      sep_char_meta="\t"
    }else{
      sep_char_meta=input$sep_meta
    }
    meta_format(metadata=input$file_meta$datapath,metadata_sep=sep_char_meta,meta_sample_name_col=input$meta_sample_name_col)
  })


  output$head_meta <- renderTable({
    head(data_meta(), input$n_meta)
  },rownames = TRUE)

  observeEvent(input$button_meta, {
    # Show a modal when the button is pressed
    names_len=length(intersect(rownames(data_meta()),colnames(data_raw())))
    if( names_len <1){
      showNotification(paste("There are",names_len,"matched sample names between the data table and the metadata table."),duration = 5, type = "error")
    }
    if( names_len > 2){
      showNotification(paste("There are",names_len,"matched sample names between the data table and the metadata table."),duration = 5, type = "message")
    }
  })



  # output metadata variables for test
  observe({updateSelectInput(session, "test_metadata_mds",
                    choices = colnames(data_meta()),
                    selected = c(""))})
  output$head_test_mds <- renderText({
    head(unique(na.omit(data_meta()[,input$test_metadata_mds])),n=15)
  })

  observe({updateSelectInput(session, "test_metadata_alpha",
                             choices = colnames(data_meta()),
                             selected = c(""))})

  output$head_test_alpha <- renderText({
    head(unique(na.omit(data_meta()[,input$test_metadata_alpha])),n=15)
  })

  observe({updateSelectInput(session, "test_metadata_bar",
                             choices = colnames(data_meta()),
                             selected = c(""))})
  output$head_test_bar <- renderText({
    head(unique(na.omit(data_meta()[,input$test_metadata_bar])),n=15)
  })

  observe({updateSelectInput(session, "test_metadata_stat",
                             choices = colnames(data_meta()),
                             selected = c(""))})
  output$head_test_stat <- renderText({
    head(unique(na.omit(data_meta()[,input$test_metadata_stat])),n=15)
  })


  observe({updateSelectInput(session, "col_metadata_cor",
                             choices = colnames(data_meta()),
                             selected = c(""))})

  output$head_metadata_cor_col <- renderText({
    head(unique(na.omit(data_meta()[,input$col_metadata_cor])),n=15)
  })



  # output metadata variables for correlation
  observe({updateSelectInput(session, "test_metadata_cor",
                             choices = colnames(data_meta()),
                             selected = c(""))})
  output$head_metadata_cor_test <- renderText({
    head(unique(na.omit(data_meta()[,input$test_metadata_cor])),n=15)
  })

  # output metadata variables for stratify
  observe({updateSelectInput(session, "stratify_by_metadata_filter",
                             choices = c("none",colnames(data_meta())),
                             selected = c("none"))})
  output$head_stratify_metadata_filter <- renderText({
    if (input$stratify_by_metadata_filter=="none"){
      ""
    }else{
      head(unique(na.omit(data_meta()[,input$stratify_by_metadata_filter])),n=15)
    }
  })

  observe({updateSelectInput(session, "stratify_by_metadata_mds",
                             choices = c("none",colnames(data_meta())),
                             selected = c("none"))})

  output$head_stratify_mds <- renderText({
    if (input$stratify_by_metadata_mds=="none"){
      ""
    }else{
      head(unique(na.omit(data_meta()[,input$stratify_by_metadata_mds])),n=15)
    }
  })


  stratify_by_value_mds_outVar = eventReactive(input$stratify_by_metadata_mds,{
    if (input$stratify_by_metadata_mds=="none"){
      ""
    }else{
      unique(na.omit(data_meta()[,input$stratify_by_metadata_mds]))
    }
  })

  observe({
    updateSelectInput(session, "stratify_by_value_mds",
                     choices = stratify_by_value_mds_outVar()
    )})

  stratify_by_value_filter_outVar = eventReactive(input$stratify_by_metadata_filter,{
    if (input$stratify_by_metadata_filter=="none"){
      ""
    }else{
      unique(na.omit(data_meta()[,input$stratify_by_metadata_filter]))
    }
  })

  observe({
    updateSelectInput(session, "stratify_by_value_filter",
                      choices = stratify_by_value_filter_outVar()
    )})



  observe({updateSelectInput(session, "stratify_by_metadata_alpha",
                             choices = c("none",colnames(data_meta())),
                             selected = c("none"))})
  output$head_stratify_alpha <- renderText({
    if (input$stratify_by_metadata_alpha=="none"){
      ""
    }else{
      head(unique(na.omit(data_meta()[,input$stratify_by_metadata_alpha])),n=15)
    }
  })


  stratify_by_value_alpha_outVar = eventReactive(input$stratify_by_metadata_alpha,{
    if (input$stratify_by_metadata_alpha=="none"){
      ""
    }else{
      unique(na.omit(data_meta()[,input$stratify_by_metadata_alpha]))
    }
  })


  observe({
    updateSelectInput(session, "stratify_by_value_alpha",
                      choices = stratify_by_value_alpha_outVar()
    )})


  observe({updateSelectInput(session, "stratify_by_metadata_bar",
                             choices = c("none",colnames(data_meta())),
                             selected = c("none")
    )})

  output$head_stratify_bar <- renderText({
    if (input$stratify_by_metadata_bar=="none"){
      ""
    }else{
      head(unique(na.omit(data_meta()[,input$stratify_by_metadata_bar])),n=30)
    }
  })

  stratify_by_value_bar_outVar = eventReactive(input$stratify_by_metadata_bar,{
    if (input$stratify_by_metadata_bar=="none"){
      ""
    }else{
      head(unique(na.omit(data_meta()[,input$stratify_by_metadata_bar])),n=30)
    }
  })

  observe({
    updateSelectInput(session, "stratify_by_value_bar",
                      choices = stratify_by_value_bar_outVar()
    )})

  stratify_by_value_p1_outVar = eventReactive(input$stratify_by_metadata_p1,{
    if (input$stratify_by_metadata_p1=="none"){
      ""
    }else{
      unique(na.omit(data_meta()[,input$stratify_by_metadata_p1]))
    }
  })

  observe({
    updateSelectInput(session, "stratify_by_value_p1",
                      choices = stratify_by_value_p1_outVar()
    )})

  stratify_by_value_p2_outVar = eventReactive(input$stratify_by_metadata_p2,{
    if (input$stratify_by_metadata_p2=="none"){
      ""
    }else{
      unique(na.omit(data_meta()[,input$stratify_by_metadata_p2]))
    }
  })

  observe({
    updateSelectInput(session, "stratify_by_value_p2",
                      choices = stratify_by_value_p2_outVar()
    )})

  #multiple selections of metadata variables to change the order of metadata in figures

  test_metadata_order_alpha_default = eventReactive(input$test_metadata_order_alpha_checkbox,{
    input$test_metadata_order_alpha_checkbox
  })


  test_metadata_order_alpha_outVar = eventReactive({
    input$test_metadata_alpha
    input$test_metadata_order_alpha_checkbox
    },{
    if(test_metadata_order_alpha_default ()){
      "default"
    }else{
      unique(na.omit(data_meta()[,input$test_metadata_alpha]))
    }
  })

  observe({
    updateSelectInput(session, "test_metadata_order_alpha",
                      choices = test_metadata_order_alpha_outVar()
  )})


  #multiple selections of metadata variables to change the order of metadata in figures
  #add default as an option

  test_metadata_order_bar_default = eventReactive(input$test_metadata_order_bar_checkbox,{
    input$test_metadata_order_bar_checkbox
  })


  test_metadata_order_bar_outVar = eventReactive({
    input$test_metadata_bar
    input$test_metadata_order_bar_checkbox
  },{
    if(test_metadata_order_bar_default()){
      "default"
    }else{
      unique(na.omit(data_meta()[,input$test_metadata_bar]))
    }
  })

  observe({
    updateSelectInput(session, "test_metadata_order_bar",
                      choices = test_metadata_order_bar_outVar()
  )})

  test_metadata_order_box_default = eventReactive(input$test_metadata_order_box_checkbox,{
    input$test_metadata_order_box_checkbox
  })


  test_metadata_order_box_outVar = eventReactive({
    input$test_metadata_stat
    input$test_metadata_order_box_checkbox
  },{
    if(test_metadata_order_box_default()){
      "default"
    }else{
      unique(na.omit(data_meta()[,input$test_metadata_stat]))
    }
  })


  observe({
    updateSelectInput(session, "test_metadata_order_box",
                      choices = test_metadata_order_box_outVar()
  )})

  #show metadata for PvP
  observe({updateSelectInput(session, "test_metadata_p1",
                             choices = colnames(data_meta()),
                             selected = c(""))})
  output$head_test_p1 <- renderText({
    head(unique(na.omit(data_meta()[,input$test_metadata_p1])),n=15)
  })

  observe({updateSelectInput(session, "stratify_by_metadata_p1",
                             choices = c("none",colnames(data_meta())),
                             selected = c("none"))})
  output$head_stratify_p1 <- renderText({
    if (input$stratify_by_metadata_p1=="none"){
      ""
    }else{
      head(unique(na.omit(data_meta()[,input$stratify_by_metadata_p1])),n=15)
    }
  })

  observe({updateSelectInput(session, "test_metadata_p2",
                             choices = colnames(data_meta()),
                             selected = c(""))})
  output$head_test_p2 <- renderText({
    head(unique(na.omit(data_meta()[,input$test_metadata_p2])),n=15)
  })

  observe({updateSelectInput(session, "stratify_by_metadata_p2",
                             choices = c("none",colnames(data_meta())),
                             selected = c("none"))})
  output$head_stratify_p2 <- renderText({
    if (input$stratify_by_metadata_p2=="none"){
      ""
    }else{
      head(unique(na.omit(data_meta()[,input$stratify_by_metadata_p2])),n=15)
    }
  })

  #MDS plot
  output$taxa_level_mdsUI <- renderUI({
    if (input$data_type!="Table without taxonomic structure") {
      selectInput("taxa_level_mds", "Select the taxonomic level shown (Kraken2 does not have strain level)",c("phylum","class","order","family","genus","species","ASV_or_strain"),selected="genus")
    } else {
      return(NULL)
    }
  })

  plotMDS <- eventReactive(input$button_mds,{
    dataf1=table_subset(taxa_table = data_raw(),metadata=data_meta(),stratify_by_metadata=input$stratify_by_metadata_mds,stratify_by_value=input$stratify_by_value_mds,one_level=as.logical(one_level_all()))
    mds_plot(taxa_table =dataf1,metadata=data_meta(),test_metadata=input$test_metadata_mds,taxa_level=input$taxa_level_mds,method_mds=input$method_mds,one_level=as.logical(one_level_all()),log_norm=as.logical(input$log_normalization_mds),palette_group=strsplit(input$palette_group_mds, ",\\s*")[[1]],distance_type=input$distance_type,dot_transparency=as.numeric(input$dot_transparency_mds),dot_size=as.numeric(input$dot_size_mds),show_sample_name=as.logical(input$show_sample_name_mds),ellipse_label_size=as.numeric(input$ellipse_label_size_mds),label_size=as.numeric(input$label_size_mds),legend_label_size=as.numeric(input$legend_label_size_mds))
  })

  plotMDS1 <- function(){
    dataf1=table_subset(taxa_table = data_raw(),metadata=data_meta(),stratify_by_metadata=input$stratify_by_metadata_mds,stratify_by_value=input$stratify_by_value_mds,one_level=as.logical(one_level_all()))
    mds_plot(taxa_table =dataf1,metadata=data_meta(),test_metadata=input$test_metadata_mds,taxa_level=input$taxa_level_mds,method_mds=input$method_mds,one_level=as.logical(one_level_all()),log_norm=as.logical(input$log_normalization_mds),palette_group=strsplit(input$palette_group_mds, ",\\s*")[[1]],distance_type=input$distance_type,dot_transparency=as.numeric(input$dot_transparency_mds),dot_size=as.numeric(input$dot_size_mds),show_sample_name=as.logical(input$show_sample_name_mds),ellipse_label_size=as.numeric(input$ellipse_label_size_mds),label_size=as.numeric(input$label_size_mds),legend_label_size=as.numeric(input$legend_label_size_mds))
  }

  output$plotMDS <- renderPlot({

    plotMDS()

  })

  #make figure height and width dynamic
  mds_h <- reactive({
    req(input$mds_h)
    as.numeric(input$mds_h)
  })

  mds_w <- reactive({
    req(input$mds_w)
    as.numeric(input$mds_w)
  })

  output$plotMDS.ui <- renderUI({
    if(input$method_mds=="pcoa"){
      plotOutput("plotMDS", height = 1200*mds_h(),width=400*mds_w())
    }else{
      plotOutput("plotMDS", height = 400*mds_h(),width=400*mds_w())
    }
  })

  output$plotMDSDownload <- downloadHandler(
    filename = "MDS.pdf",
    content = function(file) {
      if(input$method_mds=="pcoa"){
        pdf(file, height = 18*mds_h(),width=6*mds_w(),onefile=T)
      }else{
        pdf(file, height = 6*mds_h(),width=6*mds_w())
      }
      plotMDS1()
      dev.off()
    },contentType = "image/pdf")

  #Alpha diversity plot

  #hide option for metadata order if it is default
  output$test_metadata_order_alpha <- renderUI({
    if (input$method_stat%in%c("lr","glm","lme")) {
      selectInput(inputId = "glm_anova", "Run ANOVA on linear models?",c("False","True"),selected="False")
    } else {
      return(NULL)
    }
  })

  plotAlpha <- eventReactive(input$button_alpha,{
    dataf2=table_subset(taxa_table = data_raw(),metadata=data_meta(),stratify_by_metadata=input$stratify_by_metadata_alpha,stratify_by_value=input$stratify_by_value_alpha,one_level=as.logical(one_level_all()))
    alpha_plot(taxa_table =dataf2,metadata=data_meta(),test_metadata=input$test_metadata_alpha,one_level=as.logical(one_level_all()),test_metadata_order=input$test_metadata_order_alpha,method=input$method_alpha,xlab_direction=as.integer(input$x_dir_alpha),palette_group=strsplit(input$palette_group_alpha, ",\\s*")[[1]],xlab=input$xlab_alpha)
  })

  output$plotAlpha <- renderPlot({

    plotAlpha()

  })

  #make figure height and width dynamic
  alpha_h <- reactive({
    req(input$alpha_h)
    as.numeric(input$alpha_h)
  })

  alpha_w <- reactive({
    req(input$alpha_w)
    as.numeric(input$alpha_w)
  })

  output$plotAlpha.ui <- renderUI({
      plotOutput("plotAlpha", height = 2000*alpha_h(),width=800*alpha_w())
  })

  plotAlpha1 <- function(){
    dataf2=table_subset(taxa_table = data_raw(),metadata=data_meta(),stratify_by_metadata=input$stratify_by_metadata_alpha,stratify_by_value=input$stratify_by_value_alpha,one_level=as.logical(one_level_all()))
    alpha_plot(taxa_table =dataf2,metadata=data_meta(),test_metadata=input$test_metadata_alpha,one_level=as.logical(one_level_all()),test_metadata_order=input$test_metadata_order_alpha,method=input$method_alpha,xlab_direction=as.integer(input$x_dir_alpha),palette_group=strsplit(input$palette_group_alpha, ",\\s*")[[1]],xlab=input$xlab_alpha)
  }

  output$plotAlphaDownload <- downloadHandler(
    filename ="alpha_diversity.pdf",
    content = function(file) {
      pdf(file, height = 30*alpha_h(),width=15*alpha_w())
      plotAlpha1()
      dev.off()
    },contentType = "image/pdf")


  data_alpha <- eventReactive(input$button_alpha,{
    alpha_data(taxa_table =data_raw(),one_level=as.logical(one_level_all()))
  })


  output$downloadAlpha <- downloadHandler(
    filename = "alpha_diversity.csv",
    content = function(file) {
      write.csv(data_alpha(), file, row.names = TRUE)
    }
  )

  #Taxa barplot

  plotBar <- eventReactive(input$button_bar,{
    dataf2=table_subset(taxa_table = data_raw(),metadata=data_meta(),stratify_by_metadata=input$stratify_by_metadata_bar,stratify_by_value=input$stratify_by_value_bar,one_level=as.logical(one_level_all()))
    taxa_barplot(taxa_table =dataf2,metadata=data_meta(),test_metadata=input$test_metadata_bar,one_level=as.logical(one_level_all()),num_taxa=as.integer(input$num_taxa_bar),test_metadata_order=input$test_metadata_order_bar,taxa_level=input$taxa_level_bar,xlab_direction=as.integer(input$x_dir_bar),legend_size=as.numeric(input$legend_size_bar),palette_group=strsplit(input$palette_group_bar, ",\\s*")[[1]])
  })

  output$plotBar <- renderPlot({

    plotBar()

  })

  #make figure height and width dynamic
  bar_h <- reactive({
    req(input$bar_h)
    as.numeric(input$bar_h)
  })

  bar_w <- reactive({
    req(input$bar_w)
    as.numeric(input$bar_w)
  })


  output$plotBar.ui <- renderUI({
    plotOutput("plotBar", height = 1000*bar_h(),width=1000*bar_w())
  })


  plotBar1 <- function(){
    dataf2=table_subset(taxa_table = data_raw(),metadata=data_meta(),stratify_by_metadata=input$stratify_by_metadata_bar,stratify_by_value=input$stratify_by_value_bar,one_level=as.logical(one_level_all()))
    taxa_barplot(taxa_table =dataf2,metadata=data_meta(),test_metadata=input$test_metadata_bar,one_level=as.logical(one_level_all()),num_taxa=as.integer(input$num_taxa_bar),test_metadata_order=input$test_metadata_order_bar,taxa_level=input$taxa_level_bar,xlab_direction=as.integer(input$x_dir_bar),legend_size=as.numeric(input$legend_size_bar),palette_group=strsplit(input$palette_group_bar, ",\\s*")[[1]])
   }

  output$plotBarDownload <- downloadHandler(
    filename = "barplot.pdf",
    content = function(file) {
      pdf(file, height = 15*bar_h(),width=15*bar_w())
      plotBar1()
      dev.off()
    },contentType = "image/pdf")

  data_bar <- eventReactive(input$button_bar,{
    dataf2=table_subset(taxa_table = data_raw(),metadata=data_meta(),stratify_by_metadata=input$stratify_by_metadata_bar,stratify_by_value=input$stratify_by_value_bar,one_level=as.logical(one_level_all()))
    taxa_bar_table(taxa_table =dataf2,metadata=data_meta(),test_metadata=input$test_metadata_bar,one_level=as.logical(one_level_all()),num_taxa=as.integer(input$num_taxa_bar),test_metadata_order=input$test_metadata_order_bar,taxa_level=input$taxa_level_bar,xlab_direction=as.integer(input$x_dir_bar),legend_size=as.numeric(input$legend_size_bar),palette_group=strsplit(input$palette_group_bar, ",\\s*")[[1]])
  })


  output$downloadBar <- downloadHandler(
    filename = "barplot_table.csv",
    content = function(file) {
      write.csv(data_bar(), file, row.names = TRUE)
    }
  )


  #data filter/subset
  data_filtered <- eventReactive(input$button_filter,{
    table_subset(taxa_table = data_raw(),metadata=data_meta(),stratify_by_metadata=input$stratify_by_metadata_filter,stratify_by_value=input$stratify_by_value_filter,prevalence_cutoff=as.numeric(input$prevalence_cutoff), abundance_cutoff=as.numeric(input$abundance_cutoff),one_level=as.logical(one_level_all()),exclude_ASV_strain=as.logical(input$exclude_ASV_filter))
   })


  output$head_filtered <- renderTable({
    if(ncol(data_filtered())<as.numeric(input$n_filtered_col)){
      n_col=ncol(data_filtered())
    }else{
      n_col=as.numeric(input$n_filtered_col)
    }
    head(data_filtered()[,1:n_col], input$n_filtered)
  },rownames = TRUE)


  output$filtered_dim <- renderText({
    paste("Number of rows :",nrow(data_filtered()),"\nNumber of columns :",ncol(data_filtered()))
  })


  #statistical test
  output$method_statUI <- renderUI({
    if (input$test_metadata_continuous_stat=="False") {
      selectInput(inputId = "method_stat", "Select method (lr:logistic regression; glm:generalized linear model; lme:linear mixed-effects models)",c("wilcoxon","t.test","kruskal-wallis","anova","glm","lr","lme"),selected ="wilcoxon")
    } else {
      selectInput(inputId = "method_stat", "Select method (lr:logistic regression; glm:generalized linear model; lme:linear mixed-effects models)",c("pearson","spearman","kendall","glm","lr","lme"),selected ="spearman")
    }
  })

  output$outcome_metaUI <- renderUI({
    if (input$method_stat%in%c("lr","glm","lme")) {
      selectInput("outcome_meta", "Is the metadata used as the outcome? Select True for lr test.",c("False","True"),selected="False")
    } else {
      return(NULL)
    }
  })

  output$model_glmUI <- renderUI({
    if (input$method_stat%in%c("lr","glm","lme")) {
      textAreaInput(inputId = "model_glm", "Covariates for adjusting (e.g., +age+factor(sex)+factor(batch) or *age for including interactions). Please make sure no space in covariate names.", value = "")
    } else {
      return(NULL)
    }
  })

  output$glm_anovaUI <- renderUI({
    if (input$method_stat%in%c("lr","glm","lme")) {
      selectInput(inputId = "glm_anova", "Run ANOVA on linear models?",c("False","True"),selected="False")
    } else {
      return(NULL)
    }
  })

  output$glm_distUI <- renderUI({
    if (input$method_stat=="glm") {
      textAreaInput(inputId = "glm_dist", "Family variable of glm function (e.g.,binomial, gaussian, poisson). Default is gaussian for continuous outcomes, and binomial for two-category outcomes.", value = "default")
    } else {
      return(NULL)
    }
  })

  output$random_effect_varUI <- renderUI({
    if (input$method_stat=="lme") {
      observe({updateSelectInput(session, "random_effect_var",choices = colnames(data_meta()),selected = c(""))})
      output$random_effect_var <- renderText({
        head(unique(na.omit(data_meta()[,input$random_effect_var])),n=15)
      })
      selectInput(inputId = "random_effect_var", "Random effect variable for mixed effects models",c(""))
    } else {
      return(NULL)
    }
  })


  data_fdrs <- eventReactive(input$button_fdrs,{
    stat_test(taxa_table =data_filtered(),metadata=data_meta(),test_metadata=input$test_metadata_stat,method=input$method_stat,log_norm=as.logical(input$log_norm_stat),outcome_meta=as.logical(input$outcome_meta),test_metadata_continuous=as.logical(input$test_metadata_continuous_stat),glm_anova=as.logical(input$glm_anova),model_glm=input$model_glm,glm_dist=input$glm_dist,random_effect_var=input$random_effect_var)
   })

  data_fdrs1 <- reactive({
    if(input$sort_fdr){
      data_fdrs()[order(as.numeric(data_fdrs()[,2])),]
    }else{
      data_fdrs()
      }
  })

  output$head_fdrs <- renderTable({
    head(data_fdrs1(), input$n_fdrs)
  },rownames = TRUE,digits=-2)

  output$downloadStat <- downloadHandler(
    filename = function() {
      "stats.csv"
    },
    content = function(file) {
      write.csv(data_fdrs1(), file, row.names = TRUE)
    }
  )

  #Histogram of P values
  plotPvalsHist <- eventReactive(input$button_fdrs,{
    hist(as.numeric(data_fdrs1()[,2]),breaks=20,xlab="P values",main="")
   })

  output$plotPvalsHist <- renderPlot({

    plotPvalsHist()

  })

  output$plotPvalsHist.ui <- renderUI({
    plotOutput("plotPvalsHist", height = 600,width=600)
  })

  #Tree plot

  plotTree <- eventReactive(input$button_tree,{
    tree_view(taxa_table = data_raw(),metadata=data_meta(),fdrs=data_fdrs(),test_metadata=input$test_metadata_stat,fdr_cutoff=input$fdr_cutoff_tree,test_metadata_continuous=as.logical(input$test_metadata_continuous_tree),
              node_size_breaks=as.numeric(strsplit(input$node_size_breaks, ",\\s*")[[1]]),palette_highlight=strsplit(input$palette_group_tree, ",\\s*")[[1]],single_parent_branch_removal=as.logical(input$single_parent_branch_removal),single_child_branch_removal=as.logical(input$single_child_branch_removal),
              prevalence_cutoff=as.numeric(input$prevalence_cutoff_tree), abundance_cutoff=as.numeric(input$abundance_cutoff_tree),taxa_removal=input$taxa_removal_tree,branch=input$branch_tree)
  })

  output$plotTree <- renderPlot({

    plotTree()

  })

  #make figure height and width dynamic
  tree_h <- reactive({
    req(input$tree_h)
    as.numeric(input$tree_h)
  })

  tree_w <- reactive({
    req(input$tree_w)
    as.numeric(input$tree_w)
  })


  output$plotTree.ui <- renderUI({
    plotOutput("plotTree", height = 600*tree_h(),width=900*tree_w())
  })


  plotTree1 <- function(){
    tree_view(taxa_table = data_raw(),metadata=data_meta(),fdrs=data_fdrs(),test_metadata=input$test_metadata_stat,fdr_cutoff=input$fdr_cutoff_tree,test_metadata_continuous=as.logical(input$test_metadata_continuous_tree),
              node_size_breaks=as.numeric(strsplit(input$node_size_breaks, ",\\s*")[[1]]),palette_highlight=strsplit(input$palette_group_tree, ",\\s*")[[1]],single_parent_branch_removal=as.logical(input$single_parent_branch_removal),single_child_branch_removal=as.logical(input$single_child_branch_removal),
              prevalence_cutoff=as.numeric(input$prevalence_cutoff_tree), abundance_cutoff=as.numeric(input$abundance_cutoff_tree),taxa_removal=input$taxa_removal_tree)
    }

  output$plotTreeDownload <- downloadHandler(
    filename = "tree.pdf",
    content = function(file) {
      pdf(file,height=9*tree_h(),width=12*tree_w())
      print(plotTree1())
      dev.off()
    })

  observeEvent(input$button_tree, {
    #notification when there is one level
    if(as.logical(one_level_all())){
      showNotification(paste("The taxa table does not have a taxonomic structure that can be used to build tree. Please check https://github.com/ssun6/plotmicrobiome for sample files."),duration = 20, type = "error")
    }
  })


  #Box plot

  plotBox <- eventReactive(input$button_box,{
    taxa_boxplot(taxa_table =data_filtered(),metadata=data_meta(),fdrs=data_fdrs(),test_metadata=input$test_metadata_stat,test_metadata_order=input$test_metadata_order_box,one_level=as.logical(one_level_all()),
                 log_norm=input$log_norm_box,taxa_shown=input$taxa_shown_box,cutoff=input$fdr_cutoff_box,page=input$page_box,xlab=input$xlab_box,ylab=input$ylab_box,
                 xlab_direction=as.integer(input$x_dir_box),palette_group=strsplit(input$palette_group_box, ",\\s*")[[1]])
  })

  output$plotBox <- renderPlot({

    plotBox()

  })

  #make figure height and width dynamic
  box_h <- reactive({
    req(input$box_h)
    as.numeric(input$box_h)
  })

  box_w <- reactive({
    req(input$box_w)
    as.numeric(input$box_w)
  })

  output$plotBox.ui <- renderUI({
    plotOutput("plotBox", height = 1000*box_h(),width=1000*box_w())
  })


  plotBox1 <- function(){
    taxa_boxplot_download(taxa_table =data_filtered(),metadata=data_meta(),fdrs=data_fdrs(),test_metadata=input$test_metadata_stat,test_metadata_order=input$test_metadata_order_box,one_level=as.logical(one_level_all()),
                 log_norm=input$log_norm_box,taxa_shown=input$taxa_shown_box,cutoff=input$fdr_cutoff_box,xlab=input$xlab_box,ylab=input$ylab_box,
                 xlab_direction=as.integer(input$x_dir_box),palette_group=strsplit(input$palette_group_box, ",\\s*")[[1]])
  }

  output$plotBoxDownload <- downloadHandler(
    filename = "boxplots.pdf",
    content = function(file) {
      pdf(file,height=10*box_h(),width=10*box_w(),onefile = T)
      plotBox1()
      dev.off()
    })


  #correlation plot

  plotCor <- eventReactive(input$button_cor,{
    meta_corplot(taxa_table =data_filtered(),metadata=data_meta(),test_metadata=input$test_metadata_cor,cor_method=input$cor_method,one_level=as.logical(one_level_all()),
                 col_metadata=input$col_metadata_cor,log_norm=input$log_norm_cor,taxa_shown=input$taxa_shown_cor,page=input$page_cor,ylab=input$ylab_cor,
                 fdr_cutoff=input$fdr_cutoff_cor,palette_group=strsplit(input$palette_group_cor, ",\\s*")[[1]])
  })


  output$plotCor <- renderPlot({

    plotCor()

  })

  #make figure height and width dynamic
  cor_h <- reactive({
    req(input$cor_h)
    as.numeric(input$cor_h)
  })

  cor_w <- reactive({
    req(input$cor_w)
    as.numeric(input$cor_w)
  })

  output$plotCor.ui <- renderUI({
    plotOutput("plotCor", height = 1000*cor_h(),width=1000*cor_w())
  })


  plotCor1 <- function(){
    meta_corplot_download(taxa_table =data_filtered(),metadata=data_meta(),test_metadata=input$test_metadata_cor,cor_method=input$cor_method,one_level=as.logical(one_level_all()),
                 col_metadata=input$col_metadata_cor,log_norm=input$log_norm_cor,taxa_shown=input$taxa_shown_cor,ylab=input$ylab_cor,
                 fdr_cutoff=input$fdr_cutoff_cor,palette_group=strsplit(input$palette_group_cor, ",\\s*")[[1]])
  }

  output$plotCorDownload <- downloadHandler(
    filename = "correlation.pdf",
    content = function(file) {
      pdf(file,height=12*cor_h(),width=12*cor_w(),onefile = T)
      print(plotCor1())
      dev.off()
    })


  data_p1 <- eventReactive(input$button_p1_2,{
    if(input$upload_p1=="Yes"){
      if(input$sep_p1=="tab"){
        sep_p1="\t"
      }else{
        sep_p1=input$sep_p1
      }
      read.table(file =input$file_p1$datapath,sep=input$sep_p1,row.names=1,header=T,check.names = F)
    }else{
      dataf1_p1=table_subset(taxa_table = data_raw(),metadata=data_meta(),stratify_by_metadata=input$stratify_by_metadata_p1,stratify_by_value=input$stratify_by_value_p1,prevalence_cutoff=as.numeric(input$prevalence_cutoff_p1), abundance_cutoff=as.numeric(input$abundance_cutoff_p1),one_level=as.logical(one_level_all()),exclude_ASV_strain=as.logical(input$exclude_ASV_p1))
      stat_test(taxa_table =dataf1_p1,metadata=data_meta(),test_metadata=input$test_metadata_p1,method=input$method_stat_p1,test_metadata_continuous=as.logical(input$test_metadata_continuous_p1))
    }
  })

  output$head_p1 <- renderTable({
    head(data_p1(), 5)
  },rownames = TRUE)

  data_p2 <- eventReactive(input$button_p2_2,{
    if(input$upload_p2=="Yes"){
      if(input$sep_p2=="tab"){
        sep_p2="\t"
      }else{
        sep_p2=input$sep_p2
      }
      read.table(file =input$file_p2$datapath,sep=input$sep_p2,row.names=1,header=T,check.names = F)
    }else{
      dataf1_p2=table_subset(taxa_table = data_raw(),metadata=data_meta(),stratify_by_metadata=input$stratify_by_metadata_p2,stratify_by_value=input$stratify_by_value_p2,prevalence_cutoff=as.numeric(input$prevalence_cutoff_p2), abundance_cutoff=as.numeric(input$abundance_cutoff_p2),one_level=as.logical(one_level_all()),exclude_ASV_strain=as.logical(input$exclude_ASV_p2))
      stat_test(taxa_table =dataf1_p2,metadata=data_meta(),test_metadata=input$test_metadata_p2,method=input$method_stat_p2,test_metadata_continuous=as.logical(input$test_metadata_continuous_p2))
    }
  })

  output$head_p2 <- renderTable({
    head(data_p2(), 5)
  },rownames = TRUE)

  plotPvals <- eventReactive(input$button_cor_p,{
    p_compare(data_p1(),data_p2(),p_col1=as.numeric(input$p1_col),p_col2=as.numeric(input$p2_col),indicator1=as.numeric(input$ind1_col),indicator2=as.numeric(input$ind2_col),point_color=input$point_color,lab_cutoff=as.numeric(input$lab_cutoff),cor_method=input$cor_method_p,x.reverse = as.logical(input$x_reverse),y.reverse = as.logical(input$y_reverse),exclude_unclassified=as.logical(input$exclude_unclassified),one_level=as.logical(one_level_all()),direction=as.logical(input$direction_pvp))
  })

  output$plotPvals <- renderPlot({

    plotPvals()

  })

  #make figure height and width dynamic
  pvp_h <- reactive({
    req(input$pvp_h)
    as.numeric(input$pvp_h)
  })

  pvp_w <- reactive({
    req(input$pvp_w)
    as.numeric(input$pvp_w)
  })

  output$plotPvals.ui <- renderUI({
    plotOutput("plotPvals", height = 800*pvp_h(),width=800*pvp_w())
  })

  plotPvals1 <- function(){
    p_compare(data_p1(),data_p2(),p_col1=as.numeric(input$p1_col),p_col2=as.numeric(input$p2_col),indicator1=as.numeric(input$ind1_col),indicator2=as.numeric(input$ind2_col),point_color=input$point_color,lab_cutoff=as.numeric(input$lab_cutoff),cor_method=input$cor_method_p,x.reverse = as.logical(input$x_reverse),y.reverse = as.logical(input$y_reverse),exclude_unclassified=as.logical(input$exclude_unclassified),one_level=as.logical(one_level_all()),direction=as.logical(input$direction_pvp))
  }

  output$plotPvalsDownload <- downloadHandler(
    filename = "PvsP.pdf",
    content = function(file) {
      pdf(file, height = 10*pvp_h(),width=10*pvp_w())
      plotPvals1()
      dev.off()
    },contentType = "image/pdf")

  #add link to code
  #16S file format code
  url_code_link_16s <- a("Link to code of function", href="https://github.com/ssun6/plotmicrobiome/blob/main/R/format_asv.R")
  output$code_link_16s <- renderUI({
    tagList(url_code_link_16s)
  })

  #wgs file format code
  url_code_link_wgs <- a("Link to code of function", href="https://github.com/ssun6/plotmicrobiome/blob/main/R/format_wgs.R")
  output$code_link_wgs <- renderUI({
    tagList(url_code_link_wgs)
  })

  #one level file format code
  url_code_link_tab <- a("Link to code of function", href="https://github.com/ssun6/plotmicrobiome/blob/main/R/format_tabs.R")
  output$code_link_path <- renderUI({
    tagList(url_code_link_path)
  })

  #metadata file format code
  url_code_link_meta <- a("Link to code of function", href="https://github.com/ssun6/plotmicrobiome/blob/main/R/meta_format.R")
  output$code_link_meta <- renderUI({
    tagList(url_code_link_meta)
  })

  #mds plot code
  url_code_link_mds <- a("Link to code of function", href="https://github.com/ssun6/plotmicrobiome/blob/main/R/mds_plot.R")
  output$code_link_mds <- renderUI({
    tagList(url_code_link_mds)
  })

  #alpha div code
  url_code_link_alpha <- a("Link to code of function", href="https://github.com/ssun6/plotmicrobiome/blob/main/R/alpha_plot.R")
  output$code_link_alpha <- renderUI({
    tagList(url_code_link_alpha)
  })

  #barplot code
  url_code_link_bar <- a("Link to code of function", href="https://github.com/ssun6/plotmicrobiome/blob/main/R/taxa_barplot.R")
  output$code_link_bar <- renderUI({
    tagList(url_code_link_bar)
  })

  #subset table code
  url_code_link_subset <- a("Link to code of function", href="https://github.com/ssun6/plotmicrobiome/blob/main/R/table_subset.R")
  output$code_link_subset <- renderUI({
    tagList(url_code_link_subset)
  })

  #statistical test code
  url_code_link_stat <- a("Link to code of function", href="https://github.com/ssun6/plotmicrobiome/blob/main/R/stat_test.R")
  output$code_link_stat <- renderUI({
    tagList(url_code_link_stat)
  })

  #boxplot code
  url_code_link_box <- a("Link to code of function 1 (display figures by page)", href="https://github.com/ssun6/plotmicrobiome/blob/main/R/taxa_boxplot.R")
  output$code_link_box <- renderUI({
    tagList(url_code_link_box)
  })

  url_code_link_box_download <- a("Link to code of function 2 (all figures in one file)", href="https://github.com/ssun6/plotmicrobiome/blob/main/R/taxa_boxplot_download.R")
  output$code_link_box_download <- renderUI({
    tagList(url_code_link_box_download)
  })

  #tree plot code
  url_code_link_tree <- a("Link to code of function", href="https://github.com/ssun6/plotmicrobiome/blob/main/R/tree_view.R")
  output$code_link_tree <- renderUI({
    tagList(url_code_link_tree)
  })

  #scatter plot code
  url_code_link_cor <- a("Link to code of function 1 (display figures by page)", href="https://github.com/ssun6/plotmicrobiome/blob/main/R/meta_corplot.R")
  output$code_link_cor <- renderUI({
    tagList(url_code_link_cor)
  })

  url_code_link_cor_download <- a("Link to code of function 2 (all figures in one file)", href="https://github.com/ssun6/plotmicrobiome/blob/main/R/meta_corplot_download.R")
  output$code_link_cor_download <- renderUI({
    tagList(url_code_link_cor_download)
  })

  #P-value vs P-value plot
  url_code_link_pvp <- a("Link to code of function", href="https://github.com/ssun6/plotmicrobiome/blob/main/R/p_compare.R")
  output$code_link_pvp <- renderUI({
    tagList(url_code_link_pvp)
  })

  #code download
  #data input
  output$code_datainput_Download <- downloadHandler(
    filename = "data_input_code.txt",
    content = function(file) {
      code_library_text=paste0("library(plotmicrobiome)")


      if(input$data_type == "Amplicon with taxonomic structure"){
        if(input$sep_16s=="tab"){
          sep_char_16s="\\t"
        }else{
          sep_char_16s=input$sep_16s
        }
        code_input_text=paste0("tab1 = format_asv(taxa_file = \"",input$file_16s$name,"\", biom = ",input$biom,", ASV = ",input$ASV,", sep = \"",sep_char_16s,"\", reads_cutoff = ",input$n_reads_16s,", normalization = ",input$norm_16s,", rarefy = ",input$rarefy_16s,", rarefy_num = ",input$rarefy_reads_16s,")")
      }else if(input$data_type == "Metagenomics with taxonomic structure"){
        if(input$sep_wgs=="tab"){
          sep_char_wgs="\\t"
        }else{
          sep_char_wgs=input$sep_wgs
        }
        code_input_text=paste0("tab1 = format_wgs(taxa_file = \"",input$file_wgs$name,"\", sep = \"",sep_char_wgs,"\", reads_cutoff = ",input$n_reads_wgs,", method = \"",input$method_wgs,"\", normalization = ",input$norm_wgs,", rarefy = ",input$rarefy_wgs,", rarefy_num = ",input$rarefy_reads_wgs,")")
      }else{
        if(input$sep_tab=="tab"){
          sep_char_tab="\\t"
        }else{
          sep_char_tab=input$sep_tab
        }
        code_input_text=paste0("tab1 = format_tabs(taxa_file = \"",input$file_tab$name,"\", sep = \"",sep_char_tab,"\", reads_cutoff = ",input$n_reads_tab,", normalization = ",input$norm_tab,", rarefy = ",input$rarefy_tab,", rarefy_num = ",input$rarefy_reads_tab,")")
      }

      line_datainput=paste(code_library_text,code_input_text,sep="\n")
      writeLines(line_datainput, file)
    }
  )

  #metdadata input
  output$code_metadata_Download <- downloadHandler(
    filename = "metadata_input_code.txt",
    content = function(file) {
      code_library_text=paste0("library(plotmicrobiome)")
      if(input$sep_meta=="tab"){
        sep_char_meta="\\t"
      }else{
        sep_char_meta=input$sep_meta
      }
      code_meta_text=paste0("metadata1 = meta_format(metadata = \"",input$file_meta$name,"\", metadata_sep = \"",sep_char_meta,"\", meta_sample_name_col = ",input$meta_sample_name_col,")")
      line_meta_input=paste(code_library_text,code_meta_text,sep="\n")
      writeLines(line_meta_input, file)
    }
  )


  #MDS
  output$codeMDSDownload <- downloadHandler(
    filename = "MDS_code.txt",
    content = function(file) {
      code_library_text=paste0("library(plotmicrobiome)")


      if(input$data_type == "Amplicon with taxonomic structure"){
        if(input$sep_16s=="tab"){
          sep_char_16s="\\t"
        }else{
          sep_char_16s=input$sep_16s
        }
        code_input_text=paste0("tab1 = format_asv(taxa_file = \"",input$file_16s$name,"\", biom = ",input$biom,", ASV = ",input$ASV,", sep = \"",sep_char_16s,"\", reads_cutoff = ",input$n_reads_16s,", normalization = ",input$norm_16s,", rarefy = ",input$rarefy_16s,", rarefy_num = ",input$rarefy_reads_16s,")")
      }else if(input$data_type == "Metagenomics with taxonomic structure"){
        if(input$sep_wgs=="tab"){
          sep_char_wgs="\\t"
        }else{
          sep_char_wgs=input$sep_wgs
        }
        code_input_text=paste0("tab1 = format_wgs(taxa_file = \"",input$file_wgs$name,"\", sep = \"",sep_char_wgs,"\", reads_cutoff = ",input$n_reads_wgs,", method = \"",input$method_wgs,"\", normalization = ",input$norm_wgs,", rarefy = ",input$rarefy_wgs,", rarefy_num = ",input$rarefy_reads_wgs,")")
      }else{
        if(input$sep_tab=="tab"){
          sep_char_tab="\\t"
        }else{
          sep_char_tab=input$sep_tab
        }
        code_input_text=paste0("tab1 = format_tabs(taxa_file = \"",input$file_tab$name,"\", sep = \"",sep_char_tab,"\", reads_cutoff = ",input$n_reads_tab,", normalization = ",input$norm_tab,", rarefy = ",input$rarefy_tab,", rarefy_num = ",input$rarefy_reads_tab,")")
      }

      if(input$sep_meta=="tab"){
          sep_char_meta="\\t"
      }else{
          sep_char_meta=input$sep_meta
      }
      code_meta_text=paste0("metadata1 = meta_format(metadata = \"",input$file_meta$name,"\", metadata_sep = \"",sep_char_meta,"\", meta_sample_name_col = ",input$meta_sample_name_col,")")
      code_mds_text=paste(paste0("data_mds = table_subset(taxa_table = tab1, metadata = metadata1, stratify_by_metadata = \"",input$stratify_by_metadata_mds,"\", stratify_by_value = \"",input$stratify_by_value_mds,"\", one_level = ",as.logical(one_level_all()),")"),
                                  paste0("mds_plot(taxa_table = data_mds, metadata = metadata1, test_metadata = \"",input$test_metadata_mds,"\", taxa_level = \"",input$taxa_level_mds,"\", method_mds = ",input$method_mds,"\", one_level = ",as.logical(one_level_all()),", log_norm = ",input$log_normalization_mds,
                                         ", palette_group=c(\"",paste(strsplit(input$palette_group_mds, ",\\s*")[[1]],collapse = "\", \""),"\"), distance_type = \"",input$distance_type,"\", dot_transparency = ",as.numeric(input$dot_transparency_mds),", dot_size = ",as.numeric(input$dot_size_mds),
                                         ", show_sample_name = ",as.logical(input$show_sample_name_mds),", ellipse_label_size = ",as.numeric(input$ellipse_label_size_mds),", label_size = ",as.numeric(input$label_size_mds),", legend_label_size = ",as.numeric(input$legend_label_size_mds),")"),sep="\n")

      line_mds=paste(code_library_text,code_input_text,code_meta_text,code_mds_text,sep="\n")
      writeLines(line_mds, file)
    }
  )

  #alpha diversity
  output$code_alpha_Download <- downloadHandler(
    filename = "Alpha_diversity_code.txt",
    content = function(file) {
      code_library_text=paste0("library(plotmicrobiome)")


      if(input$data_type == "Amplicon with taxonomic structure"){
        if(input$sep_16s=="tab"){
          sep_char_16s="\\t"
        }else{
          sep_char_16s=input$sep_16s
        }
        code_input_text=paste0("tab1 = format_asv(taxa_file = \"",input$file_16s$name,"\", biom = ",input$biom,", ASV = ",input$ASV,", sep = \"",sep_char_16s,"\", reads_cutoff = ",input$n_reads_16s,", normalization = ",input$norm_16s,", rarefy = ",input$rarefy_16s,", rarefy_num = ",input$rarefy_reads_16s,")")
      }else if(input$data_type == "Metagenomics with taxonomic structure"){
        if(input$sep_wgs=="tab"){
          sep_char_wgs="\\t"
        }else{
          sep_char_wgs=input$sep_wgs
        }
        code_input_text=paste0("tab1 = format_wgs(taxa_file = \"",input$file_wgs$name,"\", sep = \"",sep_char_wgs,"\", reads_cutoff = ",input$n_reads_wgs,", method = \"",input$method_wgs,"\", normalization = ",input$norm_wgs,", rarefy = ",input$rarefy_wgs,", rarefy_num = ",input$rarefy_reads_wgs,")")
      }else{
        if(input$sep_tab=="tab"){
          sep_char_tab="\\t"
        }else{
          sep_char_tab=input$sep_tab
        }
        code_input_text=paste0("tab1 = format_tabs(taxa_file = \"",input$file_tab$name,"\", sep = \"",sep_char_tab,"\", reads_cutoff = ",input$n_reads_tab,", normalization = ",input$norm_tab,", rarefy = ",input$rarefy_tab,", rarefy_num = ",input$rarefy_reads_tab,")")
      }

      if(input$sep_meta=="tab"){
        sep_char_meta="\\t"
      }else{
        sep_char_meta=input$sep_meta
      }
      code_meta_text=paste0("metadata1 = meta_format(metadata = \"",input$file_meta$name,"\", metadata_sep = \"",sep_char_meta,"\", meta_sample_name_col = ",input$meta_sample_name_col,")")

      #dataf2=table_subset(taxa_table = data_raw(),metadata=data_meta(),stratify_by_metadata=input$stratify_by_metadata_alpha,stratify_by_value=input$stratify_by_value_alpha,one_level=as.logical(one_level_all()))
      #alpha_plot(taxa_table =dataf2,metadata=data_meta(),test_metadata=input$test_metadata_alpha,one_level=as.logical(one_level_all()),test_metadata_order=input$test_metadata_order_alpha,method=input$method_alpha,xlab_direction=as.integer(input$x_dir_alpha),palette_group=strsplit(input$palette_group_alpha, ",\\s*")[[1]],xlab=input$xlab_alpha)

      code_alpha_text=paste(paste0("data_alpha = table_subset(taxa_table = tab1, metadata = metadata1, stratify_by_metadata = \"",input$stratify_by_metadata_alpha,"\", stratify_by_value = \"",input$stratify_by_value_alpha,"\", one_level = ",as.logical(one_level_all()),")"),
                          paste0("alpha_plot(taxa_table = data_alpha, metadata = metadata1, test_metadata = \"",input$test_metadata_alpha,"\", method = \"",input$method_alpha,"\", one_level = ",as.logical(one_level_all()),", test_metadata_order = c(\"",paste(input$test_metadata_order_alpha,collapse = "\",\""),
                                 "\"), palette_group=c(\"",paste(strsplit(input$palette_group_alpha, ",\\s*")[[1]],collapse = "\", \""),"\")",", xlab_direction = ",as.integer(input$x_dir_alpha),", xlab = \"", input$xlab_alpha,"\")"),sep="\n")

      line_alpha=paste(code_library_text,code_input_text,code_meta_text,code_alpha_text,sep="\n")
      writeLines(line_alpha, file)
    }
  )

  #barplot
  output$code_bar_Download <- downloadHandler(
    filename = "Barplot_code.txt",
    content = function(file) {
      code_library_text=paste0("library(plotmicrobiome)")


      if(input$data_type == "Amplicon with taxonomic structure"){
        if(input$sep_16s=="tab"){
          sep_char_16s="\\t"
        }else{
          sep_char_16s=input$sep_16s
        }
        code_input_text=paste0("tab1 = format_asv(taxa_file = \"",input$file_16s$name,"\", biom = ",as.logical(input$biom),", ASV = ",as.logical(input$ASV),", sep = \"",sep_char_16s,"\", reads_cutoff = ",input$n_reads_16s,", normalization = ",as.logical(input$norm_16s),", rarefy = ",as.logical(input$rarefy_16s),", rarefy_num = ",input$rarefy_reads_16s,")")
      }else if(input$data_type == "Metagenomics with taxonomic structure"){
        if(input$sep_wgs=="tab"){
          sep_char_wgs="\\t"
        }else{
          sep_char_wgs=input$sep_wgs
        }
        code_input_text=paste0("tab1 = format_wgs(taxa_file = \"",input$file_wgs$name,"\", sep = \"",sep_char_wgs,"\", reads_cutoff = ",input$n_reads_wgs,", method = \"",input$method_wgs,"\", normalization = ",as.logical(input$norm_wgs),", rarefy = ",as.logical(input$rarefy_wgs),", rarefy_num = ",input$rarefy_reads_wgs,")")
      }else{
        if(input$sep_tab=="tab"){
          sep_char_tab="\\t"
        }else{
          sep_char_tab=input$sep_tab
        }
        code_input_text=paste0("tab1 = format_tabs(taxa_file = \"",input$file_tab$name,"\", sep = \"",sep_char_tab,"\", reads_cutoff = ",input$n_reads_tab,", normalization = ",as.logical(input$norm_tab),", rarefy = ",as.logical(input$rarefy_tab),", rarefy_num = ",input$rarefy_reads_tab,")")
      }

      if(input$sep_meta=="tab"){
        sep_char_meta="\\t"
      }else{
        sep_char_meta=input$sep_meta
      }
      code_meta_text=paste0("metadata1 = meta_format(metadata = \"",input$file_meta$name,"\", metadata_sep = \"",sep_char_meta,"\", meta_sample_name_col = ",input$meta_sample_name_col,")")

      #dataf2=table_subset(taxa_table = data_raw(),metadata=data_meta(),stratify_by_metadata=input$stratify_by_metadata_bar,stratify_by_value=input$stratify_by_value_bar,one_level=as.logical(one_level_all()))
      #taxa_barplot(taxa_table =dataf2,metadata=data_meta(),test_metadata=input$test_metadata_bar,one_level=as.logical(one_level_all()),num_taxa=as.integer(input$num_taxa_bar),test_metadata_order=input$test_metadata_order_bar,taxa_level=input$taxa_level_bar,xlab_direction=as.integer(input$x_dir_bar),legend_size=as.numeric(input$legend_size_bar),palette_group=strsplit(input$palette_group_bar, ",\\s*")[[1]])

      code_bar_text=paste(paste0("data_bar = table_subset(taxa_table = tab1, metadata = metadata1, stratify_by_metadata = \"",input$stratify_by_metadata_bar,"\", stratify_by_value = c(\"",paste(input$stratify_by_value_bar,collapse = "\", \""),"\"), one_level = ",as.logical(one_level_all()),")"),
                            paste0("taxa_barplot(taxa_table = data_bar, metadata = metadata1, test_metadata = \"",input$test_metadata_bar,"\", one_level = ",as.logical(one_level_all()),", num_taxa = ",as.integer(input$num_taxa_bar),", taxa_level = \"",input$taxa_level_bar,"\", test_metadata_order = c(\"",paste(input$test_metadata_order_bar,collapse = "\",\""),
                                   "\"), palette_group=c(\"",paste(strsplit(input$palette_group_alpha, ",\\s*")[[1]],collapse = "\", \""),"\")",", xlab_direction = ",as.integer(input$x_dir_alpha),", legend_size = ", as.numeric(input$legend_size_bar),")"),sep="\n")

      line_bar=paste(code_library_text,code_input_text,code_meta_text,code_bar_text,sep="\n")
      writeLines(line_bar, file)
    }
  )

  #filter
  output$code_filter_Download <- downloadHandler(
    filename = "Data_filter_code.txt",
    content = function(file) {
      code_library_text=paste0("library(plotmicrobiome)")


      if(input$data_type == "Amplicon with taxonomic structure"){
        if(input$sep_16s=="tab"){
          sep_char_16s="\\t"
        }else{
          sep_char_16s=input$sep_16s
        }
        code_input_text=paste0("tab1 = format_asv(taxa_file = \"",input$file_16s$name,"\", biom = ",input$biom,", ASV = ",input$ASV,", sep = \"",sep_char_16s,"\", reads_cutoff = ",input$n_reads_16s,", normalization = ",input$norm_16s,", rarefy = ",input$rarefy_16s,", rarefy_num = ",input$rarefy_reads_16s,")")
      }else if(input$data_type == "Metagenomics with taxonomic structure"){
        if(input$sep_wgs=="tab"){
          sep_char_wgs="\\t"
        }else{
          sep_char_wgs=input$sep_wgs
        }
        code_input_text=paste0("tab1 = format_wgs(taxa_file = \"",input$file_wgs$name,"\", sep = \"",sep_char_wgs,"\", reads_cutoff = ",input$n_reads_wgs,", method = \"",input$method_wgs,"\", normalization = ",input$norm_wgs,", rarefy = ",input$rarefy_wgs,", rarefy_num = ",input$rarefy_reads_wgs,")")
      }else{
        if(input$sep_tab=="tab"){
          sep_char_tab="\\t"
        }else{
          sep_char_tab=input$sep_tab
        }
        code_input_text=paste0("tab1 = format_tabs(taxa_file = \"",input$file_tab$name,"\", sep = \"",sep_char_tab,"\", reads_cutoff = ",input$n_reads_tab,", normalization = ",input$norm_tab,", rarefy = ",input$rarefy_tab,", rarefy_num = ",input$rarefy_reads_tab,")")
      }

      if(input$sep_meta=="tab"){
        sep_char_meta="\\t"
      }else{
        sep_char_meta=input$sep_meta
      }
      code_meta_text=paste0("metadata1 = meta_format(metadata = \"",input$file_meta$name,"\", metadata_sep = \"",sep_char_meta,"\", meta_sample_name_col = ",input$meta_sample_name_col,")")

      #table_subset(taxa_table = data_raw(),metadata=data_meta(),stratify_by_metadata=input$stratify_by_metadata_filter,stratify_by_value=input$stratify_by_value_filter,prevalence_cutoff=as.numeric(input$prevalence_cutoff), abundance_cutoff=as.numeric(input$abundance_cutoff),one_level=as.logical(one_level_all()),exclude_ASV_strain=as.logical(input$exclude_ASV_filter))
      code_filter_text=paste0("data_filtered = table_subset(taxa_table = tab1, metadata = metadata1, stratify_by_metadata = \"",input$stratify_by_metadata_filter,"\", stratify_by_value = \"",input$stratify_by_value_filter,"\", one_level = ",as.logical(one_level_all()),", prevalence_cutoff = ",as.numeric(input$prevalence_cutoff),", abundance_cutoff = ",as.numeric(input$abundance_cutoff),", exclude_ASV_strain = ",as.logical(input$exclude_ASV_filter),")")

      line_filter=paste(code_library_text,code_input_text,code_meta_text,code_filter_text,sep="\n")
      writeLines(line_filter, file)
    }
  )

  #stats
  output$code_stat_Download <- downloadHandler(
    filename = "Statistical_test_code.txt",
    content = function(file) {
      code_library_text=paste0("library(plotmicrobiome)")


      if(input$data_type == "Amplicon with taxonomic structure"){
        if(input$sep_16s=="tab"){
          sep_char_16s="\\t"
        }else{
          sep_char_16s=input$sep_16s
        }
        code_input_text=paste0("tab1 = format_asv(taxa_file = \"",input$file_16s$name,"\", biom = ",input$biom,", ASV = ",input$ASV,", sep = \"",sep_char_16s,"\", reads_cutoff = ",input$n_reads_16s,", normalization = ",input$norm_16s,", rarefy = ",input$rarefy_16s,", rarefy_num = ",input$rarefy_reads_16s,")")
      }else if(input$data_type == "Metagenomics with taxonomic structure"){
        if(input$sep_wgs=="tab"){
          sep_char_wgs="\\t"
        }else{
          sep_char_wgs=input$sep_wgs
        }
        code_input_text=paste0("tab1 = format_wgs(taxa_file = \"",input$file_wgs$name,"\", sep = \"",sep_char_wgs,"\", reads_cutoff = ",input$n_reads_wgs,", method = \"",input$method_wgs,"\", normalization = ",input$norm_wgs,", rarefy = ",input$rarefy_wgs,", rarefy_num = ",input$rarefy_reads_wgs,")")
      }else{
        if(input$sep_tab=="tab"){
          sep_char_tab="\\t"
        }else{
          sep_char_tab=input$sep_tab
        }
        code_input_text=paste0("tab1 = format_tabs(taxa_file = \"",input$file_tab$name,"\", sep = \"",sep_char_tab,"\", reads_cutoff = ",input$n_reads_tab,", normalization = ",input$norm_tab,", rarefy = ",input$rarefy_tab,", rarefy_num = ",input$rarefy_reads_tab,")")
      }

      if(input$sep_meta=="tab"){
        sep_char_meta="\\t"
      }else{
        sep_char_meta=input$sep_meta
      }
      code_meta_text=paste0("metadata1 = meta_format(metadata = \"",input$file_meta$name,"\", metadata_sep = \"",sep_char_meta,"\", meta_sample_name_col = ",input$meta_sample_name_col,")")

      #dataf2=table_subset(taxa_table = data_raw(),metadata=data_meta(),stratify_by_metadata=input$stratify_by_metadata_bar,stratify_by_value=input$stratify_by_value_bar,one_level=as.logical(one_level_all()))
      #stat_test(taxa_table =data_filtered(),metadata=data_meta(),test_metadata=input$test_metadata_stat,method=input$method_stat,log_norm=as.logical(input$log_norm_stat),outcome_meta=as.logical(input$outcome_meta),test_metadata_continuous=as.logical(input$test_metadata_continuous_stat),glm_anova=as.logical(input$glm_anova),model_glm=input$model_glm,glm_dist=input$glm_dist,random_effect_var=input$random_effect_var)

      code_stat_text_filter=paste0("data_filtered = table_subset(taxa_table = tab1, metadata = metadata1, stratify_by_metadata = \"",input$stratify_by_metadata_filter,"\", stratify_by_value = \"",input$stratify_by_value_filter,"\", one_level = ",as.logical(one_level_all()),", prevalence_cutoff = ",as.numeric(input$prevalence_cutoff),", abundance_cutoff = ",as.numeric(input$abundance_cutoff),", exclude_ASV_strain = ",as.logical(input$exclude_ASV_filter),")")


      if (input$method_stat %in% c("wilcoxon","t.test","kruskal-wallis","anova","pearson","spearman","kendall")){
        code_stat_text=paste(code_stat_text_filter,
                              paste0("stat_results=stat_test(taxa_table = data_filtered, metadata = metadata1, test_metadata = \"",input$test_metadata_stat,"\", method = \"",input$method_stat,"\", log_norm = ",as.logical(input$log_norm_stat),
                                     ", test_metadata_continuous = ",as.logical(input$test_metadata_continuous_stat),")"),sep="\n")
      }else if (input$method_stat %in% c("glm","lr")){
        code_stat_text=paste(code_stat_text_filter,
                              paste0("stat_results=stat_test(taxa_table = data_filtered, metadata = metadata1, test_metadata = \"",input$test_metadata_stat,"\", method = \"",input$method_stat,"\", log_norm = ",as.logical(input$log_norm_stat),", model_glm = \"",input$model_glm,"\", glm_dist = \"",input$glm_dist,
                                     "\", outcome_meta = ",as.logical(input$outcome_meta),", test_metadata_continuous = ",as.logical(input$test_metadata_continuous_stat),", glm_anova = ",as.logical(input$glm_anova),")"),sep="\n")

      }else if (input$method_stat %in% c("lme")){
        code_stat_text=paste(code_stat_text_filter,
                              paste0("stat_results=stat_test(taxa_table = data_filtered, metadata = metadata1, test_metadata = \"",input$test_metadata_stat,"\", method = \"",input$method_stat,"\", log_norm = ",as.logical(input$log_norm_stat),", model_glm = \"",input$model_glm,"\", glm_dist = \"",input$glm_dist,"\", random_effect_var = \"",input$random_effect_var,
                                     "\", outcome_meta = ",as.logical(input$outcome_meta),", test_metadata_continuous = ",as.logical(input$test_metadata_continuous_stat),", glm_anova = ",as.logical(input$glm_anova),")"),sep="\n")
      }

      line_stat=paste(code_library_text,code_input_text,code_meta_text,code_stat_text,sep="\n")
      writeLines(line_stat, file)
    }
  )

  #tree
  output$code_tree_Download <- downloadHandler(
    filename = "treeplot_code.txt",
    content = function(file) {
      code_library_text=paste0("library(plotmicrobiome)")


      if(input$data_type == "Amplicon with taxonomic structure"){
        if(input$sep_16s=="tab"){
          sep_char_16s="\\t"
        }else{
          sep_char_16s=input$sep_16s
        }
        code_input_text=paste0("tab1 = format_asv(taxa_file = \"",input$file_16s$name,"\", biom = ",input$biom,", ASV = ",input$ASV,", sep = \"",sep_char_16s,"\", reads_cutoff = ",input$n_reads_16s,", normalization = ",input$norm_16s,", rarefy = ",input$rarefy_16s,", rarefy_num = ",input$rarefy_reads_16s,")")
      }else if(input$data_type == "Metagenomics with taxonomic structure"){
        if(input$sep_wgs=="tab"){
          sep_char_wgs="\\t"
        }else{
          sep_char_wgs=input$sep_wgs
        }
        code_input_text=paste0("tab1 = format_wgs(taxa_file = \"",input$file_wgs$name,"\", sep = \"",sep_char_wgs,"\", reads_cutoff = ",input$n_reads_wgs,", method = \"",input$method_wgs,"\", normalization = ",input$norm_wgs,", rarefy = ",input$rarefy_wgs,", rarefy_num = ",input$rarefy_reads_wgs,")")
      }else{
        if(input$sep_tab=="tab"){
          sep_char_tab="\\t"
        }else{
          sep_char_tab=input$sep_tab
        }
        code_input_text=paste0("tab1 = format_tabs(taxa_file = \"",input$file_tab$name,"\", sep = \"",sep_char_tab,"\", reads_cutoff = ",input$n_reads_tab,", normalization = ",input$norm_tab,", rarefy = ",input$rarefy_tab,", rarefy_num = ",input$rarefy_reads_tab,")")
      }

      if(input$sep_meta=="tab"){
        sep_char_meta="\\t"
      }else{
        sep_char_meta=input$sep_meta
      }
      code_meta_text=paste0("metadata1 = meta_format(metadata = \"",input$file_meta$name,"\", metadata_sep = \"",sep_char_meta,"\", meta_sample_name_col = ",input$meta_sample_name_col,")")

      #dataf2=table_subset(taxa_table = data_raw(),metadata=data_meta(),stratify_by_metadata=input$stratify_by_metadata_bar,stratify_by_value=input$stratify_by_value_bar,one_level=as.logical(one_level_all()))
      #stat_test(taxa_table =data_filtered(),metadata=data_meta(),test_metadata=input$test_metadata_stat,method=input$method_stat,log_norm=as.logical(input$log_norm_stat),outcome_meta=as.logical(input$outcome_meta),test_metadata_continuous=as.logical(input$test_metadata_continuous_stat),glm_anova=as.logical(input$glm_anova),model_glm=input$model_glm,glm_dist=input$glm_dist,random_effect_var=input$random_effect_var)

      code_stat_text_filter=paste0("data_filtered = table_subset(taxa_table = tab1, metadata = metadata1, stratify_by_metadata = \"",input$stratify_by_metadata_filter,"\", stratify_by_value = \"",input$stratify_by_value_filter,"\", one_level = ",as.logical(one_level_all()),", prevalence_cutoff = ",as.numeric(input$prevalence_cutoff),", abundance_cutoff = ",as.numeric(input$abundance_cutoff),", exclude_ASV_strain = ",as.logical(input$exclude_ASV_filter),")")


      if (input$method_stat %in% c("wilcoxon","t.test","kruskal-wallis","anova","pearson","spearman","kendall")){
        code_stat_text=paste(code_stat_text_filter,
                             paste0("stat_results=stat_test(taxa_table = data_filtered, metadata = metadata1, test_metadata = \"",input$test_metadata_stat,"\", method = \"",input$method_stat,"\", log_norm = ",as.logical(input$log_norm_stat),
                                    ", test_metadata_continuous = ",as.logical(input$test_metadata_continuous_stat),")"),sep="\n")
      }else if (input$method_stat %in% c("glm","lr")){
        code_stat_text=paste(code_stat_text_filter,
                             paste0("stat_results=stat_test(taxa_table = data_filtered, metadata = metadata1, test_metadata = \"",input$test_metadata_stat,"\", method = \"",input$method_stat,"\", log_norm = ",as.logical(input$log_norm_stat),", model_glm = \"",input$model_glm,"\", glm_dist = \"",input$glm_dist,
                                    "\", outcome_meta = ",as.logical(input$outcome_meta),", test_metadata_continuous = ",as.logical(input$test_metadata_continuous_stat),", glm_anova = ",as.logical(input$glm_anova),")"),sep="\n")

      }else if (input$method_stat %in% c("lme")){
        code_stat_text=paste(code_stat_text_filter,
                             paste0("stat_results=stat_test(taxa_table = data_filtered, metadata = metadata1, test_metadata = \"",input$test_metadata_stat,"\", method = \"",input$method_stat,"\", log_norm = ",as.logical(input$log_norm_stat),", model_glm = \"",input$model_glm,"\", glm_dist = \"",input$glm_dist,"\", random_effect_var = \"",input$random_effect_var,
                                    "\", outcome_meta = ",as.logical(input$outcome_meta),", test_metadata_continuous = ",as.logical(input$test_metadata_continuous_stat),", glm_anova = ",as.logical(input$glm_anova),")"),sep="\n")
      }

      #tree_view(taxa_table = data_raw(),metadata=data_meta(),fdrs=data_fdrs(),test_metadata=input$test_metadata_stat,fdr_cutoff=input$fdr_cutoff_tree,test_metadata_continuous=as.logical(input$test_metadata_continuous_tree),
      #node_size_breaks=as.numeric(strsplit(input$node_size_breaks, ",\\s*")[[1]]),palette_highlight=strsplit(input$palette_group_tree, ",\\s*")[[1]],single_parent_branch_removal=as.logical(input$single_parent_branch_removal),single_child_branch_removal=as.logical(input$single_child_branch_removal),
      #prevalence_cutoff=as.numeric(input$prevalence_cutoff_tree), abundance_cutoff=as.numeric(input$abundance_cutoff_tree),taxa_removal=input$taxa_removal_tree,branch=input$branch_tree)

      code_tree_text=paste0("tree_view(taxa_table = data_filtered, metadata = metadata1, fdrs = stat_results, test_metadata =\"",input$test_metadata_stat,"\", fdr_cutoff = ",as.numeric(input$fdr_cutoff_tree),", test_metadata_continuous = ",as.logical(input$test_metadata_continuous_tree),
                            ", node_size_breaks = c(",paste(strsplit(input$node_size_breaks, ",\\s*")[[1]],collapse = ","),"), palette_highlight = c(\"",paste(strsplit(input$palette_group_tree, ",\\s*")[[1]],collapse = "\", \""),"\"), single_parent_branch_removal = ",as.logical(input$single_parent_branch_removal),", single_child_branch_removal = ",as.logical(input$single_child_branch_removal),
                            ", prevalence_cutoff = ",as.numeric(input$prevalence_cutoff_tree), ", abundance_cutoff = ",as.numeric(input$abundance_cutoff_tree),", taxa_removal = \"",input$taxa_removal_tree,"\", branch = \"",input$branch_tree,"\")")


      line_stat=paste(code_library_text,code_input_text,code_meta_text,code_stat_text,code_tree_text,sep="\n")
      writeLines(line_stat, file)
    }
  )

  #boxplot
  output$code_box_Download <- downloadHandler(
    filename = "boxplot_code.txt",
    content = function(file) {
      code_library_text=paste0("library(plotmicrobiome)")


      if(input$data_type == "Amplicon with taxonomic structure"){
        if(input$sep_16s=="tab"){
          sep_char_16s="\\t"
        }else{
          sep_char_16s=input$sep_16s
        }
        code_input_text=paste0("tab1 = format_asv(taxa_file = \"",input$file_16s$name,"\", biom = ",input$biom,", ASV = ",input$ASV,", sep = \"",sep_char_16s,"\", reads_cutoff = ",input$n_reads_16s,", normalization = ",input$norm_16s,", rarefy = ",input$rarefy_16s,", rarefy_num = ",input$rarefy_reads_16s,")")
      }else if(input$data_type == "Metagenomics with taxonomic structure"){
        if(input$sep_wgs=="tab"){
          sep_char_wgs="\\t"
        }else{
          sep_char_wgs=input$sep_wgs
        }
        code_input_text=paste0("tab1 = format_wgs(taxa_file = \"",input$file_wgs$name,"\", sep = \"",sep_char_wgs,"\", reads_cutoff = ",input$n_reads_wgs,", method = \"",input$method_wgs,"\", normalization = ",input$norm_wgs,", rarefy = ",input$rarefy_wgs,", rarefy_num = ",input$rarefy_reads_wgs,")")
      }else{
        if(input$sep_tab=="tab"){
          sep_char_tab="\\t"
        }else{
          sep_char_tab=input$sep_tab
        }
        code_input_text=paste0("tab1 = format_tabs(taxa_file = \"",input$file_tab$name,"\", sep = \"",sep_char_tab,"\", reads_cutoff = ",input$n_reads_tab,", normalization = ",input$norm_tab,", rarefy = ",input$rarefy_tab,", rarefy_num = ",input$rarefy_reads_tab,")")
      }

      if(input$sep_meta=="tab"){
        sep_char_meta="\\t"
      }else{
        sep_char_meta=input$sep_meta
      }
      code_meta_text=paste0("metadata1 = meta_format(metadata = \"",input$file_meta$name,"\", metadata_sep = \"",sep_char_meta,"\", meta_sample_name_col = ",input$meta_sample_name_col,")")

      #dataf2=table_subset(taxa_table = data_raw(),metadata=data_meta(),stratify_by_metadata=input$stratify_by_metadata_bar,stratify_by_value=input$stratify_by_value_bar,one_level=as.logical(one_level_all()))
      #stat_test(taxa_table =data_filtered(),metadata=data_meta(),test_metadata=input$test_metadata_stat,method=input$method_stat,log_norm=as.logical(input$log_norm_stat),outcome_meta=as.logical(input$outcome_meta),test_metadata_continuous=as.logical(input$test_metadata_continuous_stat),glm_anova=as.logical(input$glm_anova),model_glm=input$model_glm,glm_dist=input$glm_dist,random_effect_var=input$random_effect_var)

      code_stat_text_filter=paste0("data_filtered = table_subset(taxa_table = tab1, metadata = metadata1, stratify_by_metadata = \"",input$stratify_by_metadata_filter,"\", stratify_by_value = \"",input$stratify_by_value_filter,"\", one_level = ",as.logical(one_level_all()),", prevalence_cutoff = ",as.numeric(input$prevalence_cutoff),", abundance_cutoff = ",as.numeric(input$abundance_cutoff),", exclude_ASV_strain = ",as.logical(input$exclude_ASV_filter),")")


      if (input$method_stat %in% c("wilcoxon","t.test","kruskal-wallis","anova","pearson","spearman","kendall")){
        code_stat_text=paste(code_stat_text_filter,
                             paste0("stat_results=stat_test(taxa_table = data_filtered, metadata = metadata1, test_metadata = \"",input$test_metadata_stat,"\", method = \"",input$method_stat,"\", log_norm = ",as.logical(input$log_norm_stat),
                                    ", test_metadata_continuous = ",as.logical(input$test_metadata_continuous_stat),")"),sep="\n")
      }else if (input$method_stat %in% c("glm","lr")){
        code_stat_text=paste(code_stat_text_filter,
                             paste0("stat_results=stat_test(taxa_table = data_filtered, metadata = metadata1, test_metadata = \"",input$test_metadata_stat,"\", method = \"",input$method_stat,"\", log_norm = ",as.logical(input$log_norm_stat),", model_glm = \"",input$model_glm,"\", glm_dist = \"",input$glm_dist,
                                    "\", outcome_meta = ",as.logical(input$outcome_meta),", test_metadata_continuous = ",as.logical(input$test_metadata_continuous_stat),", glm_anova = ",as.logical(input$glm_anova),")"),sep="\n")

      }else if (input$method_stat %in% c("lme")){
        code_stat_text=paste(code_stat_text_filter,
                             paste0("stat_results=stat_test(taxa_table = data_filtered, metadata = metadata1, test_metadata = \"",input$test_metadata_stat,"\", method = \"",input$method_stat,"\", log_norm = ",as.logical(input$log_norm_stat),", model_glm = \"",input$model_glm,"\", glm_dist = \"",input$glm_dist,"\", random_effect_var = \"",input$random_effect_var,
                                    "\", outcome_meta = ",as.logical(input$outcome_meta),", test_metadata_continuous = ",as.logical(input$test_metadata_continuous_stat),", glm_anova = ",as.logical(input$glm_anova),")"),sep="\n")
      }

      #taxa_boxplot(taxa_table =data_filtered(),metadata=data_meta(),fdrs=data_fdrs(),test_metadata=input$test_metadata_stat,test_metadata_order=input$test_metadata_order_box,one_level=as.logical(one_level_all()),
      #log_norm=input$log_norm_box,taxa_shown=input$taxa_shown_box,cutoff=input$fdr_cutoff_box,page=input$page_box,xlab=input$xlab_box,ylab=input$ylab_box,
      #xlab_direction=as.integer(input$x_dir_box),palette_group=strsplit(input$palette_group_box, ",\\s*")[[1]])
      #show each page
      code_box_text1=paste0("taxa_boxplot(taxa_table = data_filtered, metadata = metadata1, fdrs = stat_results, test_metadata =\"",input$test_metadata_stat,"\", cutoff = ",as.numeric(input$fdr_cutoff_box),", test_metadata_order = c(\"",paste(input$test_metadata_order_box,collapse = "\",\""),
                            "\"), one_level = ",as.logical(one_level_all()),", log_norm = ",as.logical(input$log_norm_box),", taxa_shown = \"",input$taxa_shown_box, "\", palette_group = c(\"",paste(strsplit(input$palette_group_box, ",\\s*")[[1]],collapse = "\", \""),"\"), xlab = \"",input$xlab_box,"\", ylab = \"",input$xlab_box,
                            "\", xlab_direction = ",as.integer(input$x_dir_box),", page = ",as.numeric(input$page_box), ")")
      #download all pages
      code_box_text2=paste0("taxa_boxplot_download(taxa_table = data_filtered, metadata = metadata1, fdrs = stat_results, test_metadata =\"",input$test_metadata_stat,"\", cutoff = ",as.numeric(input$fdr_cutoff_box),", test_metadata_order = c(\"",paste(input$test_metadata_order_box,collapse = "\",\""),
                            "\"), one_level = ",as.logical(one_level_all()),", log_norm = ",as.logical(input$log_norm_box),", taxa_shown = \"",input$taxa_shown_box, "\", palette_group = c(\"",paste(strsplit(input$palette_group_box, ",\\s*")[[1]],collapse = "\", \""),"\"), xlab = \"",input$xlab_box,"\", ylab = \"",input$xlab_box,
                            "\", xlab_direction = ",as.integer(input$x_dir_box), ")")


      line_box=paste(code_library_text,code_input_text,code_meta_text,code_stat_text,"#to show individual page",code_box_text1,"#to download all pages","pdf(\"boxplots.pdf\",onefile=True)",code_box_text2,"dev.off()",sep="\n")
      writeLines(line_box, file)
    }
  )

  #correlation plot
  output$code_cor_Download <- downloadHandler(
    filename = "Correlation_plot_code.txt",
    content = function(file) {
      code_library_text=paste0("library(plotmicrobiome)")


      if(input$data_type == "Amplicon with taxonomic structure"){
        if(input$sep_16s=="tab"){
          sep_char_16s="\\t"
        }else{
          sep_char_16s=input$sep_16s
        }
        code_input_text=paste0("tab1 = format_asv(taxa_file = \"",input$file_16s$name,"\", biom = ",input$biom,", ASV = ",input$ASV,", sep = \"",sep_char_16s,"\", reads_cutoff = ",input$n_reads_16s,", normalization = ",input$norm_16s,", rarefy = ",input$rarefy_16s,", rarefy_num = ",input$rarefy_reads_16s,")")
      }else if(input$data_type == "Metagenomics with taxonomic structure"){
        if(input$sep_wgs=="tab"){
          sep_char_wgs="\\t"
        }else{
          sep_char_wgs=input$sep_wgs
        }
        code_input_text=paste0("tab1 = format_wgs(taxa_file = \"",input$file_wgs$name,"\", sep = \"",sep_char_wgs,"\", reads_cutoff = ",input$n_reads_wgs,", method = \"",input$method_wgs,"\", normalization = ",input$norm_wgs,", rarefy = ",input$rarefy_wgs,", rarefy_num = ",input$rarefy_reads_wgs,")")
      }else{
        if(input$sep_tab=="tab"){
          sep_char_tab="\\t"
        }else{
          sep_char_tab=input$sep_tab
        }
        code_input_text=paste0("tab1 = format_tabs(taxa_file = \"",input$file_tab$name,"\", sep = \"",sep_char_tab,"\", reads_cutoff = ",input$n_reads_tab,", normalization = ",input$norm_tab,", rarefy = ",input$rarefy_tab,", rarefy_num = ",input$rarefy_reads_tab,")")
      }

      if(input$sep_meta=="tab"){
        sep_char_meta="\\t"
      }else{
        sep_char_meta=input$sep_meta
      }
      code_meta_text=paste0("metadata1 = meta_format(metadata = \"",input$file_meta$name,"\", metadata_sep = \"",sep_char_meta,"\", meta_sample_name_col = ",input$meta_sample_name_col,")")

      #dataf2=table_subset(taxa_table = data_raw(),metadata=data_meta(),stratify_by_metadata=input$stratify_by_metadata_bar,stratify_by_value=input$stratify_by_value_bar,one_level=as.logical(one_level_all()))
      #stat_test(taxa_table =data_filtered(),metadata=data_meta(),test_metadata=input$test_metadata_stat,method=input$method_stat,log_norm=as.logical(input$log_norm_stat),outcome_meta=as.logical(input$outcome_meta),test_metadata_continuous=as.logical(input$test_metadata_continuous_stat),glm_anova=as.logical(input$glm_anova),model_glm=input$model_glm,glm_dist=input$glm_dist,random_effect_var=input$random_effect_var)

      code_stat_text_filter=paste0("data_filtered = table_subset(taxa_table = tab1, metadata = metadata1, stratify_by_metadata = \"",input$stratify_by_metadata_filter,"\", stratify_by_value = \"",input$stratify_by_value_filter,"\", one_level = ",as.logical(one_level_all()),", prevalence_cutoff = ",as.numeric(input$prevalence_cutoff),", abundance_cutoff = ",as.numeric(input$abundance_cutoff),", exclude_ASV_strain = ",as.logical(input$exclude_ASV_filter),")")


      if (input$method_stat %in% c("wilcoxon","t.test","kruskal-wallis","anova","pearson","spearman","kendall")){
        code_stat_text=paste(code_stat_text_filter,
                             paste0("stat_results=stat_test(taxa_table = data_filtered, metadata = metadata1, test_metadata = \"",input$test_metadata_stat,"\", method = \"",input$method_stat,"\", log_norm = ",as.logical(input$log_norm_stat),
                                    ", test_metadata_continuous = ",as.logical(input$test_metadata_continuous_stat),")"),sep="\n")
      }else if (input$method_stat %in% c("glm","lr")){
        code_stat_text=paste(code_stat_text_filter,
                             paste0("stat_results=stat_test(taxa_table = data_filtered, metadata = metadata1, test_metadata = \"",input$test_metadata_stat,"\", method = \"",input$method_stat,"\", log_norm = ",as.logical(input$log_norm_stat),", model_glm = \"",input$model_glm,"\", glm_dist = \"",input$glm_dist,
                                    "\", outcome_meta = ",as.logical(input$outcome_meta),", test_metadata_continuous = ",as.logical(input$test_metadata_continuous_stat),", glm_anova = ",as.logical(input$glm_anova),")"),sep="\n")

      }else if (input$method_stat %in% c("lme")){
        code_stat_text=paste(code_stat_text_filter,
                             paste0("stat_results=stat_test(taxa_table = data_filtered, metadata = metadata1, test_metadata = \"",input$test_metadata_stat,"\", method = \"",input$method_stat,"\", log_norm = ",as.logical(input$log_norm_stat),", model_glm = \"",input$model_glm,"\", glm_dist = \"",input$glm_dist,"\", random_effect_var = \"",input$random_effect_var,
                                    "\", outcome_meta = ",as.logical(input$outcome_meta),", test_metadata_continuous = ",as.logical(input$test_metadata_continuous_stat),", glm_anova = ",as.logical(input$glm_anova),")"),sep="\n")
      }

      #meta_corplot(taxa_table =data_filtered(),metadata=data_meta(),test_metadata=input$test_metadata_cor,cor_method=input$cor_method,one_level=as.logical(one_level_all()),
      #col_metadata=input$col_metadata_cor,log_norm=input$log_norm_cor,taxa_shown=input$taxa_shown_cor,page=input$page_cor,xlab=input$xlab_cor,ylab=input$ylab_cor,
      #fdr_cutoff=input$fdr_cutoff_cor,palette_group=strsplit(input$palette_group_cor, ",\\s*")[[1]])

      #show each page
      code_cor_text1=paste0("meta_corplot(taxa_table = data_filtered, metadata = metadata1,  test_metadata =\"",input$test_metadata_cor,"\", fdr_cutoff = ",as.numeric(input$fdr_cutoff_cor),", cor_method = \"",input$cor_method,
                            "\", one_level = ",as.logical(one_level_all()),", col_metadata = \"",input$col_metadata_cor,"\", log_norm = ",as.logical(input$log_norm_cor),", taxa_shown = \"",input$taxa_shown_cor, "\", palette_group = c(\"",paste(strsplit(input$palette_group_cor, ",\\s*")[[1]],collapse = "\", \""),"\"), ylab = \"",input$ylab_cor,
                            "\", page = ",as.numeric(input$page_cor), ")")
      #download all pages
      code_cor_text2=paste0("meta_corplot_download(taxa_table = data_filtered, metadata = metadata1,  test_metadata =\"",input$test_metadata_cor,"\", fdr_cutoff = ",as.numeric(input$fdr_cutoff_cor),", cor_method = \"",input$cor_method,
                            "\", one_level = ",as.logical(one_level_all()),", col_metadata = \"",input$col_metadata_cor,"\", log_norm = ",as.logical(input$log_norm_cor),", taxa_shown = \"",input$taxa_shown_cor, "\", palette_group = c(\"",paste(strsplit(input$palette_group_cor, ",\\s*")[[1]],collapse = "\", \""),"\"), ylab = \"",input$ylab_cor, "\")")


      line_cor=paste(code_library_text,code_input_text,code_meta_text,code_stat_text,"#to show individual page",code_cor_text1,"#to download all pages","pdf(\"correlation_plots.pdf\",onefile=True)",code_cor_text2,"dev.off()",sep="\n")
      writeLines(line_cor, file)
    }
  )

  #P vs P plot
  output$code_pvp_Download <- downloadHandler(
    filename = "Pvalues_compare_code.txt",
    content = function(file) {
      code_library_text=paste0("library(plotmicrobiome)")


      if(input$data_type == "Amplicon with taxonomic structure"){
        if(input$sep_16s=="tab"){
          sep_char_16s="\\t"
        }else{
          sep_char_16s=input$sep_16s
        }
        code_input_text=paste0("tab1 = format_asv(taxa_file = \"",input$file_16s$name,"\", biom = ",input$biom,", ASV = ",input$ASV,", sep = \"",sep_char_16s,"\", reads_cutoff = ",input$n_reads_16s,", normalization = ",input$norm_16s,", rarefy = ",input$rarefy_16s,", rarefy_num = ",input$rarefy_reads_16s,")")
      }else if(input$data_type == "Metagenomics with taxonomic structure"){
        if(input$sep_wgs=="tab"){
          sep_char_wgs="\\t"
        }else{
          sep_char_wgs=input$sep_wgs
        }
        code_input_text=paste0("tab1 = format_wgs(taxa_file = \"",input$file_wgs$name,"\", sep = \"",sep_char_wgs,"\", reads_cutoff = ",input$n_reads_wgs,", method = \"",input$method_wgs,"\", normalization = ",input$norm_wgs,", rarefy = ",input$rarefy_wgs,", rarefy_num = ",input$rarefy_reads_wgs,")")
      }else{
        if(input$sep_tab=="tab"){
          sep_char_tab="\\t"
        }else{
          sep_char_tab=input$sep_tab
        }
        code_input_text=paste0("tab1 = format_tabs(taxa_file = \"",input$file_tab$name,"\", sep = \"",sep_char_tab,"\", reads_cutoff = ",input$n_reads_tab,", normalization = ",input$norm_tab,", rarefy = ",input$rarefy_tab,", rarefy_num = ",input$rarefy_reads_tab,")")
      }

      if(input$sep_meta=="tab"){
        sep_char_meta="\\t"
      }else{
        sep_char_meta=input$sep_meta
      }
      code_meta_text=paste0("metadata1 = meta_format(metadata = \"",input$file_meta$name,"\", metadata_sep = \"",sep_char_meta,"\", meta_sample_name_col = ",input$meta_sample_name_col,")")


      if(input$upload_p1=="Yes"){
        if(input$sep_p1=="tab"){
          sep_p1="\t"
        }else{
          sep_p1=input$sep_p1
        }
        code_p1_text=paste0("p1 = read.table(file = \"",input$file_p1$name,"\", sep = ",sep_p1,", row.names=1,header=T,check.names = F)")
      }else{
        #dataf1_p1=table_subset(taxa_table = data_raw(),metadata=data_meta(),stratify_by_metadata=input$stratify_by_metadata_p1,stratify_by_value=input$stratify_by_value_p1,prevalence_cutoff=as.numeric(input$prevalence_cutoff_p1), abundance_cutoff=as.numeric(input$abundance_cutoff_p1),one_level=as.logical(one_level_all()),exclude_ASV_strain=as.logical(input$exclude_ASV_p1))
        #stat_test(taxa_table =dataf1_p1,metadata=data_meta(),test_metadata=input$test_metadata_p1,method=input$method_stat_p1,test_metadata_continuous=as.logical(input$test_metadata_continuous_p1))
        code_p1_text1=paste0("dataf1_p1 = table_subset(taxa_table = tab1, metadata = metadata1, stratify_by_metadata = \"",input$stratify_by_metadata_p1,"\", stratify_by_value = \"",input$stratify_by_value_p1,"\", one_level = ",as.logical(one_level_all()),", prevalence_cutoff = ",as.numeric(input$prevalence_cutoff_p1),", abundance_cutoff = ",as.numeric(input$abundance_cutoff_p1),", exclude_ASV_strain = ",as.logical(input$exclude_ASV_p1),")")
        code_p1_text2=paste0("p1 = stat_test(taxa_table = dataf1_p1, metadata = metadata1, test_metadata = \"",input$test_metadata_p1,"\", method = \"",input$method_stat_p1,"\",test_metadata_continuous = ",as.logical(input$test_metadata_continuous_p1),")")
        code_p1_text=paste(code_p1_text1,code_p1_text2,sep="\n")
      }

      if(input$upload_p2=="Yes"){
        if(input$sep_p2=="tab"){
          sep_p2="\t"
        }else{
          sep_p2=input$sep_p2
        }
        code_p2_text=paste0("p2 = read.table(file = \"",input$file_p1$name,"\", sep = ",sep_p2,", row.names=1,header=T,check.names = F)")
      }else{
        #dataf1_p1=table_subset(taxa_table = data_raw(),metadata=data_meta(),stratify_by_metadata=input$stratify_by_metadata_p1,stratify_by_value=input$stratify_by_value_p1,prevalence_cutoff=as.numeric(input$prevalence_cutoff_p1), abundance_cutoff=as.numeric(input$abundance_cutoff_p1),one_level=as.logical(one_level_all()),exclude_ASV_strain=as.logical(input$exclude_ASV_p1))
        #stat_test(taxa_table =dataf1_p1,metadata=data_meta(),test_metadata=input$test_metadata_p1,method=input$method_stat_p1,test_metadata_continuous=as.logical(input$test_metadata_continuous_p1))
        code_p2_text1=paste0("dataf1_p2 = table_subset(taxa_table = tab1, metadata = metadata1, stratify_by_metadata = \"",input$stratify_by_metadata_p2,"\", stratify_by_value = \"",input$stratify_by_value_p2,"\", one_level = ",as.logical(one_level_all()),", prevalence_cutoff = ",as.numeric(input$prevalence_cutoff_p2),", abundance_cutoff = ",as.numeric(input$abundance_cutoff_p2),", exclude_ASV_strain = ",as.logical(input$exclude_ASV_p2),")")
        code_p2_text2=paste0("p2 = stat_test(taxa_table = dataf1_p2, metadata = metadata1, test_metadata = \"",input$test_metadata_p2,"\", method = \"",input$method_stat_p2,"\",test_metadata_continuous = ",as.logical(input$test_metadata_continuous_p2),")")
        code_p2_text=paste(code_p2_text1,code_p2_text2,sep="\n")
      }

      #p_compare(data_p1(),data_p2(),p_col1=as.numeric(input$p1_col),p_col2=as.numeric(input$p2_col),indicator1=as.numeric(input$ind1_col),indicator2=as.numeric(input$ind2_col),point_color=input$point_color,lab_cutoff=as.numeric(input$lab_cutoff),cor_method=input$cor_method_p,
      #x.reverse = as.logical(input$x_reverse),y.reverse = as.logical(input$y_reverse),exclude_unclassified=as.logical(input$exclude_unclassified),one_level=as.logical(one_level_all()),direction=as.logical(input$direction_pvp))
      code_pvp_text=paste0("p_compare(p1, p2, p_col1 = ",as.numeric(input$p1_col),", p_col2 = ",as.numeric(input$p2_col),", indicator1 = ",as.numeric(input$ind1_col),", indicator2 = ",as.numeric(input$ind2_col),", point_color = \"",input$point_color, "\", lab_cutoff=",as.numeric(input$lab_cutoff),
                           ", cor_method = \"",input$cor_method_p,"\", x.reverse = ",as.logical(input$x_reverse),", y.reverse = ",as.logical(input$y_reverse),", exclude_unclassified = ",as.logical(input$exclude_unclassified),", one_level = ",as.logical(one_level_all()),", direction = ",as.logical(input$direction_pvp),")")

      line_pvp=paste(code_library_text,code_input_text,code_meta_text,code_p1_text,code_p2_text,code_pvp_text,sep="\n")
      writeLines(line_pvp, file)
    }
  )


}

# Run the application
shinyApp(ui = ui, server = server)


