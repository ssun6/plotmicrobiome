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
  navbarPage(
    "plotmicrobiome",
    id = "main_navbar",

    tabPanel(
      "Data input",
      sidebarLayout(
        sidebarPanel(
          useShinyjs(),
          width = 3,
          br(),
          h4("Data input"),
          selectInput("data_type", "Is this 16S taxonomic data, metagenomics taxonomic data or other one level data?", c("16S","Metagenomics","One level")),
          br(),
          br(),
          div(id = "16S",
            h4("16S taxonomic data"),
            h5("Samples should be in rows and metadata should be in columns for .csv and .tsv files."),
            fileInput("file_16s", "Choose 16S taxa file (accept .csv, .tsv and .biom files.)"),
            selectInput("onefile", "Is the data in one file?", c("True","False")),
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
            br()
          ),
          div(id = "wgs",
            h4("Metagenomics taxonomic data"),
            h5("Samples should be in rows and metadata should be in columns."),
            fileInput("file_wgs", "Choose metagenomics taxa file"),
            selectInput("method_wgs", "Which tool was used for taxonomic classification?", c("kraken2","metaphlan2")),
            selectInput("sep_wgs", "What is the delimiter?", c(",","tab")),
            numericInput("n_reads_wgs", "Exclude samples with the number of reads lower than", value = 1000),
            selectInput("norm_wgs", "Normalize the data? (Supports proportion scaled by average sequencing depth and rarefaction)", c("True","False")),
            selectInput("rarefy_wgs", "Use rarefaction for normalization? Default is proportion scaled by average sequencing depth.", c("False","True")),
            numericInput("rarefy_reads_wgs", "Rarefy samples to how many reads?", value = 1000),
            numericInput("n_raw_wgs", "Preview rows", value = 5, min = 1, step = 1),
            numericInput("n_raw_wgs_col", "Preview columns", value = 10, min = 1, step = 1),
            br(),
            br(),
            br()
          ),
          div(id = "one_level",
            h4("One level data"),
            h5("Samples should be in rows and metadata should be in columns."),
            fileInput("file_path", "Choose file"),
            selectInput("sep_path", "What is the delimiter?", c(",","tab")),
            numericInput("n_reads_path", "Exclude samples with the number of reads lower than", value = 0),
            selectInput("norm_path", "Normalize the data? (Supports proportion scaled by average sequencing depth and rarefaction)", c("False","True")),
            selectInput("rarefy_path", "Use rarefaction for normalization? Default is proportion scaled by average sequencing depth.", c("False","True")),
            numericInput("rarefy_reads_path", "Rarefy samples to how many reads?", value = 1000),
            numericInput("n_raw_path", "Preview rows", value = 5, min = 1, step = 1),
            numericInput("n_raw_path_col", "Preview columns", value = 10, min = 1, step = 1),
            br(),
            br(),
            br()
          ),
          actionButton("button_raw", "Run"),
          textOutput("data_dim"),
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
          h5("Samples should be in rows and metadata should be in columns."),
          fileInput("file_meta", "Choose metadata file"),
          #textInput("dir_meta", "Metadata directory",value="./data-raw/metadata_cafe.csv"),
          selectInput("sep_meta", "What is the delimiter?", c(",","tab")),
          numericInput("meta_sample_name_col", "Which column are the sample names in?", value = 0),
          numericInput("n_meta", "Preview rows", value = 5, min = 1, step = 1),
          actionButton("button_meta", "Run"),
          br(),
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
          selectInput("taxa_level_mds", "Select the taxonomic level shown (Kraken2 does not have strain level)",c("phylum","class","order","family","genus","species","ASV_or_strain")),
          selectInput("method_mds", "Which method should be used for ordination?", c("pcoa","nmds")),
          selectInput("distance_type", "Distance type",c("bray","euclidean","manhattan","jaccard")),
          selectInput("log_normalization_mds", "Should the data be log10 normalized?", c("False","True")),
          textAreaInput("palette_group_mds", "Colors for plot", value = "red,blue,orange,green"),
          actionButton("button_mds", "Run"),
          br(),
          br(),
          br(),
          sliderInput("mds_h", label ="Figure height", min = 0.5, max = 5, value = 1),
          sliderInput("mds_w", label ="Figure width", min = 0.5, max = 5, value = 1),
          br(),
          h5("Download figure:"),
          downloadButton("plotMDSDownload", "Download"),
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
          textInput("test_metadata_order_alpha", "Type in the order of metadata separated with comma to change those in the figure ",value="default"),
          selectInput("method_alpha", "Select statistical test method",c("wilcoxon","t.test","kruskal-wallis","anova")),
          textAreaInput("palette_group_alpha", "Colors for plot", value = "red,blue,orange,green"),
          selectInput("x_dir_alpha", "Direction of X axis labels (1 is horizontal, 2 is vertical)", c(1,2)),
          actionButton("button_alpha", "Run"),
          br(),
          br(),
          sliderInput("alpha_h", label ="Figure height", min = 0.5, max = 5, value = 1),
          sliderInput("alpha_w", label ="Figure width", min = 0.5, max = 5, value = 1),
          br(),
          br(),
          h5("Download figure:"),
          downloadButton("plotAlphaDownload", "Download"),
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
          textInput("test_metadata_order_bar", "Type in the order of metadata separated with comma to change those in the figure ",value="default"),
          textInput("num_taxa_bar", "Select the number of taxa shown",value=8),
          selectInput("taxa_level_bar", "Select the taxonomic level shown",c("phylum","class","order","family","genus")),
          textAreaInput("palette_group_bar", "Colors for plot", value = "default"),
          selectInput("x_dir_bar", "Direction of X axis labels (1 is horizontal, 2 is vertical)",c(1,2)),
          textInput("legend_size_bar", "Select the legend size", value = 1.5),
          actionButton("button_bar", "Run"),
          br(),
          br(),
          sliderInput("bar_h", label ="Figure height", min = 0.5, max = 5, value = 1),
          sliderInput("bar_w", label ="Figure width", min = 0.5, max = 5, value = 1),
          br(),
          br(),
          h5("Download figure:"),
          downloadButton("plotBarDownload", "Download"),
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
          selectInput("exclude_ASV_filter", "Should ASV be excluded from the analysis (this changes the P-value distribution and FDR)?", c("True","False")),
          numericInput("prevalence_cutoff", "Prevalence cutoff", value = 0.25),
          numericInput("abundance_cutoff", "Abundance cutoff", value = 0),
          numericInput("n_filtered", "Preview rows", value = 5, min = 1, step = 1),
          numericInput("n_filtered_col", "Preview columns", value = 10, min = 1, step = 1),
          actionButton("button_filter", "Run"),
          textOutput("filtered_dim"),
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
          selectInput("test_metadata_continuous_stat", "Is the metadata for testing continuous?",c("False","True"),selected ="False"),
          uiOutput(outputId = 'method_statUI'),
          selectInput(inputId = "log_norm_stat", "Should the data be log10 normalization?",c("True","False"),selected="True"),
          h5("Variable preview:"),
          textOutput("head_test_stat"),
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
          br(),
          br(),
          selectInput("test_metadata_continuous_tree", "Is the test metadata continuous?",c("False","True")),
          numericInput("fdr_cutoff_tree", "FDR cutoff", value = 0.1, min = 0, step = 0.01),
          textAreaInput("node_size_breaks", "Breaks for node size", value = "0,0.01,0.05,0.5,5"),
          textAreaInput("palette_group_tree", "Colors for plot", value = "red,blue,orange,green"),
          textInput("taxa_removal_tree", "Remove taxa from plot?",value=NULL),
          selectInput("single_parent_branch_removal","Remove the parent branch if there is only one child branch within the parent branch?",c("False","True")),
          selectInput("single_child_branch_removal","Remove the child branch if there is only one child branch within the parent branch?",c("False","True")),
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
          textInput("test_metadata_order_box", "Type in the order of metadata separated with comma to change those in the figure ",value="default"),
          numericInput("page_box", "Page number", value = 1, min = 1, step = 1),
          textInput("taxa_shown_box", "Select specific taxa", value = ""),
          textAreaInput("palette_group_box", "Colors for plot", value = "red,blue,orange,green"),
          selectInput("x_dir_box", "Direction of X axis labels (1 is horizontal, 2 is vertical)", c(1,2)),
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
          selectInput("cor_method", "Correlation method",c("spearman","pearson","kendall")),
          selectInput("log_norm_cor", "Log normalization?",c("True","False")),
          numericInput("fdr_cutoff_cor", "FDR cutoff \n(Try increasing the cutoff if there is no taxa shown)", value = 0.1, min = 0, step = 0.01),
          numericInput("page_cor", "Page number", value = 1, min = 1, step = 1),
          textInput("taxa_shown_cor", "Select taxa", value = ""),
          textAreaInput("palette_group_cor", "Colors for plot", value = "red,blue,orange,green"),
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
          selectInput("exclude_ASV_p1", "Should ASV be excluded from the analysis (this changes the P-value distribution and FDR)?", c("True","False")),
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
          selectInput("exclude_ASV_p2", "Should ASV be excluded from the analysis (this changes the P-value distribution and FDR)?", c("True","False")),
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
        ),
        mainPanel(
          tabsetPanel(
            tabPanel("Data1", tableOutput("head_p1")),
            tabPanel("Data2", tableOutput("head_p2")),
            tabPanel("Correlation plot", uiOutput("plotPvals.ui"))
          )
        )
      )
    )
  )
)

server <- function(input, output, session) {
  #taxonomic table input
  observeEvent(input$data_type, {

    if(input$data_type == "16S"){
      shinyjs::show(id = "16S")
      shinyjs::hide(id = "wgs")
      shinyjs::hide(id = "one_level")
    }else if(input$data_type == "Metagenomics"){
      shinyjs::show(id = "wgs")
      shinyjs::hide(id = "16S")
      shinyjs::hide(id = "one_level")
    }else{
      shinyjs::show(id = "one_level")
      shinyjs::hide(id = "16S")
      shinyjs::hide(id = "wgs")
    }
  })

  data_raw <- eventReactive(input$button_raw,{
    if(input$data_type=="16S"){
      if(input$sep_16s=="tab"){
        sep_char_16s="\t"
      }else{
        sep_char_16s=input$sep_16s
      }
      format_asv(taxa_file =input$file_16s$datapath,onefile=as.logical(input$onefile),biom=as.logical(input$biom),ASV=as.logical(input$ASV),sep=sep_char_16s,reads_cutoff=as.numeric(input$n_reads_16s),normalization=as.logical(input$norm_16s),rarefy=as.logical(input$rarefy_16s),rarefy_num=as.numeric(input$rarefy_reads_16s))
    }else if (input$data_type=="Metagenomics"){
      if(input$sep_wgs=="tab"){
        sep_char_wgs="\t"
      }else{
        sep_char_wgs=input$sep_wgs
      }
      format_wgs(taxa_file =input$file_wgs$datapath,sep=sep_char_wgs,method=input$method_wgs,reads_cutoff=as.numeric(input$n_reads_wgs),normalization=as.logical(input$norm_wgs),rarefy=as.logical(input$rarefy_wgs),rarefy_num=as.numeric(input$rarefy_reads_wgs))
    }else if (input$data_type=="One level"){
      if(input$sep_path=="tab"){
        sep_char_path="\t"
      }else{
        sep_char_path=input$sep_path
      }
      format_pathway(taxa_file =input$file_path$datapath,sep=sep_char_path,reads_cutoff=as.numeric(input$n_reads_path),normalization=as.logical(input$norm_path),rarefy=as.logical(input$rarefy_path),rarefy_num=as.numeric(input$rarefy_reads_path))
    }
   })

  output$head_raw <- renderTable({
    if(input$data_type=="16S"){
      if(ncol(data_raw())<as.numeric(input$n_raw_16s_col)){
        n_col=ncol(data_raw())
      }else{
        n_col=as.numeric(input$n_raw_16s_col)
      }
      head(data_raw()[,1:n_col], input$n_raw_16s)
    }else if (input$data_type=="Metagenomics"){
      if(ncol(data_raw())<as.numeric(input$n_raw_wgs_col)){
        n_col=ncol(data_raw())
      }else{
        n_col=as.numeric(input$n_raw_wgs_col)
      }
      head(data_raw()[,1:n_col], input$n_raw_wgs)
      }else if (input$data_type=="One level"){
        if(ncol(data_raw())<as.numeric(input$n_raw_path_col)){
          n_col=ncol(data_raw())
        }else{
          n_col=as.numeric(input$n_raw_path_col)
        }
        head(data_raw()[,1:n_col], input$n_raw_path)
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
    if(input$data_type=="16S" | input$data_type=="Metagenomics"){
      "FALSE"
    }else{
      "TRUE"
    }
  })

  # output data dimensions
  output$data_dim <- renderText({
    paste("Number of rows :",nrow(data_raw()),"\nNumber of columns :",ncol(data_raw()))
  })



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
                             choices = c("",colnames(data_meta())),
                             selected = c(""))})
  output$head_stratify_metadata_filter <- renderText({
    if (input$stratify_by_metadata_filter==""){
      ""
    }else{
      head(unique(na.omit(data_meta()[,input$stratify_by_metadata_filter])),n=15)
    }
  })

  observe({updateSelectInput(session, "stratify_by_metadata_mds",
                             choices = c("",colnames(data_meta())),
                             selected = c(""))})
  output$head_stratify_mds <- renderText({
    if (input$stratify_by_metadata_mds==""){
      ""
    }else{
      head(unique(na.omit(data_meta()[,input$stratify_by_metadata_mds])),n=15)
    }
  })


  stratify_by_value_mds_outVar = eventReactive(input$stratify_by_metadata_mds,{
    unique(na.omit(data_meta()[,input$stratify_by_metadata_mds]))
  })

  observe({
    updateSelectInput(session, "stratify_by_value_mds",
                     choices = stratify_by_value_mds_outVar()
    )})

  stratify_by_value_filter_outVar = eventReactive(input$stratify_by_metadata_filter,{
    unique(na.omit(data_meta()[,input$stratify_by_metadata_filter]))
  })

  observe({
    updateSelectInput(session, "stratify_by_value_filter",
                      choices = stratify_by_value_filter_outVar()
    )})



  observe({updateSelectInput(session, "stratify_by_metadata_alpha",
                             choices = c("",colnames(data_meta())),
                             selected = c(""))})
  output$head_stratify_alpha <- renderText({
    if (input$stratify_by_metadata_alpha==""){
      ""
    }else{
      head(unique(na.omit(data_meta()[,input$stratify_by_metadata_alpha])),n=15)
    }
  })


  stratify_by_value_alpha_outVar = eventReactive(input$stratify_by_metadata_alpha,{
    unique(na.omit(data_meta()[,input$stratify_by_metadata_alpha]))
  })

  observe({
    updateSelectInput(session, "stratify_by_value_alpha",
                      choices = stratify_by_value_alpha_outVar()
    )})

  observe({updateSelectInput(session, "stratify_by_metadata_bar",
                             choices = c("",colnames(data_meta())),
                             selected = c(""))})
  output$head_stratify_bar <- renderText({
    if (input$stratify_by_metadata_bar==""){
      ""
    }else{
      head(unique(na.omit(data_meta()[,input$stratify_by_metadata_bar])),n=15)
    }
  })

  stratify_by_value_bar_outVar = eventReactive(input$stratify_by_metadata_bar,{
    unique(na.omit(data_meta()[,input$stratify_by_metadata_bar]))
  })

  observe({
    updateSelectInput(session, "stratify_by_value_bar",
                      choices = stratify_by_value_bar_outVar()
    )})

  stratify_by_value_p1_outVar = eventReactive(input$stratify_by_metadata_p1,{
    unique(na.omit(data_meta()[,input$stratify_by_metadata_p1]))
  })

  observe({
    updateSelectInput(session, "stratify_by_value_p1",
                      choices = stratify_by_value_p1_outVar()
    )})

  stratify_by_value_p2_outVar = eventReactive(input$stratify_by_metadata_p2,{
    unique(na.omit(data_meta()[,input$stratify_by_metadata_p2]))
  })

  observe({
    updateSelectInput(session, "stratify_by_value_p2",
                      choices = stratify_by_value_p2_outVar()
    )})


  #show metadata for PvP
  observe({updateSelectInput(session, "test_metadata_p1",
                             choices = colnames(data_meta()),
                             selected = c(""))})
  output$head_test_p1 <- renderText({
    head(unique(na.omit(data_meta()[,input$test_metadata_p1])),n=15)
  })

  observe({updateSelectInput(session, "stratify_by_metadata_p1",
                             choices = c("",colnames(data_meta())),
                             selected = c(""))})
  output$head_stratify_p1 <- renderText({
    if (input$stratify_by_metadata_p1==""){
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
                             choices = c("",colnames(data_meta())),
                             selected = c(""))})
  output$head_stratify_p2 <- renderText({
    if (input$stratify_by_metadata_p2==""){
      ""
    }else{
      head(unique(na.omit(data_meta()[,input$stratify_by_metadata_p2])),n=15)
    }
  })

  #MDS plot

  plotMDS <- eventReactive(input$button_mds,{
    dataf1=table_subset(taxa_table = data_raw(),metadata=data_meta(),stratify_by_metadata=input$stratify_by_metadata_mds,stratify_by_value=input$stratify_by_value_mds,one_level=as.logical(one_level_all()))
    mds_plot(taxa_table =dataf1,metadata=data_meta(),test_metadata=input$test_metadata_mds,taxa_level=input$taxa_level_mds,method_mds=input$method_mds,one_level=as.logical(one_level_all()),log_norm=as.logical(input$log_normalization_mds),palette_group=strsplit(input$palette_group_mds, ",\\s*")[[1]],distance_type=input$distance_type)
  })

  plotMDS1 <- function(){
    dataf1=table_subset(taxa_table = data_raw(),metadata=data_meta(),stratify_by_metadata=input$stratify_by_metadata_mds,stratify_by_value=input$stratify_by_value_mds,one_level=as.logical(one_level_all()))
    mds_plot(taxa_table =dataf1,metadata=data_meta(),test_metadata=input$test_metadata_mds,taxa_level=input$taxa_level_mds,method_mds=input$method_mds,one_level=as.logical(one_level_all()),log_norm=as.logical(input$log_normalization_mds),palette_group=strsplit(input$palette_group_mds, ",\\s*")[[1]],distance_type=input$distance_type)
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

  plotAlpha <- eventReactive(input$button_alpha,{
    dataf2=table_subset(taxa_table = data_raw(),metadata=data_meta(),stratify_by_metadata=input$stratify_by_metadata_alpha,stratify_by_value=input$stratify_by_value_alpha,one_level=as.logical(one_level_all()))
    alpha_plot(taxa_table =dataf2,metadata=data_meta(),test_metadata=input$test_metadata_alpha,one_level=as.logical(one_level_all()),test_metadata_order=strsplit(input$test_metadata_order_alpha, ",\\s*")[[1]],method=input$method_alpha,xlab_direction=as.integer(input$x_dir_alpha),palette_group=strsplit(input$palette_group_alpha, ",\\s*")[[1]])
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
    alpha_plot(taxa_table =dataf2,metadata=data_meta(),test_metadata=input$test_metadata_alpha,one_level=as.logical(one_level_all()),test_metadata_order=strsplit(input$test_metadata_order_alpha, ",\\s*")[[1]],method=input$method_alpha,xlab_direction=as.integer(input$x_dir_alpha),palette_group=strsplit(input$palette_group_alpha, ",\\s*")[[1]])
  }

  output$plotAlphaDownload <- downloadHandler(
    filename ="alpha_diversity.pdf",
    content = function(file) {
      pdf(file, height = 30*alpha_h(),width=15*alpha_w())
      plotAlpha1()
      dev.off()
    },contentType = "image/pdf")


  #Taxa barplot

  plotBar <- eventReactive(input$button_bar,{
    dataf2=table_subset(taxa_table = data_raw(),metadata=data_meta(),stratify_by_metadata=input$stratify_by_metadata_bar,stratify_by_value=input$stratify_by_value_bar,one_level=as.logical(one_level_all()))
    #boxplot(dataf2[1,]~data_meta()[,19],main=input$test_metadata_bar)
    taxa_barplot(taxa_table =dataf2,metadata=data_meta(),test_metadata=input$test_metadata_bar,one_level=as.logical(one_level_all()),num_taxa=as.integer(input$num_taxa_bar),test_metadata_order=strsplit(input$test_metadata_order_bar, ",\\s*")[[1]],taxa_level=input$taxa_level_bar,xlab_direction=as.integer(input$x_dir_bar),legend_size=as.numeric(input$legend_size_bar),palette_group=strsplit(input$palette_group_bar, ",\\s*")[[1]])
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
    taxa_barplot(taxa_table =dataf2,metadata=data_meta(),test_metadata=input$test_metadata_bar,one_level=as.logical(one_level_all()),num_taxa=as.integer(input$num_taxa_bar),test_metadata_order=strsplit(input$test_metadata_order_bar, ",\\s*")[[1]],taxa_level=input$taxa_level_bar,xlab_direction=as.integer(input$x_dir_bar),legend_size=as.numeric(input$legend_size_bar),palette_group=strsplit(input$palette_group_bar, ",\\s*")[[1]])
   }

  output$plotBarDownload <- downloadHandler(
    filename = "barplot.pdf",
    content = function(file) {
      pdf(file, height = 12*bar_h(),width=12*bar_w())
      plotBar1()
      dev.off()
    },contentType = "image/pdf")


  #data filter/subset
  data_filtered <- eventReactive(input$button_filter,{
    table_subset(taxa_table = data_raw(),metadata=data_meta(),stratify_by_metadata=input$stratify_by_metadata_filter,stratify_by_value=input$stratify_by_value_filter,prevalence_cutoff=as.numeric(input$prevalence_cutoff), abundance_cutoff=as.numeric(input$abundance_cutoff),one_level=as.logical(one_level_all()),exclude_ASV=as.logical(input$exclude_ASV_filter))
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
    tree_view(taxa_table =data_filtered(),metadata=data_meta(),fdrs=data_fdrs(),test_metadata=input$test_metadata_stat,fdr_cutoff=input$fdr_cutoff_tree,test_metadata_continuous=as.logical(input$test_metadata_continuous_tree),
              node_size_breaks=as.numeric(strsplit(input$node_size_breaks, ",\\s*")[[1]]),palette_highlight=strsplit(input$palette_group_tree, ",\\s*")[[1]],single_parent_branch_removal=as.logical(input$single_parent_branch_removal),single_child_branch_removal=as.logical(input$single_child_branch_removal),
              prevalence_cutoff=as.numeric(input$prevalence_cutoff_tree), abundance_cutoff=as.numeric(input$abundance_cutoff_tree),taxa_removal=input$taxa_removal_tree)
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
    tree_view(taxa_table =data_filtered(),metadata=data_meta(),fdrs=data_fdrs(),test_metadata=input$test_metadata_stat,fdr_cutoff=input$fdr_cutoff_tree,test_metadata_continuous=as.logical(input$test_metadata_continuous_tree),
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


  #Box plot

  plotBox <- eventReactive(input$button_box,{
    taxa_boxplot(taxa_table =data_filtered(),metadata=data_meta(),fdrs=data_fdrs(),test_metadata=input$test_metadata_stat,test_metadata_order=strsplit(input$test_metadata_order_box, ",\\s*")[[1]],one_level=as.logical(one_level_all()),
                 log_norm=input$log_norm_box,taxa_shown=input$taxa_shown_box,cutoff=input$fdr_cutoff_box,page=input$page_box,
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
    plotOutput("plotBox", height = 800*box_h(),width=1000*box_w())
  })


  plotBox1 <- function(){
    taxa_boxplot_download(taxa_table =data_filtered(),metadata=data_meta(),fdrs=data_fdrs(),test_metadata=input$test_metadata_stat,test_metadata_order=strsplit(input$test_metadata_order_box, ",\\s*")[[1]],one_level=as.logical(one_level_all()),
                 log_norm=input$log_norm_box,taxa_shown=input$taxa_shown_box,cutoff=input$fdr_cutoff_box,
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
                 col_metadata=input$col_metadata_cor,log_norm=input$log_norm_cor,taxa_shown=input$taxa_shown_cor,page=input$page_cor,
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
                 col_metadata=input$col_metadata_cor,log_norm=input$log_norm_cor,taxa_shown=input$taxa_shown_cor,
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
      dataf1_p1=table_subset(taxa_table = data_raw(),metadata=data_meta(),stratify_by_metadata=input$stratify_by_metadata_p1,stratify_by_value=input$stratify_by_value_p1,prevalence_cutoff=as.numeric(input$prevalence_cutoff_p1), abundance_cutoff=as.numeric(input$abundance_cutoff_p1),one_level=as.logical(one_level_all()),exclude_ASV=as.logical(input$exclude_ASV_p1))
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
      dataf1_p2=table_subset(taxa_table = data_raw(),metadata=data_meta(),stratify_by_metadata=input$stratify_by_metadata_p2,stratify_by_value=input$stratify_by_value_p2,prevalence_cutoff=as.numeric(input$prevalence_cutoff_p2), abundance_cutoff=as.numeric(input$abundance_cutoff_p2),one_level=as.logical(one_level_all()),exclude_ASV=as.logical(input$exclude_ASV_p2))
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
}

# Run the application
shinyApp(ui = ui, server = server)


