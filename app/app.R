library(shiny)
library(dplyr)
library(stringr)
library(XML)
library(rvest)
library(httr)

# Shiny UI -------
ui <- fluidPage(

  navbarPage(
    "plotmicrobiome",
    id = "main_navbar",

    tabPanel(
      "Data input",
      sidebarLayout(
        sidebarPanel(
          width = 3,
          br(),
          selectInput("data_type", "Is this 16S data or metagenomics data?", c("16S","Metagenomics")),
          br(),
          textInput("dir_16s", "16S data directory",value="/Users/shansun/git/plotmicrobiome/data-raw/biom_taxonomy.biom"),
          selectInput("onefile", "Is the data in one file?", c("True","False")),
          selectInput("biom", "Is the data in biom format?", c("True","False")),
          selectInput("ASV", "Is the data a DADA2 ASV table?", c("True","False")),
          textInput("sep_16s", "What is the delimiter?",value=","),
          numericInput("n_raw", "Rows", value = 5, min = 1, step = 1),
          br(),
          br(),
          textInput("dir_wgs", "Metagenomics data directory",value=NULL),
          selectInput("method_wgs", "Which tool was used for taxonomic classification?", c("kraken2","metaphlan2")),
          textInput("wgs_sep", "What is the delimiter of the table?", value=","),
          br(),
          actionButton("button_raw", "Run"),
          br()
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
          br(),
          textInput("dir_meta", "Metadata directory",value="/Users/shansun/git/plotmicrobiome/data-raw/metadata_cafe.csv"),
          textInput("metadata_sep", "What is the delimiter?",value=","),
          numericInput("meta_sample_name_col", "Which column are the sample names in?", value = 0),
          numericInput("n_meta", "Rows", value = 5, min = 1, step = 1),
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
          br(),
          textInput("stratify_by_metadata_mds", "Select the variable for stratification",value=NULL),
          textInput("stratify_by_value_mds", "Select the value for stratification",value=NULL),
          textInput("test_metadata_mds", "Select the metadata for testing",value=NULL),
          selectInput("method_mds", "Which method should be used for ordination?", c("pcoa","nmds")),
          textInput("distance_type", "Distance type",value="bray"),
          textInput("xaxis_mds", "Which MDS for x axis",value=1),
          textInput("yaxis_mds", "Which MDS for y axis",value=2),
          textAreaInput("palette_group_mds", "Colors for plot", value = "red,blue,orange,green"),
          actionButton("button_mds", "Run"),
          br(),
          br()
        ),
        mainPanel(plotOutput(outputId = "plotMDS",height=700,width=800))
      )
    ),

    tabPanel(
      "Alpha diversity plot",
      sidebarLayout(
        sidebarPanel(
          width = 3,
          br(),
          br(),
          textInput("stratify_by_metadata_alpha", "Select the variable for stratification",value=NULL),
          textInput("stratify_by_value_alpha", "Select the value for stratification",value=NULL),
          textInput("test_metadata_alpha", "Select the metadata for testing",value=NULL),
          textInput("taxa_level_alpha", "Taxonomic level",value="species"),
          selectInput("method_alpha", "Select method",c("wilcoxon","t.test","kruskal-wallis","anova")),
          textAreaInput("palette_group_alpha", "Colors for plot", value = "red,blue,orange,green"),
          textInput("x_dir_alpha", "Direction of X axis labels", value = 1),
          actionButton("button_alpha", "Run"),
          br(),
          br()
        ),
        mainPanel(plotOutput(outputId = "plotAlpha",height=700,width=800))
      )
    ),

    tabPanel(
      "Data filter",
      sidebarLayout(
        sidebarPanel(
          width = 3,
          br(),
          br(),
          textInput("stratify_by_metadata", "Select the variable for stratification",value="Study"),
          textInput("stratify_by_value", "Select the value for stratification",value="Sugar"),
          numericInput("prevalence_cutoff", "Prevalence cutoff", value = 0.25),
          numericInput("abundance_cutoff", "Abundance_cutoff", value = 0),
          numericInput("n_filtered", "Rows", value = 5, min = 1, step = 1),
          actionButton("button_filter", "Run"),
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
          width = 3,
          br(),
          br(),
          selectInput("method_stat", "Select method",c("wilcoxon","t.test","kruskal-wallis","anova")),
          textInput("test_metadata_stat", "Select the metadata for testing",value="Timepoint"),
          selectInput("sort_fdr", "Sort by FDR",c("True","False")),
          numericInput("n_fdrs", "Rows", value = 5, min = 1, step = 1),
          actionButton("button_fdrs", "Run"),
          br(),
          br()
        ),
        mainPanel(tableOutput("head_fdrs"))
      )
    ),

    tabPanel(
      "Boxplot plot",
      sidebarLayout(
        sidebarPanel(
          width = 3,
          br(),
          br(),
          selectInput("log_norm_box", "Log normalization?",c("True","False")),
          numericInput("fdr_cutoff_box", "FDR cutoff", value = 0.1, min = 0, step = 0.01),
          numericInput("page_box", "Page number", value = 1, min = 1, step = 1),
          textInput("taxa_shown_box", "Select taxa", value = ""),
          textAreaInput("palette_group_box", "Colors for plot", value = "red,blue,orange,green"),
          textInput("x_dir_box", "Direction of X axis labels", value = 1),
          actionButton("button_box", "Run"),
          br(),
          br(),
          br(),
          br(),
          br(),
          br(),
          br(),
          br(),
          br(),
          br(),
          br(),
          br()
        ),
        mainPanel(
          plotOutput(outputId = "plotBox",height=700,width=800))
      )
    ),

    tabPanel(
      "Tree plot",
      sidebarLayout(
        sidebarPanel(
          width = 3,
          br(),
          br(),
          numericInput("fdr_cutoff_tree", "FDR cutoff", value = 0.1, min = 0, step = 0.01),
          textAreaInput("node_size_breaks", "Breaks for node size", value = "0,0.01,0.05,0.5,5"),
          textAreaInput("palette_group_tree", "Colors for plot", value = "red,blue,orange,green"),
          textInput("taxa_removal_tree", "Remove taxa from plot?",value=NULL),
          actionButton("button_tree", "Run"),
          br(),
          br()
        ),
        #mainPanel(tableOutput("head_tree"))
        mainPanel(plotOutput(outputId = "plotTree"))
      )
    ),

    tabPanel(
      "Correlation plot",
      sidebarLayout(
        sidebarPanel(
          width = 3,
          br(),
          br(),
          textInput("test_metadata_cor", "Select the metadata for testing",value="test_score"),
          textInput("col_metadata_cor", "Select the metadata for color",value="Timepoint"),
          selectInput("cor_method", "Correlation method",c("spearman","pearson","kendall")),
          selectInput("log_norm_cor", "Log normalization?",c("True","False")),
          numericInput("fdr_cutoff_cor", "FDR cutoff", value = 0.1, min = 0, step = 0.01),
          numericInput("page_cor", "Page number", value = 1, min = 1, step = 1),
          textInput("taxa_shown_cor", "Select taxa", value = ""),
          textAreaInput("palette_group_cor", "Colors for plot", value = "red,blue,orange,green"),
          actionButton("button_cor", "Run"),
          br(),
          br(),
          br(),
          br(),
          br(),
          br(),
          br(),
          br(),
          br(),
          br(),
          br(),
          br()
        ),
        mainPanel(
          plotOutput(outputId = "plotCor",height=700,width=800))
      )
    )
  )
)

server <- function(input, output, session) {
  #taxonomic table input
  data_raw <- eventReactive(input$button_raw,{
    if(input$data_type=="16S"){
      format_asv(taxa_file =input$dir_16s,onefile=as.logical(input$onefile),biom=as.logical(input$biom),ASV=as.logical(input$ASV),sep=input$sep_16s)
    }else if (input$data_type=="Metagenomics"){
      format_wgs(taxa_file =input$dir_wgs,sep=input$wgs_sep,method=input$method_wgs)
    }
   })

  output$head_raw <- renderTable({
    head(data_raw(), input$n_raw)
  },rownames = TRUE)

  #metadata input
  data_meta <- eventReactive(input$button_meta,{
    meta_format(metadata=input$dir_meta,metadata_sep=input$metadata_sep,meta_sample_name_col=input$meta_sample_name_col)
  })

  output$head_meta <- renderTable({
    head(data_meta(), input$n_meta)
  },rownames = TRUE)

  #MDS plot

  plotMDS <- eventReactive(input$button_mds,{
    dataf1=table_subset(taxa_table = data_raw(),metadata=data_meta(),stratify_by_metadata=input$stratify_by_metadata_mds,stratify_by_value=input$stratify_by_value_mds)
    mds_plot(taxa_table =dataf1,metadata=data_meta(),test_metadata=input$test_metadata_mds,method_mds=input$method_mds,xaxis=as.numeric(input$xaxis_mds),yaxis=as.numeric(input$yaxis_mds),palette_group=strsplit(input$palette_group_mds, ",\\s*")[[1]],distance_type=input$distance_type)
  })

  output$plotMDS <- renderPlot({

    plotMDS()

  })

  #Alpha diversity plot

  plotAlpha <- eventReactive(input$button_alpha,{
    dataf2=table_subset(taxa_table = data_raw(),metadata=data_meta(),stratify_by_metadata=input$stratify_by_metadata_alpha,stratify_by_value=input$stratify_by_value_alpha)
    alpha_plot(taxa_table =dataf2,metadata=data_meta(),test_metadata=input$test_metadata_alpha,taxa_level=input$taxa_level_alpha,method=input$method_alpha,xlab_direction=input$x_dir_alpha,palette_group=strsplit(input$palette_group_alpha, ",\\s*")[[1]])
  })

  output$plotAlpha <- renderPlot({

    plotAlpha()

  })

  #data filter/subset
  data_filtered <- eventReactive(input$button_filter,{
    table_subset(taxa_table = data_raw(),metadata=data_meta(),stratify_by_metadata=input$stratify_by_metadata,stratify_by_value=input$stratify_by_value,prevalence_cutoff=input$prevalence_cutoff, abundance_cutoff=input$abundance_cutoff)
   })
  output$head_filtered <- renderTable({
    head(data_filtered(), input$n_filtered)
  },rownames = TRUE)


  #statistical test
  data_fdrs <- eventReactive(input$button_fdrs,{
    stat_test(taxa_table =data_filtered(),metadata=data_meta(),test_metadata=input$test_metadata_stat,method=input$method_stat)
   })

  data_fdrs1 <- reactive({
    if(input$sort_fdr){
      data_fdrs()[order(data_fdrs()[,2]),]
    }else{
      data_fdrs()
      }
  })

  output$head_fdrs <- renderTable({
    head(data_fdrs1(), input$n_fdrs)
  },rownames = TRUE,digits=-2)

  #Box plot

  plotBox <- eventReactive(input$button_box,{
    taxa_boxplot(taxa_table =data_filtered(),metadata=data_meta(),fdrs=data_fdrs(),test_metadata=input$test_metadata_stat,log_norm=input$log_norm_box,taxa_shown=input$taxa_shown_box,page=input$page_box,cutoff=input$fdr_cutoff_box,xlab_direction=input$x_dir_box,palette_group=strsplit(input$palette_group_box, ",\\s*")[[1]])
  })

  output$plotBox <- renderPlot({

    plotBox()

  })

  #Tree plot

  plotTree <- eventReactive(input$button_tree,{
    tree_view(taxa_table =data_filtered(),metadata=data_meta(),fdrs=data_fdrs(),test_metadata=input$test_metadata_stat,fdr_cutoff=input$fdr_cutoff_tree,node_size_breaks=as.numeric(strsplit(input$node_size_breaks, ",\\s*")[[1]]),palette_highlight=strsplit(input$palette_group_tree, ",\\s*")[[1]],taxa_removal=input$taxa_removal_tree)
  })

  output$plotTree <- renderPlot({

    plotTree()

  })

  #correlation plot

  plotCor <- eventReactive(input$button_cor,{
    meta_corplot(taxa_table =data_filtered(),metadata=data_meta(),test_metadata=input$test_metadata_cor,cor_method=input$cor_method,
                 col_metadata=input$col_metadata_cor,log_norm=input$log_norm_cor,taxa_shown=input$taxa_shown_cor,page=input$page_cor,
                 fdr_cutoff=input$fdr_cutoff_cor,palette_group=strsplit(input$palette_group_cor, ",\\s*")[[1]])
  })

  output$plotCor<- renderPlot({

    plotCor()

  })

}



# Run the application
shinyApp(ui = ui, server = server)
