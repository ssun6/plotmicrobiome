library(shiny)
#library(dplyr)
#library(stringr)
#library(rvest)
#library(httr)
#library(XML)


# Shiny UI -------
#increase upload size
options(shiny.maxRequestSize = 100*1024^2)
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
          fileInput("file_16s", "Choose file (accept .csv, .tsv and .biom files.)"),
          selectInput("onefile", "Is the data in one file?", c("True","False")),
          selectInput("biom", "Is the data in biom format?", c("True","False")),
          selectInput("ASV", "Is the data a DADA2 ASV table?", c("True","False")),
          textInput("sep_16s", "What is the delimiter?",value=","),
          numericInput("n_raw_16s", "Preview rows", value = 5, min = 1, step = 1),
          br(),
          br(),
          fileInput("file_wgs", "Choose file"),
          #textInput("dir_wgs", "Metagenomics data directory",value=NULL),
          selectInput("method_wgs", "Which tool was used for taxonomic classification?", c("kraken2","metaphlan2")),
          textInput("sep_wgs", "What is the delimiter of the table?", value="\t"),
          numericInput("n_raw_wgs", "Preview rows", value = 5, min = 1, step = 1),
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
          fileInput("file_meta", "Choose metadata file"),
          #textInput("dir_meta", "Metadata directory",value="./data-raw/metadata_cafe.csv"),
          textInput("metadata_sep", "What is the delimiter?",value=","),
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
          br(),
          selectInput("stratify_by_metadata_mds", "Select the variable for stratification",c("")),
          h5("Variable preview:"),
          textOutput("head_stratify_mds"),
          br(),
          textInput("stratify_by_value_mds", "Select the value for stratification",value=NULL),
          selectInput("test_metadata_mds", "Select the metadata for testing",c("")),
          h5("Variable preview:"),
          textOutput("head_test_mds"),
          br(),
          selectInput("method_mds", "Which method should be used for ordination?", c("pcoa","nmds")),
          selectInput("distance_type", "Distance type",c("bray","euclidean","manhattan","jaccard")),
          textAreaInput("palette_group_mds", "Colors for plot", value = "red,blue,orange,green"),
          actionButton("button_mds", "Run"),
          br(),
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
          br(),
          selectInput("stratify_by_metadata_alpha", "Select the variable for stratification",c("")),
          h5("Variable preview:"),
          textOutput("head_stratify_alpha"),
          br(),
          textInput("stratify_by_value_alpha", "Select the value for stratification",value=NULL),
          selectInput("test_metadata_alpha", "Select the metadata for testing",c("")),
          h5("Variable preview:"),
          textOutput("head_test_alpha"),
          br(),
          selectInput("method_alpha", "Select method",c("wilcoxon","t.test","kruskal-wallis","anova")),
          textAreaInput("palette_group_alpha", "Colors for plot", value = "red,blue,orange,green"),
          textInput("x_dir_alpha", "Direction of X axis labels", value = 1),
          actionButton("button_alpha", "Run"),
          br(),
          br(),
          h5("Download figure:"),
          downloadButton("plotAlphaDownload", "Download"),
          br(),
          br()
        ),
        mainPanel(plotOutput(outputId = "plotAlpha",height=1500,width=1000))
      )
    ),

    tabPanel(
      "Data filter",
      sidebarLayout(
        sidebarPanel(
          width = 3,
          br(),
          br(),
          selectInput("stratify_by_metadata", "Select the variable for stratification",c("")),
          h5("Variable preview:"),
          textOutput("head_stratify_metadata"),
          br(),
          textInput("stratify_by_value", "Select the value for stratification",value = NULL),
          numericInput("prevalence_cutoff", "Prevalence cutoff", value = 0.25),
          numericInput("abundance_cutoff", "Abundance cutoff", value = 0),
          numericInput("n_filtered", "Preview rows", value = 5, min = 1, step = 1),
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
          selectInput("test_metadata_stat", "Select the metadata for testing",c("")),
          h5("Variable preview:"),
          textOutput("head_test_stat"),
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
          textInput("taxa_shown_box", "Select specific taxa", value = ""),
          textAreaInput("palette_group_box", "Colors for plot", value = "red,blue,orange,green"),
          textInput("x_dir_box", "Direction of X axis labels", value = 1),
          actionButton("button_box", "Run"),
          br(),
          br(),
          h5("Download figure:"),
          downloadButton("plotBoxDownload", "Download"),
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
          h5("Parameters for pruning tree:"),
          numericInput("prevalence_cutoff_tree", "Prevalence cutoff", value = 0.1,min=0),
          numericInput("abundance_cutoff_tree", "Abundance cutoff", value = 0,min=0),
          br(),
          br(),
          numericInput("fdr_cutoff_tree", "FDR cutoff", value = 0.1, min = 0, step = 0.01),
          textAreaInput("node_size_breaks", "Breaks for node size", value = "0,0.01,0.05,0.5,5"),
          textAreaInput("palette_group_tree", "Colors for plot", value = "red,blue,orange,green"),
          textInput("taxa_removal_tree", "Remove taxa from plot?",value=NULL),
          actionButton("button_tree", "Run"),
          br(),
          br(),
          h5("Download figure:"),
          downloadButton("plotTreeDownload", "Download"),
          br(),
          br()
        ),
        #mainPanel(tableOutput("head_tree"))
        mainPanel(plotOutput(outputId = "plotTree",width=800,height=500))
      )
    ),

    tabPanel(
      "Correlation plot",
      sidebarLayout(
        sidebarPanel(
          width = 3,
          br(),
          br(),
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
          h5("Download figure:"),
          downloadButton("plotCorDownload", "Download"),
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
      if(input$sep_16s=="\\t"){
        sep_char_16s="\t"
      }else{
        sep_char_16s=input$sep_16s
      }
      format_asv(taxa_file =input$file_16s$datapath,onefile=as.logical(input$onefile),biom=as.logical(input$biom),ASV=as.logical(input$ASV),sep=sep_char_16s)
    }else if (input$data_type=="Metagenomics"){
      if(input$sep_wgs=="\\t"){
        sep_char_wgs="\t"
      }else{
        sep_char_wgs=input$sep_wgs
      }
      format_wgs(taxa_file =input$file_wgs$datapath,sep=sep_char_wgs,method=input$method_wgs)
    }
   })

  output$head_raw <- renderTable({
    if(input$data_type=="16S"){
      head(data_raw(), input$n_raw_16s)
    }else if (input$data_type=="Metagenomics"){
      head(data_raw(), input$n_raw_wgs)
      }},rownames = TRUE)

  #metadata input
  data_meta <- eventReactive(input$button_meta,{
    meta_format(metadata=input$file_meta$datapath,metadata_sep=input$metadata_sep,meta_sample_name_col=input$meta_sample_name_col)
  })


  output$head_meta <- renderTable({
    head(data_meta(), input$n_meta)
  },rownames = TRUE)

  # output metadata variables for test
  observe({updateSelectInput(session, "test_metadata_mds",
                    choices = colnames(data_meta()),
                    selected = c("Case_Ctrl"))})
  output$head_test_mds <- renderText({
    head(unique(na.omit(data_meta()[,input$test_metadata_mds])),n=15)
  })

  observe({updateSelectInput(session, "test_metadata_alpha",
                             choices = colnames(data_meta()),
                             selected = c("Case_Ctrl"))})
  output$head_test_alpha <- renderText({
    head(unique(na.omit(data_meta()[,input$test_metadata_alpha])),n=15)
  })

  observe({updateSelectInput(session, "test_metadata_stat",
                             choices = colnames(data_meta()),
                             selected = c("Case_Ctrl"))})
  output$head_test_stat <- renderText({
    head(unique(na.omit(data_meta()[,input$test_metadata_stat])),n=15)
  })

  observe({updateSelectInput(session, "col_metadata_cor",
                             choices = colnames(data_meta()),
                             selected = c("Case_Ctrl"))})
  output$head_metadata_cor_col <- renderText({
    head(unique(na.omit(data_meta()[,input$col_metadata_cor])),n=15)
  })

  # output metadata variables for correlation
  observe({updateSelectInput(session, "test_metadata_cor",
                             choices = colnames(data_meta()),
                             selected = c("age"))})
  output$head_metadata_cor_test <- renderText({
    head(unique(na.omit(data_meta()[,input$test_metadata_cor])),n=15)
  })

  # output metadata variables for stratify
  observe({updateSelectInput(session, "stratify_by_metadata",
                             choices = c("",colnames(data_meta())),
                             selected = c(""))})
  output$head_stratify_metadata <- renderText({
    if (input$stratify_by_metadata_mds==""){
      ""
    }else{
      head(unique(na.omit(data_meta()[,input$stratify_by_metadata])),n=15)
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

  #MDS plot

  plotMDS <- eventReactive(input$button_mds,{
    dataf1=table_subset(taxa_table = data_raw(),metadata=data_meta(),stratify_by_metadata=input$stratify_by_metadata_mds,stratify_by_value=input$stratify_by_value_mds)
    mds_plot(taxa_table =dataf1,metadata=data_meta(),test_metadata=input$test_metadata_mds,method_mds=input$method_mds,palette_group=strsplit(input$palette_group_mds, ",\\s*")[[1]],distance_type=input$distance_type)
  })

  plotMDS1 <- function(){
    dataf1=table_subset(taxa_table = data_raw(),metadata=data_meta(),stratify_by_metadata=input$stratify_by_metadata_mds,stratify_by_value=input$stratify_by_value_mds)
    mds_plot(taxa_table =dataf1,metadata=data_meta(),test_metadata=input$test_metadata_mds,method_mds=input$method_mds,palette_group=strsplit(input$palette_group_mds, ",\\s*")[[1]],distance_type=input$distance_type)
  }

  output$plotMDS <- renderPlot({

    plotMDS()

  })

  output$plotMDS.ui <- renderUI({
    if(input$method_mds=="pcoa"){
      plotOutput("plotMDS", height = 1500,width=600)
    }else{
      plotOutput("plotMDS", height = 500,width=600)
    }
  })

  output$plotMDSDownload <- downloadHandler(
    filename = paste0(input$test_metadata_mds,"_MDS.pdf"),
    content = function(file) {
      if(input$method_mds=="pcoa"){
        pdf(file, height = 23,width=8,onefile=T)
      }else{
        pdf(file, height = 8,width=8)
      }
      plotMDS1()
      dev.off()
    },contentType = "image/pdf")

  #Alpha diversity plot

  plotAlpha <- eventReactive(input$button_alpha,{
    dataf2=table_subset(taxa_table = data_raw(),metadata=data_meta(),stratify_by_metadata=input$stratify_by_metadata_alpha,stratify_by_value=input$stratify_by_value_alpha)
    alpha_plot(taxa_table =dataf2,metadata=data_meta(),test_metadata=input$test_metadata_alpha,method=input$method_alpha,xlab_direction=input$x_dir_alpha,palette_group=strsplit(input$palette_group_alpha, ",\\s*")[[1]])
  })

  output$plotAlpha <- renderPlot({

    plotAlpha()

  })

  plotAlpha1 <- function(){
    dataf2=table_subset(taxa_table = data_raw(),metadata=data_meta(),stratify_by_metadata=input$stratify_by_metadata_alpha,stratify_by_value=input$stratify_by_value_alpha)
    alpha_plot(taxa_table =dataf2,metadata=data_meta(),test_metadata=input$test_metadata_alpha,method=input$method_alpha,xlab_direction=input$x_dir_alpha,palette_group=strsplit(input$palette_group_alpha, ",\\s*")[[1]])
  }

  output$plotAlphaDownload <- downloadHandler(
    filename = paste0(input$test_metadata_alpha,"_alpha_diversity.pdf"),
    content = function(file) {
      pdf(file, height = 18,width=12)
      plotAlpha1()
      dev.off()
    },contentType = "image/pdf")

  #data filter/subset
  data_filtered <- eventReactive(input$button_filter,{
    table_subset(taxa_table = data_raw(),metadata=data_meta(),stratify_by_metadata=input$stratify_by_metadata,stratify_by_value=input$stratify_by_value,prevalence_cutoff=as.numeric(input$prevalence_cutoff), abundance_cutoff=as.numeric(input$abundance_cutoff))
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

  output$downloadStat <- downloadHandler(
    filename = function() {
      paste(input$test_metadata_stat,"_",input$method_stat, ".csv", sep = "")
    },
    content = function(file) {
      write.csv(data_fdrs1(), file, row.names = FALSE)
    }
  )

  #Box plot

  plotBox <- eventReactive(input$button_box,{
    taxa_boxplot(taxa_table =data_filtered(),metadata=data_meta(),fdrs=data_fdrs(),test_metadata=input$test_metadata_stat,
                 log_norm=input$log_norm_box,taxa_shown=input$taxa_shown_box,cutoff=input$fdr_cutoff_box,page=input$page_box,
                 xlab_direction=input$x_dir_box,palette_group=strsplit(input$palette_group_box, ",\\s*")[[1]])
  })

  output$plotBox <- renderPlot({

    plotBox()

  })

  plotBox1 <- function(){
    taxa_boxplot_download(taxa_table =data_filtered(),metadata=data_meta(),fdrs=data_fdrs(),test_metadata=input$test_metadata_stat,
                 log_norm=input$log_norm_box,taxa_shown=input$taxa_shown_box,cutoff=input$fdr_cutoff_box,
                 xlab_direction=input$x_dir_box,palette_group=strsplit(input$palette_group_box, ",\\s*")[[1]])
  }

  output$plotBoxDownload <- downloadHandler(
    filename = paste0(input$test_metadata_stat,"_boxplots.pdf"),
    content = function(file) {
      pdf(file,height=10,width=10,onefile = T)
      plotBox1()
      dev.off()
    })

  #Tree plot

  plotTree <- eventReactive(input$button_tree,{
    tree_view(taxa_table =data_raw(),metadata=data_meta(),fdrs=data_fdrs(),test_metadata=input$test_metadata_stat,fdr_cutoff=input$fdr_cutoff_tree,
              node_size_breaks=as.numeric(strsplit(input$node_size_breaks, ",\\s*")[[1]]),palette_highlight=strsplit(input$palette_group_tree, ",\\s*")[[1]],
              prevalence_cutoff=as.numeric(input$prevalence_cutoff_tree), abundance_cutoff=as.numeric(input$abundance_cutoff_tree),taxa_removal=input$taxa_removal_tree)
  })


  output$plotTree <- renderPlot({

    plotTree()

  })

  plotTree1 <- function(){
    tree_view(taxa_table =data_raw(),metadata=data_meta(),fdrs=data_fdrs(),test_metadata=input$test_metadata_stat,fdr_cutoff=input$fdr_cutoff_tree,
              node_size_breaks=as.numeric(strsplit(input$node_size_breaks, ",\\s*")[[1]]),palette_highlight=strsplit(input$palette_group_tree, ",\\s*")[[1]],
              prevalence_cutoff=as.numeric(input$prevalence_cutoff_tree), abundance_cutoff=as.numeric(input$abundance_cutoff_tree),taxa_removal=input$taxa_removal_tree)
  }

  output$plotTreeDownload <- downloadHandler(
    filename = paste0(input$test_metadata_stat,"_tree.pdf"),
    content = function(file) {
      pdf(file,height=12,width=9)
      print(plotTree1 ())
      dev.off()
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

  plotCor1 <- function(){
    meta_corplot_download(taxa_table =data_filtered(),metadata=data_meta(),test_metadata=input$test_metadata_cor,cor_method=input$cor_method,
                 col_metadata=input$col_metadata_cor,log_norm=input$log_norm_cor,taxa_shown=input$taxa_shown_cor,
                 fdr_cutoff=input$fdr_cutoff_cor,palette_group=strsplit(input$palette_group_cor, ",\\s*")[[1]])
  }

  output$plotCorDownload <- downloadHandler(
    filename = paste0("correlation.pdf"),
    content = function(file) {
      pdf(file,height=15,width=15,onefile = T)
      print(plotCor1())
      dev.off()
    })
}

# Run the application
shinyApp(ui = ui, server = server)


