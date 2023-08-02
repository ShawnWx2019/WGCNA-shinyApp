###############################
#	prj: shiny app
#	Assignment: WGCNA by click shiny app
#	Author: Shawn Wang
#	Date: Jan 21, 2023
# Version: V0.1.0
###############################
options("repos" = c(CRAN="https://mirrors.tuna.tsinghua.edu.cn/CRAN/"))
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
if (!require('devtools')) install.packages('devtools');
if (!require('DESeq2')) BiocManager::install('DESeq2',update = FALSE);
if (!require('shinyjs')) install.packages('shinyjs');
if (!require('dashboardthemes')) install.packages('dashboardthemes');
if (!require('shinydashboard')) install.packages('shinydashboard');
if (!require("DT")) install.packages('DT');
if (!require('shiny')) install.packages('shiny');
if (!require('ggpmisc')) install.packages('ggpmisc');
if (!require('dplyr')) install.packages('dplyr');
if (!require('GO.db')) BiocManager::install('GO.db',update = FALSE);
if (!require('WGCNA')) BiocManager::install('WGCNA',update = FALSE);
if (!require('ComplexHeatmap')) BiocManager::install('ComplexHeatmap',update = FALSE);
if (!require('circlize')) BiocManager::install('circlize',update = FALSE);
if (!require('stringr')) install.packages('stringr');
if (!require('ape')) install.packages('ape');
if (!require('reshape2')) install.packages('reshape2');
if (!require('edgeR')) BiocManager::install('edgeR',update = FALSE);
if (!require('shinythemes')) install.packages('shinythemes');
if (!require('ggplotify')) install.packages('ggplotify');
if (!require('ggprism')) install.packages('ggprism');
if (!require('ggpubr')) install.packages('ggpubr');
if (!require('patchwork')) install.packages('patchwork');
if (!require('tidyverse')) install.packages('tidyverse');
if (!require('shinyjqui')) install.packages('shinyjqui');
if (!require('colourpicker')) install.packages('colourpicker');
if (!require('conflicted')) BiocManager::install('conflicted',update = FALSE);
suppressMessages(library(devtools))
if (!require('ShinyWGCNA')) devtools::install_github("ShawnWx2019/WGCNAShinyFun",ref = "master");
suppressMessages(library(ShinyWGCNA))
suppressMessages(library(shinyjs))
suppressMessages(library(dashboardthemes))
suppressMessages(library(shinydashboard))
suppressMessages(library(DT))
suppressMessages(library(shiny))
suppressMessages(library(DESeq2))
suppressMessages(library(ggplot2))
suppressMessages(library(dplyr))
suppressMessages(library(WGCNA))
suppressMessages(library(stringr))
suppressMessages(library(ape))
suppressMessages(library(ComplexHeatmap))
suppressMessages(library(circlize))
suppressMessages(library(reshape2))
suppressMessages(library(edgeR))
suppressMessages(library(shinythemes))
suppressMessages(library(ggplotify))
suppressMessages(library(ggprism))
suppressMessages(library(patchwork))
suppressMessages(library(tidyverse))
suppressMessages(library(shinyjqui))
suppressMessages(library(ggpubr))
suppressMessages(library(conflicted))
options(shiny.maxRequestSize = 300*1024^2)
options(scipen = 6)
conflict_prefer("select","dplyr")
conflict_prefer("filter","dplyr")
conflict_prefer("rename","dplyr")
conflict_prefer("desc","dplyr")
conflict_prefer("cor","stats")
WGCNA::allowWGCNAThreads(nThreads = 5)
testInteger <- function(x){
  test <- all.equal(x, as.integer(x), check.attributes = FALSE)
  if(test == TRUE){ return(TRUE) }
  else { return(FALSE) }
}
# patch function ----------------------------------------------------------

# 01. UI =========================
## logo
customLogo <- shinyDashboardLogoDIY(
  
  boldText = "ShawnLearnBioinfo"
  ,mainText = "WGCNA by click mouse"
  ,textSize = 14
  ,badgeText = "V0.1.0"
  ,badgeTextColor = "white"
  ,badgeTextSize = 2
  ,badgeBackColor = "#40E0D0"
  ,badgeBorderRadius = 3
  
)

ui <- shinyUI(
  navbarPage(theme = shinytheme("spacelab"),
             customLogo,
             tabPanel(
               useShinyjs(),
               title = "Data import and cleaning",
               icon = icon("file-upload"),
               sidebarLayout(
                 div(id = "Sidebar",
                     sidebarPanel(
                       width = 2,
                       fileInput(
                         inputId = "ExpMat",
                         label = "Upload expression matrix",
                         accept = c(".txt",".csv",".xls")
                       ),
                       p("only accecpt Tab-delimited .txt, .csv and .xls file",
                         style = "color: #7a8788;font-size: 12px; font-style:Italic"),
                       radioButtons(
                         inputId = "format",
                         label = "Format",
                         choices = c("count","expected count","normalized count","peak area (metabolomics)","protein abundance"),
                         selected = "count"
                       ),
                       p("Normalized included: DEseq2::vst (for count and expected count) \nraw data \nlog10(x + 1) (for normalized data)  \nlog(x) for metabolomics or proteomics data",
                         style = "color: #7a8788;font-size: 12px; font-style:Italic"),
                       selectInput(
                         inputId = "method1",
                         label = "Normalized method",
                         choices = c(vst = "vst",
                                     raw = "raw",
                                     logarithm = "logarithm"),
                         selected = "vst"
                       ),
                       HTML('<font color = #FF6347  size = 3.2><b>First Time filter</b></font>'),
                       textInput(
                         inputId = "SamPer",
                         label = "Sample percentage",
                         value = "0.9"
                       ),
                       textInput(
                         inputId = "RCcut",
                         label = "Expression Cutoff",
                         value = "10"
                       ),
                       p("Noise removal, for example, removing all features that have a count of less than say 10 in more than 90% of the samples",
                         style = "color: #7a8788;font-size: 12px; font-style:Italic"),
                       br(),
                       HTML('<font color = #FF6347 size = 3.2><b>Second Time filter</b></font>'),
                       radioButtons(
                         inputId = "CutMethod",
                         label = "Filter Method",
                         choices = c("MAD","Var"),
                         selected = "MAD"
                       ),
                       textInput(
                         inputId = "remain",
                         label = "Reserved genes Num.",
                         value = "8000"
                       ),
                       p("Probesets or genes may be filtered by mean expression or variance (or their robust analogs such as median and median absolute deviation, MAD) since low-expressed or non-varying genes usually represent noise. Whether it is better to filter by mean expression or variance is a matter of debate; both have advantages and disadvantages, but more importantly, they tend to filter out similar sets of genes since mean and variance are usually related.",
                         style = "color: #7a8788;font-size: 12px; font-style:Italic"),
                     )# sidebarPanel
                 ),# div
                 mainPanel(
                   fluidPage(
                     actionButton("toggleSidebar",
                                  "Toggle sidebar"),
                     actionButton("action1", "Update information!"),
                     tabsetPanel(
                       tabPanel(title = "Input file check",height = "500px",width = "100%",
                                icon = icon("check-circle"),
                                htmlOutput("Inputcheck"),
                                htmlOutput("filter1"),
                                htmlOutput("filter2"),
                       ),
                       tabPanel(title = "Preview of Input",height = "500px",width = "100%",
                                icon = icon("table"),
                                DT::dataTableOutput("Inputbl"),
                       ),
                       tabPanel(title = "SampleCluster",height = "500px",width = "100%",
                                icon = icon("tree"),
                                jqui_resizable(
                                  plotOutput("clustPlot")
                                ),
                                downloadButton("downfig1","Download")
                       ),# tabPanel
                       tabPanel(title = "RemoveOutlier (option)",height = "500px",width = "100%",
                                icon = icon("trash"),
                                selectInput(
                                  inputId = "outlier",
                                  label = "Select samples",
                                  choices = c("select a sample","sample2","..."),
                                  selected = "select a sample",
                                  multiple = T
                                ),
                                actionButton("Startremove","remove outlier"),
                                downloadButton("downtabout","Download"),
                                DT::dataTableOutput("new_mat_preview")
                       )
                     )
                     
                   )# fluidPage
                 )#mainPanel
               ) # sidebarLayout
             ),##tabPanel
             tabPanel(
               useShinyjs(),
               title = "SFT and Power Select",
               icon = icon("play-circle"),
               sidebarLayout(
                 div(id = "Sidebar2",
                     sidebarPanel(
                       width = 2,
                       sliderInput(
                         inputId = "CutoffR",
                         label = HTML('R<sup>2</sup> cutoff'),
                         min = 0,
                         max = 1,
                         value = 0.8
                       ),
                       br(),
                       HTML('<font size = 2.5 color = #7a8788><i>WGCNA will generate a recommended power value. If it does not match, a power will be given according to the experience list in the WGCNA FAQ. I don’t like this experience power very much. <font color = blue>If you find that the R <sup>2</sup> value corresponding to experience power given by the software lower than your setting Threshold </font>,<font color = purple><b> please select a customized power based on the SFT plot.</b></i></font></font>'),
                       radioButtons(
                         inputId = "PowerTorF",
                         label = "Power type",
                         choices = c("Recommended","Customized"),
                         selected = "Recommended"
                       ),
                       
                       sliderInput(
                         inputId = "PowerSelect",
                         label = "Final Power Selection",
                         min = 1,
                         max = 33,
                         value = 6
                       )
                     )
                 ),
                 mainPanel(
                   fluidPage(
                     #### output field
                     actionButton("toggleSidebar2",
                                  "Toggle sidebar"),
                     tabsetPanel(
                       tabPanel(title = "Select Power",height = "500px",width = "100%",
                                icon = icon("th"),
                                actionButton("Startsft","Start analysis"),
                                htmlOutput("powerout"),
                                jqui_resizable(
                                  plotOutput("sftplot")
                                ),
                                textInput(inputId = "width2",
                                          label = "width",
                                          value = 10),
                                textInput(inputId = "height2",
                                          label = "height",
                                          value = 10),
                                actionButton("adjust2","Set fig size"),
                                downloadButton("downfig2","Download")
                                
                                
                       ),
                       tabPanel(title = "Information of sft table",height = "500px",width = "100%",
                                icon = icon("table"),
                                DT::dataTableOutput("sfttbl")
                       ),
                       tabPanel(title = "scale free estimate",height = "500px",width = "100%",
                                icon = icon("chart-bar"),
                                actionButton("Startcheck","Check Scale-free"),
                                jqui_resizable(
                                  plotOutput("sfttest")
                                ),
                                textInput(inputId = "width3",
                                          label = "width",
                                          value = 10),
                                textInput(inputId = "height3",
                                          label = "height",
                                          value = 10),
                                actionButton("adjust3","Set fig size"),
                                downloadButton("downfig3","Download")
                                
                       )## tabPanel
                     )## tabsetPanel
                   )## fluidPage
                 )
               )
             ),##tabPanel
             tabPanel(
               useShinyjs(),
               title = "Module-net",
               icon = icon("play-circle"),
               sidebarLayout(
                 div(id = "Sidebar3",
                     sidebarPanel(
                       width = 2,
                       sliderInput(
                         inputId = "minMsize",
                         label = "min Module Size",min = 0,max = 200,value = 30
                       ),
                       sliderInput(
                         inputId = "mch",
                         label = "module cuttree height",
                         min = 0, max = 1, value = 0.25
                       ),
                       textInput(inputId = "blocksize",
                                 label = "select max blocksize",
                                 value = 5000),
                       p("MaxBlockSize, The default was 5000, 4GB memory could handle 8000-10000 genes, for 16GB memory you can select at most of 24000 genes in one block, 32GB should be enough for 30000-40000. Try to keep all selected genes in one block",
                         style = "color: #7a8788;font-size: 12px; font-style:Italic"),
                     )
                 ),
                 mainPanel(
                   fluidPage(
                     actionButton("toggleSidebar3",
                                  "Toggle sidebar"),
                     tabsetPanel(
                       tabPanel(
                         title = "Cluster",height = "500px",width = "100%",
                         icon = icon("table"),
                         actionButton("Startnet","Start"),
                         jqui_resizable(
                           plotOutput("cluster")
                         ),
                         textInput(inputId = "width4",
                                   label = "width",
                                   value = 10),
                         textInput(inputId = "height4",
                                   label = "height",
                                   value = 10),
                         actionButton("adjust4","Set fig size"),
                         downloadButton("downfig4","Download"),
                         
                         br(),
                         br(),
                         tableOutput("m2num")
                       ),
                       tabPanel(
                         title = "Eigengene adjacency heatmap",height = "500px",width = "100%",
                         icon = icon("th"),
                         jqui_resizable(
                           plotOutput("eah")
                         ),
                         textInput(inputId = "width5",
                                   label = "width",
                                   value = 10),
                         textInput(inputId = "height5",
                                   label = "height",
                                   value = 10),
                         actionButton("adjust5","Set fig size"),
                         downloadButton("downfig5","Download")
                       ),
                       tabPanel(
                         title = "Gene to module",height = "500px",width = "100%",
                         icon = icon("table"),
                         DT::dataTableOutput("g2m"),
                         downloadButton("downtbl2","download")
                       )
                     )
                     
                   )
                 )
               )
             ),##tabPanel
             tabPanel(
               useShinyjs(),
               title = "Module-trait",
               icon = icon("star-of-david"),
               sidebarLayout(
                 div(id = "Sidebar4",
                     sidebarPanel(
                       width = 2,
                       fileInput(
                         inputId = "traitData",
                         label = "Upload trait data",
                         accept = c(".txt",".csv",".xls")
                       ),
                       colourpicker::colourInput(inputId = "colormin",
                                                 label = "Minimum",
                                                 value = "purple"),
                       colourpicker::colourInput(inputId = "colormid",
                                                 label = "Middle",
                                                 value = "white"),
                       colourpicker::colourInput(inputId = "colormax",
                                                 label = "Maxmum",
                                                 value = "yellow"),
                       textInput(
                         inputId = "xangle",
                         label = "x axis label angle",
                         value = 0
                       ),
                       actionButton("starttrait","Start analysis"),
                       
                       hr(style = 'border-top: dotted 2px green;'),
                       h4("Iterative WGCNA (option)"),
                       p('Skip this step if you do not interested.'),
                       sliderInput(
                         inputId = "kme_cutoff",
                         label = "KME cutoff",
                         min = 0,max = 1,value = 0.8
                       ),
                       radioButtons(
                         inputId = "inter_method",
                         label = "Choose method",
                         choices = c(1,2),
                         selected = 1
                       ),
                       h5("Method1"),
                       p("For iterative WGCNA expression geneset output,It is recommended that the threshold value of KME ≥ 0.8 and not too low",style = "color: #7a8788;font-size: 12px; font-style:Italic"),
                       a("Please read and cite the paper carefully, and manually execute the iteration WGCNA",href="https://www.biorxiv.org/content/10.1101/234062v1","iterativeWGCNA: iterative refinement to improve module detection from WGCNA co-expression networks",style = "color: blue;font-size: 12px; font-style:Italic"),
                       br(),
                       img(src = "https://shawnmagic-1257599720.cos.ap-chengdu.myqcloud.com/picgo/202206081007161.png",width = "100%"),
                       h5('Method2'),
                       p("A simple and fast strategy to purify the initial co-expression network, It is recommended that the threshold value of KME ≤ 0.5,We expect that the genes with low connectivity in the final co-expression network will also be preserved, which can better help us enrich biological issue",style = "color: #7a8788;font-size: 12px; font-style:Italic"),
                       br(),
                       img(src = "https://shawnmagic-1257599720.cos.ap-chengdu.myqcloud.com/picgo/202206081018077.png",width = "100%")
                     )
                 ),
                 mainPanel(
                   actionButton("toggleSidebar4",
                                "Toggle sidebar"),
                   fluidPage(
                     tabsetPanel(
                       tabPanel(
                         title = "Module to trait",height = "500px",width = "100%",
                         icon = icon("ht"),
                         jqui_resizable(
                           plotOutput("mtplot")
                         ),
                         textInput(inputId = "width6",
                                   label = "width",
                                   value = 10),
                         textInput(inputId = "height6",
                                   label = "height",
                                   value = 10),
                         actionButton("adjust6","Set fig size"),
                         downloadButton("downfig6","Download")
                       ),
                       tabPanel(
                         title = "Module-trait matrix",height = "500px",width = "100%",
                         icon = icon("table"),
                         DT::dataTableOutput("traitmat"),
                         DT::dataTableOutput("traitp")
                       ),
                       tabPanel(
                         title = "Eigengene-based connectivities,KME",height = "500px",width = "100%",
                         icon = icon("table"),
                         DT::dataTableOutput("KME"),
                         
                         downloadButton("downtbl3","download")
                       ),
                       tabPanel(
                         title = "Iterative WGCNA (option)",height = "500px",width = "100%",
                         icon = icon("table"),
                         actionButton(inputId = "run_filter",label = "Run KME filter",icon = icon("play")),
                         h3("Retained Genesets"),
                         DT::dataTableOutput("Iter_retained"),
                         br(),
                         h3("Removed Genesets"),
                         DT::dataTableOutput("Iter_removed"),
                         br(),
                         downloadButton("downtbl_Iter_1","download retained Geneset."),
                         downloadButton("downtbl_Iter_2","download removed Geneset")
                       ),
                       tabPanel(
                         title = "Env output (option)",
                         icon = icon("file"),
                         p("Data integration and download.The EnvImg.rds file can be used for data analysis in R.",style = "color: #7a8788;font-size: 16px; font-style:Italic"),
                         br(),
                         p("Run the following code in R.",style = "color: green;font-size: 14px; font-style:Italic"),
                         p("readRSD('filepath/EnvImg.rsd')",style = "color: blue;font-size: 14px;"),
                         actionButton("Integrate_data","Integrate Data"),
                         downloadButton("downtbl_EnvImg","download global environment.")
                       )
                     )
                   )
                 )
               )
             ),##tabPanel
             tabPanel(
               useShinyjs(),
               title = "Interested module",
               icon = icon("broom"),
               sidebarLayout(
                 div(id = "sidebar5",
                     sidebarPanel(
                       width = 2,
                       selectInput(
                         inputId = "strait",
                         label = "Select traits",
                         choices = c("select a trait","trait2","..."),
                         selected = "select a trait",
                         multiple = F
                       ),
                       selectInput(
                         inputId = "smodule",
                         label = "Select module",
                         choices = c("red","black","..."),
                         selected = "red",
                         multiple = F
                       ),
                       actionButton("InterMode","Start Analysis")
                     )
                 ),
                 mainPanel(
                   actionButton("toggleSidebar5",
                                "Toggle sidebar"),
                   fluidPage(
                     tabsetPanel(
                       tabPanel(
                         title = "GS-Connectivity",height = "500px",width = "100%",
                         icon = icon("chart-line"),
                         jqui_resizable(
                           plotOutput("GSCon")
                         ),
                         textInput(inputId = "width7",
                                   label = "width",
                                   value = 10),
                         textInput(inputId = "height7",
                                   label = "height",
                                   value = 10),
                         actionButton("adjust7","Set fig size"),
                         downloadButton("downfig7","Download")
                       ),
                       tabPanel(
                         title = "Heatmap",height = "500px",width = "100%",
                         icon = icon("buromobelexperte"),
                         jqui_resizable(
                           plotOutput("heatmap")
                         ),
                         textInput(inputId = "width8",
                                   label = "width",
                                   value = 10),
                         textInput(inputId = "height8",
                                   label = "height",
                                   value = 10),
                         actionButton("adjust8","Set fig size"),
                         downloadButton("downfig8","Download")
                       ),
                       tabPanel(
                         title = "MM vs GS all",height = "500px",width = "100%",
                         icon = icon("chart-line"),
                         jqui_resizable(
                           plotOutput("GSMM.all")
                         ),
                         textInput(inputId = "width10",
                                   label = "width",
                                   value = 10),
                         textInput(inputId = "height10",
                                   label = "height",
                                   value = 10),
                         actionButton("adjust10","Set fig size"),
                         downloadButton("downfig10","Download")
                       )
                     )
                   )
                 )
               ),
             ),##tabPanel
             tabPanel(
               useShinyjs(),
               title = "hub gene",
               icon = icon("star"),
               sidebarLayout(
                 div(id = "sidebar6",
                     sidebarPanel(
                       width = 2,
                       selectInput(
                         inputId = "hubtrait",
                         label = "Select trait",choices = c("Select a trait","..."),multiple = F
                       ),
                       selectInput(
                         inputId = "hubmodule",
                         label = "Select module",choices = c("Select a module","..."),multiple = F
                       ),
                       actionButton("starthub","Start Analysis")
                     )
                 ),
                 mainPanel(
                   actionButton("toggleSidebar6",
                                "Toggle sidebar"),
                   fluidPage(
                     tabsetPanel(
                       tabPanel(
                         title = " choose Top Hub In Each Module (Not recommended)",
                         icon = icon("sad-cry"),
                         DT::dataTableOutput("cthub")
                       ),
                       tabPanel(
                         title = "By kME and GS (Yes!)",
                         icon = icon("smile"),
                         sliderInput(
                           inputId = "kMEcut",
                           label = "cutoff of  absolute value of kME",
                           min = 0,max = 1,step = 0.01,
                           value = 0.5
                         ),
                         sliderInput(
                           inputId = "GScut",
                           label = "cutoff of  absolute value of GS",
                           min = 0,max = 1,step = 0.01,
                           value = 0.5
                         ),
                         DT::dataTableOutput("kMEhub"),
                         downloadButton("downtbl4","download")
                       ),
                       tabPanel(
                         title = "Cytoscape output",
                         icon = icon("dna"),
                         textInput(
                           inputId = "threshold",
                           label = "weight threshold",
                           value = 0.02
                         ),
                         actionButton("threadd","choose the threshold"),
                         DT::dataTableOutput("edgeFile"),
                         DT::dataTableOutput("nodeFile"),
                         downloadButton("downtbl5","download edgefile"),
                         downloadButton("downtbl6","download nodefile")
                       )
                     )
                   )
                 )
               )
             )##tabPanel
             
  )## navbarPage
)## UI

server <- function(input, output, session){
  observeEvent(input$toggleSidebar, {
    shinyjs::toggle(id = "Sidebar")
  })
  observeEvent(input$toggleSidebar2, {
    shinyjs::toggle(id = "Sidebar2")
  })
  observeEvent(input$toggleSidebar3, {
    shinyjs::toggle(id = "Sidebar3")
  })
  observeEvent(input$toggleSidebar4, {
    shinyjs::toggle(id = "Sidebar4")
  })
  observeEvent(input$toggleSidebar5, {
    shinyjs::toggle(id = "Sidebar5")
  })
  observeEvent(input$toggleSidebar6, {
    shinyjs::toggle(id = "Sidebar6")
  })
  data <- reactive({
    file1 <- input$ExpMat
    if(is.null(file1)){return()}
    read.delim(file = file1$datapath,
               sep="\t",
               header = T,
               stringsAsFactors = F)
  })
  data_check = reactive(
    {
      if (isTRUE(testInteger(data()[,2]))) {
        "count"
      } else {
        "non-count"
      }
    }
  )
  
  fmt_select = reactive({
    if (isTRUE(testInteger(data()[,2]))) {
      "count"
    } else {
      "normalized count, peak area (metabolomics), protein abundance or expected count"
    }
  })
  output$Inputcheck = renderUI({
    if(is.null(data())){return()}
    if(length(which(is.na(data()))) == 0) {
      HTML('<font color = red><b>
          Congratulations!,</b></font> There is no problem with your expression matrix format, please proceed to the next step','<br/>',
          '<font color = red> <b>Notice:</font> </b> It seems your input data is a:','<font color = red><b>',data_check(),'</b></font> expression matrix, it is recommend to select <font color = red> <b>',fmt_select(),'</b></font> in <font color = blue> <b>Format </font> </b>','<br/>',
          '<font color = red> <b>Notice:</font> </b> If readcount is <font color = red> <b> expected count</font> </b> generated by RSEM, Please select <font color = red> <b> expected count </font> </b>.')
    } else {
      HTML(
        '<font color = blue><b>Sorry!</b></font>
       Your expression matrix has blank (NA) values or blank (NA) rows,<font color = blue> Please double check and manually remove the blanks or rows and upload file again</font>
       '
      )
    }
    
  })
  ## count number
  fmt = reactive({
    input$format
  })
  observe({
    if(fmt() == "count" | fmt() == "expected count") {
      updateSelectInput(session, "method1",choices = c(vst = "vst"))
      updateTextInput(session,"RCcut",value = 10)
    } else {
      updateSelectInput(session, "method1",choices = c(raw = "raw",
                                                       logarithm = "logarithm"))
      updateTextInput(session,"RCcut",value = 1)
    }
    
  })
  mtd = reactive({
    input$method1
  })
  
  sampP = reactive({
    as.numeric(input$SamPer)
  })
  rccutoff = reactive({
    as.numeric(input$RCcut)
  })
  GNC = reactive({
    as.numeric(input$remain)
  })
  cutmethod = reactive({
    input$CutMethod
  })
  ## set reactiveValues
  exp.ds<-reactiveValues(data=NULL)
  downloads <- reactiveValues(data = NULL)
  observeEvent(
    input$action1,
    {
      if(is.null(data())){return()}
      if(length(which(is.na(data()))) != 0) {return()}
      exp.ds$table = data.frame()
      exp.ds$table2 = data.frame()
      exp.ds$param = list()
      exp.ds$gnccheck = list()
      exp.ds$layout = as.character(input$treelayout)
      exp.ds$GNC_check = list()
      exp.ds$GNC = list()
      output$filter1 = renderUI({
        input$action1
        p_mass = c("Processing step1, remove very low expressed genes",
                   paste("Processing step2, pick out high variation genes via",cutmethod()))
        withProgress(
          message = "Raw data normlization",
          value = 0,{
            for (i in 1:2) {
              incProgress(1/2,detail = p_mass[i])
              if(i == 1) {
                exp.ds$table = getdatExpr(rawdata = data(),
                                          RcCutoff = rccutoff(),samplePerc = sampP(),
                                          datatype = fmt(),method = mtd())
              } else if (i == 2){
                exp.ds$GNC_check = GNC() - nrow(exp.ds$table)
                if(exp.ds$GNC_check > 0 ) {
                  exp.ds$GNC = nrow(exp.ds$table)
                  exp.ds$gnccheck = "The number of Genes you want to retain is greater than the total number of genes after the first filter. The number of genes retained here is equal to the total number after the first filter."
                } else {
                  exp.ds$GNC = GNC()
                  exp.ds$gnccheck = "All going well!"
                }
                exp.ds$table2 = getdatExpr2(datExpr = exp.ds$table,
                                            GeneNumCut = 1-exp.ds$GNC/nrow(exp.ds$table),cutmethod = cutmethod())
                exp.ds$param = getsampleTree(exp.ds$table2,layout = exp.ds$layout)
              }
              Sys.sleep(0.1)
            }
          }
        )
        isolate(HTML(paste0('<font color = red> <b>After filtered by conditions:</b> </font>removing all features that have a count of less than say <font color = red><b>',rccutoff(),'</b></font> in more than <font color = red> <b>',100*sampP(),'% </b></font> of the samples','<br/>',
                            '<font color = red> <b>Remaining Gene Numbers: </b> </font>',nrow(exp.ds$table),'<br/>',
                            '<font color = red> <b>After filtered by conditions:</b> </font>Genes with <font color = red><b>',cutmethod(),'</b></font> ranked top <font color = red> <b>',exp.ds$GNC,' </b></font> of all expressed genes','<br/>',
                            '<font color = red> <b>Remaining Gene Numbers: </b> </font>',ncol(exp.ds$table2),'<br/>',
                            '<font color = red> <b>Notice: </b> </font>',exp.ds$gnccheck
                            )
                     )
                )
      })
    }
  )
  
  
  ## summary num
  output$Inputbl = DT::renderDataTable({
    if(is.null(data())){return()}
    if(length(which(is.na(data()))) != 0) {return()}
    as.data.frame(t(exp.ds$table2))
  })
  ## sample tree
  output$clustPlot = renderPlot({
    if(is.null(data())){return()}
    if(length(which(is.na(data()))) != 0) {return()}
    if(is.null(exp.ds$table2)){return()}
    plot(exp.ds$param$sampleTree,main = "Sample clustering to detect outlier", sub = "", xlab = "")
  })
  
  ##remove outlier
  s_outlier = reactive({
    colnames(data())
  })
  
  observe({
    updateSelectInput(session, "outlier",choices = s_outlier())
  })
  observeEvent(
    input$Startremove,
    {
      exp.ds$outliersamples = as.character(input$outlier)
      if(!exp.ds$outliersamples[1]%in%s_outlier()) {return()} else {
        exp.ds$rmoutlier = 
          mv_outlier(x = data(),y = exp.ds$outliersamples)
      }
      output$new_mat_preview = DT::renderDataTable({
        if(is.null(exp.ds$table2)){return()}
        if(!exp.ds$outliersamples[1]%in%s_outlier()) {return()}
        as.data.frame(exp.ds$rmoutlier)
      })
    }
  )
  
  
  
  
  ## download sample tree
  
  rscut = reactive({
    as.numeric(input$CutoffR)
  })
  
  observeEvent(
    input$Startsft,
    {
      if(is.null(exp.ds$table2)){return()}
      exp.ds$sft = list()
      output$powerout = renderUI({
        sft_mess = c("pick soft threshold in processing ...",
                     "Finish.")
        withProgress(message = 'SFT selection', value = 0,
                     expr = {
                       for (i in 1:2) {
                         incProgress(1/2, detail = sft_mess[i] )
                         if (i == 1) {
                           exp.ds$sft = getpower(datExpr = exp.ds$table2,rscut = rscut())
                         } else {
                           return()
                         }
                         
                       }
                     })
        isolate(HTML(paste0('<font color = red> <b>The power recommended by WGCNA is:</b> </font><font color = bule><b>',exp.ds$sft$power,'</b></font> ','<br/>',
                            '<font color = pink> <i>If all power values lower than the R square threshold which you set, it means that the power value is an empirical value. At this time, you need to infer a power value based on the results on this plot and check whether it can form a scale-free network. </i> </font>')))
      })
    }
  )
  
  
  ## outsft
  output$sftplot = renderPlot({
    if(is.null(exp.ds$table2)){return()}
    input$Startsft
    if(length(exp.ds$sft) == 0){return()}
    exp.ds$sft$plot
  })
  ## outtbl
  output$sfttbl = DT::renderDataTable({
    if(is.null(exp.ds$table2)){return()}
    input$Startsft
    if(length(exp.ds$sft) == 0){return()}
    as.data.frame(exp.ds$sft$sft)
  })
  ## test sft
  pcus = reactive({
    as.numeric(input$PowerSelect)
  })
  PowerTorF = reactive({
    input$PowerTorF
  })
  observeEvent(
    input$Startcheck,
    {
      if(is.null(exp.ds$table2)){return()}
      if(is.null(exp.ds$sft)){return()}
      sftcheck_mess = c("Checking scale free network ...",
                        "Finish.")
      exp.ds$power = exp.ds$sft$power
      exp.ds$cksft = list()
      withProgress(message = 'SFT selection', value = 0,
                   expr = {
                     for (i in 1:2) {
                       incProgress(1/2, detail = sftcheck_mess[i] )
                       if (i == 1) {
                         if(PowerTorF() == "Recommended"){
                           exp.ds$power = exp.ds$sft$power
                           exp.ds$cksft = powertest(power.test = exp.ds$sft$power,datExpr = exp.ds$table2,nGenes = exp.ds$param$nGenes)
                         } else if (PowerTorF() == "Customized"){
                           exp.ds$power = pcus()
                           exp.ds$cksft = powertest(power.test = pcus(),datExpr = exp.ds$table2,nGenes = exp.ds$param$nGenes)
                         }
                       } else {
                         return()
                       }
                     }
                   })
    }
  )
  
  output$sfttest = renderPlot({
    if(is.null(exp.ds$sft)){return()}
    input$Startcheck
    if(length(exp.ds$cksft ) == 0){return()}
    exp.ds$cksft
  })
  mms = reactive({
    as.numeric(input$minMsize)
  })
  mch = reactive({
    as.numeric(input$mch)
  })
  blocksize = reactive({
    as.numeric(input$blocksize)
  })
  observeEvent(
    input$Startnet,
    {
      if(is.null(exp.ds$table2)){return()}
      if(is.null(exp.ds$power)){return()}
      exp.ds$netout = list()
      network_mess = c("Start module detection ...",
                     "Finish.")
      withProgress(message = 'Module detection', value = 0,
                     expr = {
                       for (i in 1:2) {
                         incProgress(1/2, detail = network_mess[i] )
                         if (i == 1) {
                           exp.ds$netout = getnetwork(datExpr = exp.ds$table2,power = exp.ds$power,
                                                      minModuleSize = mms(),mergeCutHeight = mch(),maxBlocksize = blocksize())
                         } else {
                           return()
                         }
                       }
                     })
      exp.ds$nSamples = nrow(exp.ds$table2)
      exp.ds$net = exp.ds$netout$net
      exp.ds$moduleLabels = exp.ds$netout$moduleLabels
      exp.ds$moduleColors = exp.ds$netout$moduleColors
      exp.ds$MEs_col = exp.ds$netout$MEs_col
      exp.ds$MEs = exp.ds$netout$MEs
      exp.ds$Gene2module = exp.ds$netout$Gene2module
    }
  )
  output$cluster = renderPlot({
    input$Startnet
    if(is.null(exp.ds$net)){return()}
    plotDendroAndColors(exp.ds$net$dendrograms[[1]], exp.ds$moduleColors[exp.ds$net$blockGenes[[1]]],
                        "Module colors",
                        dendroLabels = FALSE, hang = 0.03,
                        addGuide = TRUE, guideHang = 0.05)
  })
  output$m2num = renderTable({
    input$Startnet
    if(is.null(exp.ds$net)){return()}
    table(exp.ds$moduleColors)
  })
  output$eah = renderPlot({
    input$Startnet
    if(is.null(exp.ds$net)){return()}
    plotEigengeneNetworks(exp.ds$MEs_col, "Eigengene adjacency heatmap",
                          marDendro = c(3,3,2,4),
                          marHeatmap = c(3,4,2,2), plotDendrograms = T,
                          xLabelsAngle = 90)
  })
  
  output$g2m = DT::renderDataTable({
    input$Startnet
    if(is.null(exp.ds$net)){return()}
    exp.ds$Gene2module
  })
  
  phen <- reactive({
    file2 <- input$traitData
    if(is.null(file2)){return()}
    read.delim(file = file2$datapath,
               sep="\t",
               header = T,
               stringsAsFactors = F)
  })
  
  
  observeEvent(
    input$starttrait,
    {
      if(is.null(phen())){return()}
      if (ncol(phen()) == 2) {
        x <- phen()
        Tcol = as.character(unique(x[,2]))
        b <- list()
        for (i in 1:length(Tcol)) {
          b[[i]] = data.frame(row.names = x[,1],
                              levels = ifelse(x[,2] == Tcol[i],1,0))
        }
        c <- bind_cols(b)
        c <- data.frame(row.names = x$name,
                        c)
        colnames(c) = Tcol
        rownames(c) = phen()[,1]
        exp.ds$phen<- c
      } else {
        exp.ds$phen = data.frame(row.names = phen()[,1],
                                 phen()[,-1])
      }
      exp.ds$phen =  exp.ds$phen[match(rownames(exp.ds$table2),rownames(exp.ds$phen)),]
      exp.ds$traitout = list()
      exp.ds$KME = list()
      m2t2_mess = c("Module-trait relationships ...",
                    "Calculate KME",
                   "Finish.")
      withProgress(message = 'Module-trait', value = 0,
                   expr = {
                     for (i in 1:3) {
                       incProgress(1/3, detail = m2t2_mess[i] )
                       if (i == 1) {
                         exp.ds$traitout = getMt(phenotype = exp.ds$phen,
                                                 nSamples = exp.ds$nSamples,moduleColors = exp.ds$moduleColors,datExpr = exp.ds$table2)
                       } else if (i == 2){
                         exp.ds$KME = getKME(datExpr = exp.ds$table2,moduleColors = exp.ds$moduleColors,MEs_col = exp.ds$MEs_col)
                       } else {
                         return()
                       }
                     }
                   })
      exp.ds$xangle = as.numeric(input$xangle)
      exp.ds$c_min = as.character(input$colormin)
      exp.ds$c_mid = as.character(input$colormid)
      exp.ds$c_max = as.character(input$colormax)
      exp.ds$modTraitCor = exp.ds$traitout$modTraitCor
      exp.ds$modTraitP = exp.ds$traitout$modTraitP
      exp.ds$textMatrix = exp.ds$traitout$textMatrix
      exp.ds$mod_color = gsub(pattern = "^..",replacement = "",rownames(exp.ds$modTraitCor))
      exp.ds$mod_color_anno = setNames(exp.ds$mod_color,rownames(exp.ds$modTraitCor))
      exp.ds$Left_anno = rowAnnotation(
        Module = rownames(exp.ds$modTraitCor),
        col = list(
          Module = exp.ds$mod_color_anno
        ),
        show_legend = F,
        show_annotation_name = F
      )
    }
  )
  
  output$mtplot = renderPlot({
    input$starttrait
    if(is.null(phen())){return()}
    if(is.null(exp.ds$phen)){return()}
    
    exp.ds$heatmap = list()
    m2t_mess = c("Draw heatmap ...",
                     "Finish.")
    withProgress(message = 'Module-trait', value = 0,
                 expr = {
                   for (i in 1:2) {
                     incProgress(1/2, detail = m2t_mess[i] )
                     if (i == 1) {
                       exp.ds$heatmap = Heatmap(
                         matrix = exp.ds$modTraitCor,
                         cluster_rows = F, cluster_columns = F,
                         left_annotation = exp.ds$Left_anno,
                         cell_fun = function(j,i,x,y,width,height,fill) {
                           grid.text(sprintf(exp.ds$textMatrix[i,j]),x,y,gp = gpar(fontsize = 12))
                         },
                         row_names_side = "left",
                         column_names_rot = exp.ds$xangle,
                         heatmap_legend_param = list(
                           at = c(-1,-0.5,0,0.5, 1),
                           labels = c("-1","-0.5", "0","0.5", "1"),
                           title = "",
                           legend_height = unit(9, "cm"),
                           title_position = "lefttop-rot"
                         ),
                         rect_gp = gpar(col = "black", lwd = 1.2),
                         column_title = "Module-trait relationships",
                         column_title_gp = gpar(fontsize = 15, fontface = "bold"),
                         col = colorRamp2(c(-1, 0, 1), c(exp.ds$c_min, exp.ds$c_mid, exp.ds$c_max))
                       )
                     } else {
                       return()
                     }
                   }
                 })
    draw(exp.ds$heatmap)
  })
  
  output$traitmat = DT::renderDataTable({
    input$starttrait
    if(is.null(phen())){return()}
    if(is.null(exp.ds$phen)){return()}
    as.data.frame(exp.ds$modTraitCor)
  })
  
  output$traitp = DT::renderDataTable({
    input$starttrait
    if(is.null(phen())){return()}
    if(is.null(exp.ds$phen)){return()}
    as.data.frame(exp.ds$modTraitP)
  })
  
  output$KME = DT::renderDataTable({
    input$starttrait
    if(is.null(phen())){return()}
    if(is.null(exp.ds$phen)){return()}
    as.data.frame(exp.ds$KME)
  })
  
  observeEvent(
    input$run_filter,
    {
      if(is.null(phen())){return()};
      exp.ds$kme_method = as.numeric(input$inter_method);
      exp.ds$kme_cutoff = as.numeric(input$kme_cutoff);
      exp.ds$kme_outlist = iterative_out(
        g2m = exp.ds$Gene2module,
        rawMat = data(),
        tbl = as.data.frame(exp.ds$KME),
        method = exp.ds$kme_method,
        KME_cutoff = exp.ds$kme_cutoff
      )
      message("kme analysis finish")
    }
  )
  
  output$Iter_retained = DT::renderDataTable({
    input$run_filter
    if(is.null(phen())){return()}
    as.data.frame(exp.ds$kme_outlist$retain)
  })
  
  output$Iter_removed = DT::renderDataTable({
    input$run_filter
    if(is.null(phen())){return()}
    as.data.frame(exp.ds$kme_outlist$remove)
  })
  
  s_mod = reactive({
    gsub(pattern = "^..",replacement = "",rownames(exp.ds$modTraitP))
  })
  
  observe({
    updateSelectInput(session, "smodule",choices = s_mod())
  })
  
  s_trait = reactive({
    colnames(exp.ds$modTraitP)
  })
  observe({
    updateSelectInput(session, "strait",choices = s_trait())
  })
  
  observeEvent(
    input$InterMode,
    {
      if(is.null(phen())){return()}
      if(is.null(exp.ds$phen)){return()}
      exp.ds$GSout = getMM(datExpr = exp.ds$table2,MEs_col = exp.ds$MEs_col,nSamples = exp.ds$nSamples,corType = "pearson")
      exp.ds$MM = exp.ds$GSout$MM
      exp.ds$MMP = exp.ds$GSout$MMP
      exp.ds$sml = as.character(input$smodule)
      exp.ds$st = as.character(input$strait)
      exp.ds$Heatmap = moduleheatmap(datExpr = exp.ds$table2,MEs = exp.ds$MEs_col,which.module = exp.ds$sml,
                                     moduleColors = exp.ds$moduleColors)
    }
  )
  
  output$GSCon = renderPlot({
    input$InterMode
    if(is.null(exp.ds$st)){return()}
    if(is.null(exp.ds$sml)){return()}
    getverboseplot(datExpr = exp.ds$table2,module = exp.ds$sml,pheno = exp.ds$st,MEs = exp.ds$MEs_col,
                   traitData = exp.ds$phen,moduleColors = exp.ds$moduleColors,
                   geneModuleMembership = exp.ds$MM,nSamples = exp.ds$nSamples)
  })
  
  output$heatmap = renderPlot({
    input$InterMode
    if(is.null(exp.ds$st)){return()}
    if(is.null(exp.ds$sml)){return()}
    exp.ds$Heatmap
  })
  
  output$GSMM.all = renderPlot({
    input$InterMode
    if(is.null(exp.ds$st)){return()}
    if(is.null(exp.ds$sml)){return()}
    MMvsGSall(which.trait = exp.ds$st,
              traitData = exp.ds$phen,
              datExpr = exp.ds$table2,
              moduleColors = exp.ds$moduleColors,
              geneModuleMembership = exp.ds$MM,
              MEs = exp.ds$MEs_col,
              nSamples = exp.ds$nSamples)
  })
  
  observe({
    updateSelectInput(session, "hubmodule",choices = s_mod())
  })
  
  observe({
    updateSelectInput(session, "hubtrait",choices = s_trait())
  })
  
  observeEvent(
    input$starthub,
    {
      exp.ds$hubml = as.character(input$hubmodule)
      exp.ds$hubt = as.character(input$hubtrait)
      exp.ds$kMEcut = as.numeric(input$kMEcut)
      exp.ds$GScut = as.numeric(input$GScut)
      print(exp.ds$hubml)
      exp.ds$hub.all = hubgenes(datExpr = exp.ds$table2,
                                mdl = exp.ds$hubml,
                                power = exp.ds$power,
                                trt = exp.ds$hubt,
                                KME = exp.ds$KME,
                                GS.cut = exp.ds$GScut,
                                kME.cut =exp.ds$kMEcut,
                                datTrait = exp.ds$phen,
                                g2m = exp.ds$Gene2module
      )
      
    }
  )
  
  observeEvent(
    input$threadd,
    {
      exp.ds$threshold = as.numeric(input$threshold)
      exp.ds$cyt = cytoscapeout(datExpr = exp.ds$table2,
                                power = exp.ds$power,module = exp.ds$hubml,
                                moduleColors = exp.ds$moduleColors,
                                threshold = exp.ds$threshold)
    }
  )
  # checkAdjMat
  output$cthub = DT::renderDataTable({
    input$starthub
    if(is.null(exp.ds$hubml)){return()}
    if(is.null(exp.ds$hubt)){return()}
    exp.ds$hub.all$hub1
  })
  
  output$kMEhub = DT::renderDataTable({
    input$starthub
    if(is.null(exp.ds$hubml)){return()}
    if(is.null(exp.ds$hubt)){return()}
    if(is.null(exp.ds$kMEcut)){return()}
    if(is.null(exp.ds$GScut)){return()}
    exp.ds$hub.all$hub3
  })
  
  output$edgeFile = DT::renderDataTable({
    input$threadd
    if(is.null(exp.ds$hubml)){return()}
    if(is.null(exp.ds$threshold)){return()}
    exp.ds$cyt[[1]]
  })
  
  output$nodeFile = DT::renderDataTable({
    input$threadd
    if(is.null(exp.ds$hubml)){return()}
    if(is.null(exp.ds$threshold)){return()}
    exp.ds$cyt[[2]]
  })
  
  observeEvent(
    input$Integrate_data,
    {
      exp.ds$data_interagrate = list(
        expmat = exp.ds$table2,
        traitmat = exp.ds$phen,
        expmat_format = as.character(fmt()),
        method = as.character(mtd()),
        samplePercentage = as.numeric(sampP()),
        rccutoff = as.numeric(rccutoff()),
        GNC = as.numeric(GNC()),
        cutmethod = as.character(cutmethod()),
        power = exp.ds$power,
        min_module_size = as.numeric(mms()),
        module_cuttree_height = as.numeric(mch()),
        blocksize = as.numeric(blocksize()),
        net = exp.ds$net,
        moduleLabels = exp.ds$moduleLabels,
        moduleCoolors = exp.ds$moduleColors,
        MEs = exp.ds$MEs,
        MEs_col = exp.ds$MEs_col,
        exp.ds$Gene2module,
        modTraitCor = exp.ds$modTraitCor,
        modTraitP = exp.ds$modTraitP,
        textMatrix = exp.ds$textMatrix,
        KME = exp.ds$KME
      )
      outrsd <<- exp.ds$data_interagrate
    }
  )
  
  # download ----------------------------------------------------------------
  
  
  observeEvent(
    input$adjust1,
    {
      downloads$width1 <- as.numeric(input$width1)
      downloads$height1 <-  as.numeric(input$height1)
    }
  )
  observeEvent(
    input$adjust2,
    {
      downloads$width2 <- as.numeric(input$width2)
      downloads$height2 <-  as.numeric(input$height2)
    }
  )
  observeEvent(
    input$adjust3,
    {
      downloads$width3 <- as.numeric(input$width3)
      downloads$height3 <-  as.numeric(input$height3)
    }
  )
  observeEvent(
    input$adjust4,
    {
      downloads$width4 <- as.numeric(input$width4)
      downloads$height4 <- as.numeric(input$height4)
    }
  )
  observeEvent(
    input$adjust5,
    {
      downloads$width5 <- as.numeric(input$width5)
      downloads$height5 <-  as.numeric(input$height5)
    }
  )
  observeEvent(
    input$adjust6,
    {
      downloads$width6 <- as.numeric(input$width6)
      downloads$height6 <-  as.numeric(input$height6)
    }
  )
  observeEvent(
    input$adjust7,
    {
      downloads$width7 <- as.numeric(input$width7)
      downloads$height7 <-  as.numeric(input$height7)
    }
  )
  observeEvent(
    input$adjust8,
    {
      downloads$width8 <- as.numeric(input$width8)
      downloads$height8 <-  as.numeric(input$height8)
    }
  )
  observeEvent(
    input$adjust10,
    {
      downloads$width10 <- as.numeric(input$width10)
      downloads$height10 <-  as.numeric(input$height10)
    }
  )

  output$downfig1 = downloadHandler(
    filename = function() {
      "01.SampleCluster.nwk"
    },
    content = function(file) {
      write.tree(phy = exp.ds$param$tree,file = file)
    }
  )
  output$downfig2 = downloadHandler(
    filename = function() {
      "02.SftResult.pdf"
    },
    content = function(file) {
      ggsave(plot = exp.ds$sft$plot,filename = file,width = downloads$width2,height = downloads$height2)
    }
  )
  output$downfig3 = downloadHandler(
    filename = function() {
      "03.CheckSft.pdf"
    },
    content = function(file) {
      ggsave(plot = exp.ds$cksft,filename = file,width = downloads$width3,height = downloads$height3)
    }
  )
  output$downfig4 = downloadHandler(
    filename = function() {
      "04.ClusterDendrogram.pdf"
    },
    content = function(file) {
      pdf(file = file,width = downloads$width4, height = downloads$height4 )
      plotDendroAndColors(exp.ds$net$dendrograms[[1]], exp.ds$moduleColors[exp.ds$net$blockGenes[[1]]],
                          "Module colors",
                          dendroLabels = FALSE, hang = 0.03,
                          addGuide = TRUE, guideHang = 0.05)
      dev.off()
    }
  )
  output$downfig5 = downloadHandler(
    filename = function() {
      "05.EigengeneadJacencyHeatmap.pdf"
    },
    content = function(file) {
      pdf(file = file,width = downloads$width5, height = downloads$height5)
      plotEigengeneNetworks(exp.ds$MEs_col, "Eigengene adjacency heatmap",
                            marDendro = c(3,3,2,4),
                            marHeatmap = c(3,4,2,2), plotDendrograms = T,
                            xLabelsAngle = 90)
      dev.off()
    }
  )
  output$downfig6 = downloadHandler(
    filename = function() {
      "06.Module2Trait.pdf"
    },
    content = function(file) {
      pdf(file = file,width = downloads$width6, height = downloads$height6)
      
      print(Heatmap(
        matrix = exp.ds$modTraitCor,
        cluster_rows = F, cluster_columns = F,
        left_annotation = exp.ds$Left_anno,
        cell_fun = function(j,i,x,y,width,height,fill) {
          grid.text(sprintf(exp.ds$textMatrix[i,j]),x,y,gp = gpar(fontsize = 12))
        },
        row_names_side = "left",
        column_names_rot = exp.ds$xangle,
        heatmap_legend_param = list(
          at = c(-1,-0.5,0,0.5, 1),
          labels = c("-1","-0.5", "0","0.5", "1"),
          title = "",
          legend_height = unit(9, "cm"),
          title_position = "lefttop-rot"
        ),
        rect_gp = gpar(col = "black", lwd = 1.2),
        column_title = "Module-trait relationships",
        column_title_gp = gpar(fontsize = 15, fontface = "bold"),
        col = colorRamp2(c(-1, 0, 1), c(exp.ds$c_min, exp.ds$c_mid, exp.ds$c_max))
      ))
      
      dev.off()
    }
  )
  output$downfig7 = downloadHandler(
    filename = function() {
      paste0("07.GS",exp.ds$sml,"-",exp.ds$st,"-Connectivity.pdf")
    },
    content = function(file) {
      pdf(file = file,width = downloads$width7, height = downloads$height7)
      print(getverboseplot(datExpr = exp.ds$table2,module = exp.ds$sml,pheno = exp.ds$st,MEs = exp.ds$MEs_col,
                           traitData = exp.ds$phen,moduleColors = exp.ds$moduleColors,
                           geneModuleMembership = exp.ds$MM,nSamples = exp.ds$nSamples))
      dev.off()
    }
  )
  output$downfig8 = downloadHandler(
    
    filename = function() {
      paste0("08.",exp.ds$sml,"-",exp.ds$st,"MEandGeneHeatmap.pdf")
    },
    content = function(file) {
      pdf(file = file,width = downloads$width8, height = downloads$height8)
      print(exp.ds$Heatmap)
      dev.off()
    }
  )
  output$downfig10 = downloadHandler(
    
    filename = function() {
      "09.GSvsMM.all.pdf"
    },
    content = function(file) {
      pdf(file = file,width = downloads$width8, height = downloads$height8)
      print(MMvsGSall(which.trait = exp.ds$st,
                      traitData = exp.ds$phen,nSamples = exp.ds$nSamples,
                      datExpr = exp.ds$table2,
                      moduleColors = exp.ds$moduleColors,
                      geneModuleMembership = exp.ds$MM,MEs = exp.ds$MEs_col))
      dev.off()
    }
  )
  output$downtbl2 = downloadHandler(
    
    filename = function() {
      if(is.null(exp.ds$net)){return()}
      "01.Gene2Module.xls"
    },
    content = function(file) {
      write.table(x = exp.ds$Gene2module,file = file,sep = "\t",row.names = F,quote = F)
    }
  )
  output$downtbl3 = downloadHandler(
    filename = function() {
      "02.KMEofAllGenes.xls"
    },
    content = function(file) {
      write.table(x = exp.ds$KME,file = file,sep = "\t",row.names = T,quote = F)
    }
  )
  output$downtbl4 = downloadHandler(
    filename = function() {
      paste0("03.",exp.ds$hubml,"-",exp.ds$hubt,"hubgene_by_GS_MM.xls")
    },
    content = function(file) {
      write.table(x = exp.ds$hub.all$hub3,file = file,sep = "\t",row.names = F,quote = F)
    }
  )
  output$downtbl5 = downloadHandler(
    filename = function() {
      paste0("04.",exp.ds$hubml,".edge.xls")
    },
    content = function(file) {
      write.table(x = exp.ds$cyt[[1]],file = file,sep = "\t",row.names = F,quote = F)
    }
  )
  output$downtbl6 = downloadHandler(
    filename = function() {
      paste0("04.cyt",exp.ds$hubml,".node.xls")
    },
    content = function(file) {
      write.table(x = exp.ds$cyt[[2]],file = file,sep = "\t",row.names = F,quote = F)
    }
  )
  output$downtabout = downloadHandler(
    filename = function() {
      paste0("00.Remove_outlier_Table.xls")
    },
    content = function(file) {
      write.table(x = exp.ds$rmoutlier,file = file,sep = "\t",row.names = F,quote = F)
    }
  )
  output$downtbl_Iter_1 = downloadHandler(
    filename = function() {
      paste0("00.Retained_GeneSet_for_Next_Round.xls")
    },
    content = function(file) {
      write.table(x = as.data.frame(exp.ds$kme_outlist$retain),file = file,sep = "\t",row.names = F,quote = F)
    }
  )
  output$downtbl_Iter_2 = downloadHandler(
    filename = function() {
      paste0("00.Removed_GeneSet.xls")
    },
    content = function(file) {
      write.table(x = as.data.frame(exp.ds$kme_outlist$remove),file = file,sep = "\t",row.names = F,quote = F)
    }
  )
  output$downtbl_EnvImg = downloadHandler(
    filename = function() {
      paste0("EnvImg.rds")
    },
    content = function(file) {
      saveRDS(outrsd,file = file)
    }
  )
}


shinyApp(ui,server,
         options = list(launch.browser = TRUE)
         )

