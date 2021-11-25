# install these three packages from Bioconductor. 
#************************************************
library(limma)
library(DelayedMatrixStats)
library(vsn)
##********************************

## install.packages() 
#*********************************
library(gtools)
library(naturalsort)
library("shinyWidgets")
library(data.table)
library(plyr)
library(gplots)
library(shinyHeatmaply)
library(matrixStats)
library(shinyShortcut)
library(heatmaply)
library(shiny)
library(plotly)
library(viridis)
library(jsonlite)
library(RColorBrewer)
library(readxl)
library(DT)
library("Hmisc")
library(xtable)
library(htmltools)
library(htmlwidgets)
library(dplyr)
library(tidyverse)
library(hrbrthemes)
library(viridis)
library(shinydashboard)
##***********************************

## install from local zipped files. 
#***********************************
library(ber) 
##**********************************

#Create a shortcut for the application
shinyShortcut(shinyDirectory = getwd(), OS = .Platform$OS.type, gitIgnore = FALSE)


# Define UI for dataset viewer app ----
ui <- fluidPage(
  
  # App title ----
  h1(id="big-heading", "iPADshiny"),
  tags$style(HTML("#big-heading{color: red;}")),  
  h3(id="small-heading", "integrated Protein Array Data management,analysis and visualization tools"),
  tags$style(HTML("#small-heading{color: blue;}")),
  
  ## titlePanel("iPADShiny:integrated Protein Array Data management,analysis and visualization tools"),
  
  sidebarLayout(
    sidebarPanel(
      # conditionalPanel(condition="input.tabs1=='DEapp'"
      #                  
      # ),
      conditionalPanel(condition="input.tabs1=='About'",
                       h4("About iPADShiny"),
                       h4("version 1.0")
                       
      ),
      conditionalPanel(condition="input.tabs1=='Data Upload'",
                       
                       fileInput("file2","Upload layout file",multiple =FALSE, accept=c('text/gpr','text/comma-separated-values,text/plain','.txt')),
                       radioButtons("detection_type", "Select a detection type:", list("Single Channel"= 1,"Double Channel"=2), selected = 1),
                       checkboxGroupInput("Ig_pro", "Select a protein type:",
                                          c("IgG" = "IgG",
                                            "IgM" = "IgM",
                                            "IgA" = "IgA",
                                            "IgE" = "IgE"),"IgG"),
                       radioButtons("pro_type", "Select a data type:", list("cy3(532)"= 532,"cy5(635)"= 635), selected = 532),
                       fileInput("file1","Upload GPR files",multiple = TRUE, accept=c('text/gpr','text/comma-separated-values,text/plain','.gpr'))
      ),
      conditionalPanel(condition="input.tabs1=='Data Preprocessing'",
                       radioButtons("measure_type", "Measuring statistic", list("Medium"=1,"Mean"=2,"Average"=3), selected = 1),
                       radioButtons("protein_type","Background Correction", list("substraction"=0,"division"=1, "minimum"=2, "normexp"=3, "neighborhood"=4), selected = 0),
                       radioButtons("qc_type", "Qualitity Control", list("None"=0,"YES"=1), selected = 1),
                       radioButtons("batch_type", "Batch Correction", list("alone"=0,"ANOVA"=1, "ComBat_p"=2, "ComBat_np"=3),selected = 0 ),
                       radioButtons("norm_type", "Normalization", list("Alone"=0,"RLM"=1,"BCF"=2,"Scaling"=3,"Quantile"=4,"LOWESS"=5,"VSN"=6),selected = 0 ),
                       # added slide bar for client to filter values based on SNR%
                       # filter will do job on current matrix but all based on the SNR3%
                       sliderInput("decimal_percent",
                                   "Filter by SNR",
                                   min = 0.000,
                                   max = 100,
                                   step = 3,
                                   value = 0),
                       # whether client pick which type of matrix it will work
                       radioButtons("data_type", "Select a data type", list("NSI"= "Signal","SNR"= "SNR", "Abs"= "Ab_score"),"Signal"),
                       actionButton("show", "Generate a Matrix/Clean", icon("paper-plane"), 
                                    style="color: #fff; background-color: #362BA8; border-color: #2e6da4; height:40px ;width:200px;"),
                       br(),
                       tags$b("Rate of Progress"),
                       progressBar(
                         id = "pb2",
                         value = 0,
                         total = 100,
                         title = "",
                         display_pct = TRUE
                       ),
                       radioButtons("preprocessing_download_type", "Choose the download type:", list(".csv"=0, ".xlsx"=1, ".txt"=2), selected = 0),
                       HTML('<p><strong> Click to download the matrix:</strong></p>'),
                       downloadBttn(
                         outputId = "preprocessing_data_download",
                         style = "unite",
                         color = "warning",
                         size = "sm"
                       )
      ),
      
      conditionalPanel(condition="input.tabs1=='QC reports'",
                       ##h4("Positive or negative controls"),
                       radioButtons("graph_type", "Which Graph do you want to generate?", list("Ig Controls (Line Chart)"=1,"Anti-Ig Controls (Line Chart)"=2, "PBS Controls(Boxplot)"=3), selected = 1),
                       actionButton("qc_graph", "Generate", icon("paper-plane"), 
                                    style="color: #fff; background-color: #362BA8; border-color: #2e6da4; height:33px ;width:120px;"),
                       br(),br(),
                       conditionalPanel('input.graph_type<3',
                                        radioButtons("legend_position", "Change your legend status:", list("Show"=0, "Hide"=1), selected = 0)
                       ),
                       
                       
                       #radioButtons("line_chart_download_type", "Choose the download type:", list(".png"='.png', ".jpeg"='.jpeg'), '.png'),
                       HTML('<p><strong> Click to download the plot:</strong></p>'),
                       downloadBttn(
                         outputId = "download_control_graph",
                         style = "unite",
                         color = "warning",
                         size = "sm"
                       ),
                       hr(),
                       ## h4("Antibody or Sample Summary"),
                       radioButtons("summary_table_type", "Table Summary is corresponding to:", list("Samples"=0, "Antigens"=1), selected = 0),
                       actionButton("generate_table_sum", "Generate", icon("paper-plane"), 
                                    style="color: #fff; background-color: #006400; border-color: #2e6da4; height:33px ;width:120px;"),
                       br(),br(),
                       radioButtons("antibodies_table_summary_download_type", "Choose the download type:", list(".csv"=0, ".xlsx"=1, ".txt"=2), selected = 0),
                       downloadBttn(
                         outputId = "antibodies_table_summary_data",
                         style = "unite",
                         color = "warning",
                         size = "sm"
                       )
      ),
      
      conditionalPanel(condition="input.tabs1=='Heatmap'",
                       width=4,
                       h4('Data Selection'),
                       fileInput(inputId="mydata", label = "Import Data(.csv)",multiple = T),
                       hr(),
                       uiOutput('data'),
                       #uiOutput('show2'),
                       br(),
                       checkboxInput('showSample','Subset Data'),
                       conditionalPanel('input.showSample',uiOutput('sample')),
                       hr(),h4('Data Preprocessing'),
                       column(width=4,selectizeInput('transpose','Transpose',choices = c('No'=FALSE,'Yes'=TRUE),selected = FALSE)),
                       column(width=8,selectizeInput("transform_fun", "Transform", c("Identity"=".", "Divide Antibodies by SD"='sd',"Mean Center Antibodies"='mean',"Median Center Antibodies"='median',"Normalize Antibodies"='normalize',Percentize='percentize',"Missing values"='is.na10', Correlation='cor'),selected = '.')),
                       uiOutput('annoVars'),
                       column(width=8,selectizeInput("row_or_col","Row or Column",choices = c("Row"=1,"Column"=2),selected = 1)),
                       br(),br(),hr(),h4('Row dendrogram'),
                       column(width=6,selectizeInput("distFun_row", "Distance method", c(Euclidean="euclidean",Maximum='maximum',Manhattan='manhattan',Canberra='canberra',Binary='binary',Minkowski='minkowski'),selected = 'euclidean')),
                       column(width=6,selectizeInput("hclustFun_row", "Clustering linkage", c(Complete= "complete",Single= "single",Average= "average",Mcquitty= "mcquitty",Median= "median",Centroid= "centroid",Ward.D= "ward.D",Ward.D2= "ward.D2"),selected = 'complete')),
                       
                       column(width=12,sliderInput("r", "Number of Clusters", min = 1, max = 15, value = 2)),    
                       #column(width=4,numericInput("r", "Number of Clusters", min = 1, max = 20, value = 2, step = 1)),   
                       br(),hr(),h4('Column dendrogram'),
                       column(width=6,selectizeInput("distFun_col", "Distance method", c(Euclidean="euclidean",Maximum='maximum',Manhattan='manhattan',Canberra='canberra',Binary='binary',Minkowski='minkowski'),selected = 'euclidean')),
                       column(width=6,selectizeInput("hclustFun_col", "Clustering linkage", c(Complete= "complete",Single= "single",Average= "average",Mcquitty= "mcquitty",Median= "median",Centroid= "ward.D2",Ward.D= "ward.D",Ward.D2= "ward.D2"),selected = 'complete')),
                       column(width=12,sliderInput("c", "Number of Clusters", min = 1, max = 15, value = 2)),
                       #column(width=4,numericInput("c", "Number of Clusters", min = 1, max = 20, value = 2, step = 1)),    
                       
                       br(),hr(),  h4('Additional Parameters'),
                       
                       column(3,checkboxInput('showColor','Color')),
                       column(3,checkboxInput('showMargin','Layout')),
                       column(3,checkboxInput('showDendo','Dendrogram')),
                       hr(),
                       conditionalPanel('input.showColor==1',
                                        hr(),
                                        h4('Color Manipulation'),
                                        uiOutput('colUI'),
                                        sliderInput("ncol", "Set Number of Colors", min = 1, max = 256, value = 256),
                                        checkboxInput('colRngAuto','Auto Color Range',value = T),
                                        conditionalPanel('!input.colRngAuto',uiOutput('colRng'))
                       ),
                       
                       conditionalPanel('input.showDendo==1',
                                        hr(),
                                        h4('Dendrogram Manipulation'),
                                        selectInput('dendrogram','Dendrogram Type',choices = c("both", "row", "column", "none"),selected = 'both'),
                                        selectizeInput("seriation", "Seriation", c(OLO="OLO",GW="GW",Mean="mean",None="none"),selected = 'OLO'),
                                        sliderInput('branches_lwd','Dendrogram Branch Width',value = 0.6,min=0,max=5,step = 0.1)
                       ),             
                       
                       conditionalPanel('input.showMargin==1',
                                        hr(),
                                        h4('Widget Layout'),
                                        column(4,textInput('main','Title','')),
                                        column(4,textInput('xlab','X Title','')),
                                        column(4,textInput('ylab','Y Title','')),
                                        sliderInput('row_text_angle','Row Text Angle',value = 0,min=0,max=180),
                                        sliderInput('column_text_angle','Column Text Angle',value = 45,min=0,max=180),
                                        sliderInput("l", "Set Margin Width", min = 0, max = 200, value = 130),
                                        sliderInput("b", "Set Margin Height", min = 0, max = 200, value = 40)
                       )
      ),
      
      conditionalPanel(condition="input.tabs1=='FAQ'",
                       h3("Questions ?"),
                       p("This App is developed and maintained by Chengsong Zhu at the",
                         a("Genomics and Microarray Core Facility, ", href="https://microarray.swmed.edu/"), 
                         a("Department of Immunology, ", href="https://www.utsouthwestern.edu/education/medical-school/departments/immunology/"), 
                         " University of Texas Southwestern Medical Center."),
                       p("If you have any questions, comments, or suggestions, feel free to contact our core at chengsongzhu@gmail.com"),
                       
                       br(),br(),br(),
                       # radioButtons("question_faq", "Select the field of your question", list("About"=0,"Data Upload"=1, "Data Preprocessing"=2, "Heatmap"=3, "Other"=4), 0),
                       br(),br(),
                       
                       HTML('<h4 style="color:black"> <Strong>Thanks for using our application!</Strong></h4>'),
                       HTML('<h5 style="color:black"> <Strong>How helpful do you think of our app?</Strong></h5>'),
                       actionButton("good", "Good", icon("kiss-wink-heart"), 
                                    style="color: #fff; background-color: #EB0909; border-color: #2e6da4; height:40px ;width:110px;"),
                       actionButton("fine", "Fine", icon("grin-tongue"), 
                                    style="color: #fff; background-color: #362BA8; border-color: #2e6da4; height:40px ;width:110px;"),
                       actionButton("bad", "Bad", icon("frown"), 
                                    style="color: #fff; background-color: #167A1E; border-color: #2e6da4; height:40px ;width:110px;"),
                       
                       br(),br(),
                       textOutput('evaluation')
      ), 
      width=3 
      
    ),
    
    
    mainPanel(
      tabsetPanel(id="tabs1",
                  tabPanel("About",
                           img(src='images.png', align = "left"),
                           br(),br(),br(),br(),
                           h3("Introduction"),
                           p("This application was developed to help biologists preprocess their protein array profiling and generate new hypothesis or validate their hypothesis.
                                       We hope that you find the iPADShiny useful and we welcome suggestions for additional features by our users.
                                       We would like to thank everyone who has made constructive suggestions so far. We are continuing add new features in the News tab."),
                           p("We present integrated Protein Array Data management,analysis and visualization tools developed in R with shiny (iPADShiny), the first web tool dedicated to analysing
                                         functional protein microarray data. iPAD incorporates the strengths of PAA, while eliminating
                                         their major limitations. iPADShiny is suitable for experimental biologists who want to analyse their
                                         own data without the need to write codes. "), 
                           
                           h3("Workflow"),
                           p("Step 1: Data Upload, Upload your input data layout and raw GPR files,
				                       via  Upload your layout file and  Upload your GPR files. 
				                       a summary of your input data will be presented."),  
                           p("Step 2: QC report, This application provides complehensive quality control reports."),             
                           p("Step 3: Data preprocessing , Once the GPR files and layout file have been uploaded, data matrix will be assembled by i, background correction ii, averaging the duplicate spots after multiple data filtering  
				                     iii,batch correction optional iV, normalization process "),
                           p("Step 4: Differential analysis,Differential protein reactivity levels between study conditions will be identified and visualized via Figs and tables"), 
                           p("Step 5: Heatmap"), 
                           #h3("Reference"),
                           # p(" ")
                           
                           h3("Feedback"),
                           p("This application uses the shiny package from RStudio. This App is developed and maintained by the " ,
                             a("Microarray core facility, ", href="https://genomics-microarray.swmed.edu/"), 
                             a("Department of Immunology, ", href="https://www.utsouthwestern.edu/education/medical-school/departments/immunology/"), 
                             a("University of Texas Southwestern Medical Center. ", href="http://www.utsouthwestern.edu/"),
                             "We are actively improving and expanding our pipeline features.If you have any questions, comments, or suggestions, feel free to contact us at chengsongzhu@gmail.com "),
                           br(),br(),br(),br(),br(),br()
                           
                           
                  ),
                  
                  
                  tabPanel("Data Upload",
                           h2("Basic information"),
                           tableOutput('basicInfo')
                  ),
  
                  tabPanel("Data Preprocessing",
                           conditionalPanel(condition="input.show%2==0",
                                            h2("Build your own matrix here"),
                                            
                                            HTML('<p> <br> This application provides multiple algorithms for biologists to do data cleaning and data analyse. After selecting the options on the left sidebar, you will get the result matrix by clicking on "Generate a Matrix/Clean" button. By clicking the same button again, your previous matrix displayed will be cleaned.</p>
                                           
                                            <br><p><strong>Background Correction:</strong> This application provides five options to do background correction. If method="<strong>substract</strong>", the application will generate a matrix by using foreground data to substract background data. If method="<strong>division</strong>", the application will 
                                            generate a matrix using the ratio of foreground data to background data. If method="<strong>minimum</strong>"
                                            then any intensity which is zero or negative after background subtraction is set equal to half the minimum of the positive corrected intensities for that array. 
                                            If method="<strong>normexp</strong>" a convolution of normal and exponential distributions is fitted to the foreground intensities using the background intensities as a covariate, and the expected signal given the observed foreground becomes the corrected intensity. If method="<strong>neighborhood</strong>", the application will generate a matrix by smoothing each value with its neighborhood.</p>
                                            
                                            <br><p><strong>KOC:</strong> This option allows the user to do data cleaning. By choosing this option, the bad samples will be automatically removed, and the sample data with the flag "-100" will be automatically change to NA. All negative values will be automatically change to 0.001. Duplicate spots will be averaged. Net signal intensities will be calculated by subtracting the medium of all negative control samples (PBS) except all the positive control antigens which have been moved to the bottom of the matrix . </p>
                                            
                                            <br><p><strong>Batch Correction:</strong> The batch effects can dramatically reduce the accuracy of statistical inference in genomic experiments. We provide users with multiple options to make data modeled and remove bad effects in order to accurately measure biological variability and 
                                            to obtain correct statistical inference when performing high-throughput genomic analysis.</p>
                                            
                                            <br><p><strong>Normalization:</strong> This application provides users with five different data normalization options. 1.Not normalized 2.Robust linear model 3.Backgroung correction Factor 4.Scaling 5.Quantile 6.Lowess 7.Variance Stabalization Normalization </p>
                                            
                                            <br><p><strong>Download Option:</strong> This application allows users to download matrix generated through our app. After choosing the data type, the user could download the matrix generated by clicking "Download Your Matrix" button. we provide three matrix options, 1) NSI: Net/Normalized Signal Intensity 2) SNR: Signal to Noise Ratio, SNR=(Feature signal intensity- background intensity)/standard deviation of background signal   3) Abs: Antibody score; Abs=log2(NSI*SNR+1) 
                                            <br> </p>')
                                            
                           ),
                           
                           
                           conditionalPanel(condition="input.show%2!=0",
                                            conditionalPanel(condition="!output.fileUploaded",
                                                             br(),br(),br(),br(),br(),br(),br(),br(),
                                                             HTML('<h2 style="color:red"> The application failed to generate the matrix.</h2>'),
                                                             HTML('<p> <br><strong> <font size="4">Please check the following things:</strong></font></p>
                                                                   <p><font size="4">1. Have you successfullly uploaded input files?</font></p>
                                                                   <p><font size="4">2. Was your layout file valid?</font></p>
                                                                   <p><font size="4">3. Were your gpr files valid?</font></p>
                                                                   <p><font size="4">4. Did you select correct data type in "Data Upload" part?</font></p>
                                                                   <br> </p>')
                                                             
                                            ),
                                            
                                            conditionalPanel(condition="output.fileUploaded",
                                                             DT::dataTableOutput('tables2')
                                                             
                                            )
                           )
                           
                  ),
                  
                  
                  tabPanel("QC reports",
                           conditionalPanel(condition="(input.generate_table_sum==0&&input.qc_graph==0)",
                                            h2("Quality Control report"),
                                            HTML('<p> <br> You can generate and view your Quality Control (QC) report using our application. Our application allows the users to generate three different graphs about the controls: <strong>Line chart for Ig Controls</strong>, <strong>Line chart for Anti-Ig Controls</strong>, and <strong>Boxplot for PBS Controls</strong>. 
		  		                                        Simultaneously, you could generate and view the <strong>Table Summary of the Antibodies</strong>.
		                                              <br><br><p><strong>Line chart for Ig Controls:</strong> The application uses Measuring statistic option-<strong>Mean</strong>, Background Correction option-<strong>Substraction</strong>, Quality Control option-<strong>KOC</strong>, Batch Correction option-<strong>Alone</strong>, and Normalization option-<strong>Alone</strong> as defaults to generate line chart for Ig Controls reactively. 
		                                              The Line Chart legend status option allows the user to hide Line Chart legend of the line chart.</p>
		                                              
		                                              <br><p><strong>Line chart for Anti-Ig Controls:</strong> The application uses Measuring statistic option-<strong>Mean</strong>, Background Correction option-<strong>Substraction</strong>, Quality Control option-<strong>KOC</strong>, Batch Correction option-<strong>Alone</strong>, and Normalization option-<strong>Alone</strong> as defaults to generate line chart for Anti-Ig Controls reactively. 
		                                              The Line Chart legend status option allows the user to hide Line Chart legend of the line chart.</p>
		                                              
		                                              <br><p><strong>Boxplot for PBS Controls:</strong> The application allows the users to generate and view the boxplot of PBS. Every PBS will have its own boxplot. Our application integrates all boxplots into one big graph. </p>
		                                              
		                                              <br><p><strong>Table Summary of the Antibodies:</strong> This application allows the users to generate and view table summary of the antibodies. We offer two different table summary type: <strong>Table summary corresponding to samples</strong> and <strong>Table summary corresponding to antigens</strong>.  
		                                              The attribute of the table summary includes i) Min, ii) 1st Quantile, iii) Median, iv) Mean, v) 3rd Quantile, and vi) Max.</p>
		                                              
		                                              <br><p><strong>Download Option:</strong> This application allows the users to download Line chart for Ig Controls, Line chart for Anti-Ig Controls, Boxplot for PBS Controls, and Table Summary of the Antibodies generated through our application. After choosing the download type, you could click "Download" button to download our analysis.<br> </p>')
                           ),
                           
                           
                           conditionalPanel(condition="input.qc_graph!=0||input.generate_table_sum!=0",
                                            conditionalPanel(condition="!output.fileUploaded",
                                                             br(),br(),br(),br(),br(),br(),br(),br(),
                                                             HTML('<h2 style="color:red"> The application failed to generate the matrix.</h2>'),
                                                             HTML('<p> <br><strong> <font size="4">Please check the following things:</strong></font></p>
		  		                                                               <p><font size="4">1. Have you successfullly uploaded input files?</font></p>
		                                                                     <p><font size="4">2. Was your layout file valid?</font></p>
		                                                                     <p><font size="4">3. Were your gpr files valid?</font></p>
		                                                                     <p><font size="4">4. Did you select correct data type in "Data Upload" part?</font></p>
		                                                                     <br> </p>')
                                                             
                                            ),
                                            
                                            conditionalPanel(condition="output.fileUploaded",
                                                             #tableOutput('ControlsTable'),
                                                             conditionalPanel(condition="(input.qc_graph!=0)",
                                                                              br(),br(),br(),br(),br(),br(),br(),br(),
                                                                              plotOutput("line_chart_plot")
                                                             ),
                                                             conditionalPanel(condition="(input.generate_table_sum!=0)",
                                                                              br(),br(),br(),br(),
                                                                              tableOutput("Table_Summary")
                                                             )
                                            )
                           )
                  ),
                  
                  tabPanel("Heatmap",
                           tags$a(id = 'downloadData', class = paste("btn btn-default shiny-download-link",'mybutton'), href = "", target = "_blank", download = NA, icon("clone"), 'Download Heatmap as HTML'),
                           tags$head(tags$style(".mybutton{color:white;background-color:blue;} .skin-black .sidebar .mybutton{color: green;}") ),
                           plotlyOutput("heatout",height='1800px', width = '1300px')
                  ),
                  
                  tabPanel("Differential Expression",
                           shiny::htmlOutput("ggAperturaTiendas")
                  ),
                  
                  tabPanel("FAQ",
                           HTML('<h2 style="color:black"> <Strong>Frequently Asked Questions</Strong></h2>'),
                           HTML('<br><br><p><strong> <font size="4">Q: What type of data could be uploaded and graphed in "Heatmap" part?</strong></font></p>
                                 <p><strong><font size="4">A: Any matrix or data frame with numberic values in "csv" type is accepted by the application. </font></strong></p>
                                 <br><br><br><p><strong> <font size="4">Q: Will this application add more functions in the future?</strong></font></p>
                                 <p><strong><font size="4">A: Our application will keep adding more features and functionalities in the future. It will involves more algorithms and more graphs options in the future.</font></strong></p>
                                 <br> </p>')
                  )
      )
    )
  )
)
