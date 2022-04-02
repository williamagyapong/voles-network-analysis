#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Purpose: Main application interface                                                      #
# author: William Ofosu Agyapong                                                #
# Last updated: March 22, 2022                                                    #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#        Loading Required Packages
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
library(dplyr)
library(tidyr)
library(stringr)
# library(broom)
library(ggplot2)
library(ggpubr) # interface for adding test results to plot
library(patchwork)
library(readxl) # load package for reading excel files directly into R
library(igraph)
# library(CINNA)
library(kableExtra)
library(knitr)

library(RColorBrewer)
library(shiny)
library(shinytitle)
library(plotly)


#-----------------------------------------------------------------
# Set the current working directory to the file path
# setwd(dirname(rstudioapi::getSourceEditorContext()$path)) 

# Import file dependencies
source("functions.R")


#------------------------------------------------------------------
# Set a default ggplot theme for all plots
theme_set(theme_classic())

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#        Loading and Preprocessing Data
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  load("net_data.RData")
  
  
  # prepare data for adjacency matrix
  IL_More_Social <- IL_More_Social_raw[, -1]
  KI_Less_Social <- KI_Less_Social_raw[, -1]
  
  #-- Create graph object from adjacency matrix
  IL_net <- graph.adjacency(as.matrix(IL_More_Social), mode="undirected", weighted = T)
  V(IL_net)$name <- gsub("'", "", V(IL_net)$name) # clean names
  
  KI_net <- graph.adjacency(as.matrix(KI_Less_Social), mode="undirected", weighted = T)
  V(KI_net)$name <- gsub("'", "", V(KI_net)$name)
  
  # get subregions for user input
  subregions <- sort(unique(Defined_Regions$Subregion))
  
  # subregion_df <<- Defined_Regions %>%
  #   filter(Subregion == stringr::str_to_title("midbrain")) %>%
  #   mutate(ROI = gsub("'", "", ROI))
  
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#                  Visualization of content with Shiny App
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  # Define server logic ----
  server <- function(input, output, session) {
    
    # theme for plots
    thematic::thematic_shiny()
    
    # change browser tab title
    change_window_title(title = "Voles Network Analysis")
    
    # get the subregion from the user
    # subregion <- reactive(input$subregion)
    # select region of interest 
    subregion_info <<- reactive({
      Defined_Regions %>%
        filter(Subregion == stringr::str_to_title(input$subregion)) %>%
        mutate(ROI = gsub("'", "", ROI))
    })
    
    # Get colors
    IL_col <- reactive({
      c(input$subregion_col, input$IL_col)
    })
    
    KI_col <- reactive({
      c(input$subregion_col, input$KI_col)
    })
    
    # displaying full networks
    #______________________________________________________________________________________________
    
    output$full_networks <- renderPlot({
      # subregion_df <- subregion_info()
      par(mfrow= c(1,2))
      
      pgraph(IL_net, color = IL_col(), toggle.label = input$toggle_label)
      pgraph(KI_net, color = KI_col(), toggle.label = input$toggle_label)
    })
    
    output$full_networks2 <- renderPlot({
      par(mfrow= c(1,2))
      pgraph(get_subgraph(IL_net, dvwne = F), color = IL_col(),toggle.label = input$toggle_label)
      pgraph(get_subgraph(KI_net, dvwne = F), color = KI_col(), toggle.label = input$toggle_label)
    })
    
    # displaying subgraphs
    #______________________________________________________________________________________________
    output$sub_networks<- renderPlot({
      par(mfrow = c(1,2))
      
      psubgraph(get_subgraph(IL_net, dvwne = T), centrality = str_to_lower(input$centrality_type3),
                color = IL_col(), vlabel.size = input$vlabel_size, legend = T)
      psubgraph(get_subgraph(KI_net, dvwne = T), centrality = str_to_lower(input$centrality_type3),
                color = KI_col(), vlabel.size = input$vlabel_size, legend = T)
    })
    
    # Distributions
    #______________________________________________________________________________________________
    output$distributions <- renderPlotly({
      
      compare_distrib(IL_net, KI_net, str_to_lower(input$centrality_type),
                      plt = str_to_lower(input$dist_type), col_vals=c(IL_col()[2], KI_col()[2])) 
      
    })
    
    #__________________________________________________________________________
    output$diff_test <- renderPlotly({
      
      compare_distrib(IL_net, KI_net, str_to_lower(input$centrality_type2), plt = "boxplot",
                      annotate.test = T, method = input$test_method,  col_vals = c(IL_col()[2], KI_col()[2]))
    })
    
    
    output$lookup_table <- function(){
      
      subregion_info() %>%
        relocate(Subregion, .before = "ROI") %>%
        mutate(Notation = abbr(ROI)) %>%
        mutate(d1 = degree(IL_net)[ROI],
               d2 = degree(KI_net)[ROI],
               c1 = closeness(IL_net)[ROI],
               c2 = closeness(KI_net)[ROI],
               b1 = betweenness(IL_net)[ROI],
               b2 = betweenness(KI_net)[ROI],
               e1 = eigen_centrality(IL_net)$vector[ROI],
               e2 = eigen_centrality(KI_net)$vector[ROI]
        ) %>%
        kable(format = "html", col.names = c("Subregion", "ROI", "Abbr",rep(c("IL", "KI"), times = 4)), 
              align = "lllcccccccc",
              ) %>%
        add_header_above(c("", "", "", "Degree"=2, "Closeness"=2, "Betweenness"=2, "Eigenvector"=2)) %>%
        kable_styling(bootstrap_options = c("striped", "hover"))
    }
    
    output$full_net_xtics <- function() {
      measures <- c("Number of Vertices", "Number of Edges", "Density (edges)", "Average path length", "Diameter")
      IL_measures <- c(vcount(IL_net), gsize(IL_net), edge_density(IL_net), mean_distance(IL_net, directed = F), diameter(IL_net, directed = F))
      
      KI_measures <- c(vcount(KI_net), gsize(KI_net), edge_density(KI_net), mean_distance(KI_net, directed = F), diameter(KI_net, directed = F))
      
      # display results
      data.frame(measures, IL_measures, KI_measures) %>%
        kable(format = "html", col.names = c("Characteristic", "IL Network", "KI Network"), align = "lcc",
              caption = "Characteristics of the two full networks") %>%
        kable_paper() %>%
        kable_styling(bootstrap_options = c("striped", "hover"))
    }
    
    output$sub_net_xtics <- function() {
      IL_sub_net <- get_subgraph(IL_net, dvwne = T)
      KI_sub_net <- get_subgraph(KI_net, dvwne = T)
      measures <- c("Number of Vertices", "Number of Edges", "Density (edges)", "Average path length", "Diameter")
      IL_measures <- c(vcount(IL_sub_net), gsize(IL_sub_net), edge_density(IL_sub_net), 
                       mean_distance(IL_sub_net, directed = F), diameter(IL_sub_net, directed = F))
      
      KI_measures <- c(vcount(KI_sub_net), gsize(KI_sub_net), edge_density(KI_sub_net), 
                       mean_distance(KI_sub_net, directed = F), diameter(KI_sub_net, directed = F))
      
      # display results
      data.frame(measures, IL_measures, KI_measures) %>%
        kable(format = "html", col.names = c("Characteristic", "IL Subnetwork", "KI Subnetwork"), align = "lcc",
              caption = "Characteristics of the two subregion networks") %>%
        kable_paper() %>%
        kable_styling(bootstrap_options = c("striped", "hover"))
    }
  }
  

# Define User Interface (UI) ----
ui <- fluidPage(theme=bslib::bs_theme(bootswatch="united"),
                tags$head(
                  tags$style(
                    HTML(".table-hover > tbody > tr:hover { 
                      background-color: #ffff99;
                    }")
                  )
                ),

                # used here to make change_window_title() in the server side work
                use_shiny_title(),
    titlePanel(
      # creating a title bar
      div("Voles Brain Network Analysis",
          style="text-align:center;width:100%;height:60px;margin-top:-15px;font-weight:bold;
        font-size: 20px;font-family: Tahoma; background-color: #E1EDEB; line-height:60px;
        color:gray;")

    ),

    tabsetPanel(
      tabPanel("Home",
               # shared controls area
               fluidRow(style = "margin-top:20px",
                 column(style="border-bottom:2px #E1EDEB solid;padding-bottom:1px;",
                        width=3,
                        # controller for trend
                        selectInput("subregion", "Subregion:", subregions, 
                                    selected = "Midbrain"),
                 ),
                 column(style="border-bottom:2px #E1EDEB solid;padding-bottom:1px;",
                        width=3,
                        
                        selectInput("IL_col", " IL Color:",
                                    c("gold", "orange", "orange1", "orange2","yellow")),
                       
                 ),
                 column(style="border-bottom:2px #E1EDEB solid;padding-bottom:1px;",
                        width=3,
                        selectInput("KI_col", "KI Color:",
                                    c("dodgerblue", "blue", "navy", "cyan")),
                 ),
                 column(style="border-bottom:2px #E1EDEB solid;padding-bottom:1px;",
                        width=3,
                        selectInput("subregion_col", "Subregion Color:",
                                    c("tomato", "tomato1", "red", "red1", "red2")),
                 )
               ),
               #first row
               fluidRow(
                 column(style="margin-top:0px;padding-bottom:15px;border-bottom:2px #E1EDEB solid;",
                        width=6,
                        h3("Full Networks",
                           style=" color: #B1967C;font-family: 'Arial Black', serif;font-weight: normal;font-size: 1.2em;
             border-bottom: 2px #E1EDEB solid;padding-left: 25px;"),

                        radioButtons("toggle_label", label = "", inline = T,
                                     choices = c("No label"="none", "Short Label"="short", "Full label"="full")),
                        plotOutput("full_networks"),
                        plotOutput("full_networks2")
                        
                 ),

                 column(style="margin-top:0px;border-left:2px #E1EDEB solid;border-bottom:2px #E1EDEB solid;padding-bottom:20px;",
                        width=6,
                        h3("Sub Networks",
                           style=" color: #B1967C;font-family: 'Arial Black', serif;font-weight: normal;font-size: 1.2em;
              border-bottom: 2px #E1EDEB solid;padding-left: 25px;"),
                        fluidRow(
                          column(width = 6,
                                 selectInput("centrality_type3", "Centrality",
                                             c("Degree", "Closeness", "Betweenness", "Eigenvector"))
                          ),
                          column(width = 6,
                                 sliderInput("vlabel_size", "Size of vertex label", value = 1, min = 0.6, max = 4)
                          )
                        ),
                        plotOutput("sub_networks")
                 )
               )
               ,
               fluidRow(

               ),
               #second row
               fluidRow(
                 column(style="padding-top:20px;",
                        width=6,
                        h3("Distribution of Centrality",
                           style=" color: #B1967C;font-family: 'Arial Black', serif;font-weight: normal;font-size: 1.2em;
              border-bottom: 2px #E1EDEB solid;padding-left: 25px;"),

                        fluidRow(
                          column(width = 6,
                                 radioButtons("dist_type", "Distribution", c("Density", "Boxplot", "Violin"))
                          ),
                          column(width = 6,
                                 selectInput("centrality_type", "Centrality",
                                             c("Degree", "Closeness", "Betweenness", "Eigenvector"))
                          )
                        ),
                        plotlyOutput("distributions"),
                 ),

                 column(style="padding-top:20px; margin-bottom:20px;border-left: 2px #E1EDEB solid;",
                        width=6,
                        h3("Is the difference significant?",
                           style=" color: #B1967C;font-family: 'Arial Black', serif;font-weight: normal;font-size: 1.2em;border-bottom: 2px #E1EDEB solid;padding-left: 25px;"),
                        
                        fluidRow(
                          column(6,
                                 selectInput("centrality_type2", "Centrality",
                                             c("Degree", "Closeness", "Betweenness", "Eigenvector"))),
                          column(6,
                                 selectInput("test_method", "Testing Method",
                                             c("Wilcoxon signed rank test"="wilcox", "Paired t-test"="t.test")))
                        ),
                        plotlyOutput("diff_test")
                 )
               )

      ),

      tabPanel("Vertex Label Lookup",
               fluidRow(style = "margin:auto; width:80%; padding-bottom:100px;padding-top:50px;",
                 tableOutput("lookup_table")
               )
      ),

      tabPanel("Network Characteristics",
               fluidRow(style = "padding-bottom:100px;padding-top:50px;",
                 column(width = 6, style="border-right: 4px #E1EDEB solid;",
                        
                        tableOutput("full_net_xtics")
                 ),
                 column(width = 6,
                        
                        tableOutput("sub_net_xtics")
                 )
               )
      ),
      tabPanel("About the App",
               fluidRow(
                 column(width = 12,style="height:350px;padding:100px;margin-bottom:360px;",
                        HTML("Utilizing functional magnetic resonance imaging (fMRI) data from culturally 
                        and behaviorally distinct populations of prairie voles generated by 
                        <a href='https://www.sciencedirect.com/science/article/pii/S2451902221003207' target='blank'>Ortiz et al. (2021)</a>,
                        this web application provides a platform for conducting a statistical comparison of each of the centrality measures
                        (degree, closeness, betweenness, eigenvector centrality), comparing the two 
                        populations of voles (IL and KI). IL and KI denote the prairie male voles
                        from Illinois and male cross-breed off-springs of Kansas dam and Illinois sires,
                        respectively. The goal is to determine if the two populations of prairie voles 
                        exhibit any significant differences in prosocial behavior evident from differences in 
                        the centrality measures for the midbrain subregion networks. <br><hr>"),
                        plotlyOutput("about-app")
                 )
               )
      )
    ),
  # footer
  fluidRow(class = "footer",
            style="box-sizing: border-box;width:100%; height: 130px;
            background-color: whitesmoke; padding:8px;",
    column(width = 5,
              HTML('<span style="font-weight:bold;font-style:italic;">Data Source:</span> <br>
             Provided by Richard J. Ortiz, Research Technician @ Biological Science, UTEP.')
           ),
    column(width = 4, style="border-right:2px solid #fff; border-left:2px solid #fff;pading-left:10px;",
            HTML('<span style="font-weight:bold;font-style:italic;">Developed By:</span><br>
            William O. Agyapong <br>
           PhD (Data Science) student <br>University of Texas at El Paso <br>
           <span> &copy; 2022 </span>')
           ),
    column(width = 3, style="padding-left:10px;",
            HTML('<span style="font-weight:bold;font-style:italic;">Follow on:</span> <br>
              <ul>
                <li><a href="https://gitHub.com/williamagyapong">Github</a></li>
                <li><a href="https://www.datacamp.com/profile/williamagyapong">DataCamp</a></li>
                <li><a href="https://linkedin.com/in/william-agyapong-372aa4146">LinkedIn</a></li>
              </ul>')
    )
  )

)


# Run the app ----
shinyApp(ui = ui, server = server)

############################# END OF FILE ######################
