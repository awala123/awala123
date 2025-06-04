#
# This is the user-interface definition of a Shiny web application. You can
# run the application by clicking 'Run App' above.
#
# Find out more about building applications with Shiny here:
#
#    https://shiny.posit.co/
#

# Load necessary libraries
library(shiny)
library(shinydashboard)
library(ggplot2)
library(variantspark)
library(tidyverse)
library(genomic.autocorr)
library(DT)
library(vcfR)
library(yaml)
library(ape) # For phylogenetic analysis
library(VariantAnnotation)


# Define UI for application using shinydashboard
ui <- dashboardPage(
  dashboardHeader(title = "Plant Genome Visualization"),
  dashboardSidebar(
    sidebarMenu(
      menuItem("Genome Visualization", tabName = "genome_viz", icon = icon("leaf")),
      menuItem("Phylogenetic Analysis", tabName = "phylo_analysis", icon = icon("tree"))
    )
  ),
  dashboardBody(
    tabItems(
      # Tab for Genome Visualization
      tabItem(tabName = "genome_viz",
              fluidRow(
                selectInput("server1", "Server:",
                            c("Breedbase"="https://musabase.org", "BrAPI test server" = "https://test-server.brapi.org")),
                fileInput('file1', 'Choose VCF File', accept = c('.vcf')),
                actionButton("btn_viz", "Process and Visualize")
              ),
              fluidRow(
                DTOutput('variantTable'),
                plotOutput("genomePlot")
              )
      ),
      # Tab for Phylogenetic Analysis
      tabItem(tabName = "phylo_analysis",
              fluidRow(
                fileInput('file2', 'Choose Newick File', accept = c('.nwk', '.newick')),
                actionButton("btn_phylo", "Load and Analyze")
              ),
              fluidRow(
                plotOutput("phyloTree")
              )
      )
    )
  )
)

# Server logic to process data and generate visualizations
server <- function(input, output) {
  
  observeEvent(input$btn_viz, { print(input$server1)
    req(input$server1)
    
    # Read the VCF file
    vcf <- readVcf(input$file1$datapath, genome = "plant_genome")

    #set brapi
    brapi_url  <- paste0(input$server1, "/brapi/v2")
    call_url  <- paste0(brapi_url, "/variantsets", "/811p14","/calls","?page=0&pageSize=5000")

    #make request
    req <- httr2::request(utils::URLencode(call_url))
    req <- httr2::req_method(req, "GET")
    req <- httr2::req_headers(req, "Accept-Encoding" = "gzip, deflate")
    
    #handle repsonse
    response <- httr2::req_perform(req)
    flatten_results <- jsonlite::fromJSON(httr2::resp_body_string(response), flatten = TRUE)$result$data
    fr_df = data.frame(flatten_results$variantName,flatten_results$genotype.values,flatten_results$callSetName)
    wide_df <- pivot_wider(fr_df, names_from = flatten_results.callSetName, values_from = flatten_results.genotype.values)
    dataMat <- as.matrix(wide_df[,-1])
    rownames(dataMat) <- wide_df$flatten_results.variantName

    # Convert the VCF data to a data frame
    vcf_df <- as.data.frame(info(vcf))
   
    # Output the variant table
    output$variantTable <- renderDT({
      dataMat #vcf_df
    }, options = list(pageLength = 20))
    
    ##### plot
    variant_url  <- paste0(brapi_url, "/variantsets", "/811p14","/variants","?page=0&pageSize=5000")
    
    #make request
    req <- httr2::request(utils::URLencode(variant_url))
    req <- httr2::req_method(req, "GET")
    req <- httr2::req_headers(req, "Accept-Encoding" = "gzip, deflate")
    
    #handle repsonse
    response <- httr2::req_perform(req)
    flatten_results2 <- jsonlite::fromJSON(httr2::resp_body_string(response), flatten = TRUE)$result$data

    
    # Generate the genome plot
    output$genomePlot <- renderPlot({
      req(flatten_results2)
      # Plotting code for genome visualization
      ggplot(dataMat, aes(y = flatten_results2$start, x = flatten_results2$referenceName, color = flatten_results2$referenceBases)) +
        geom_point() +
        theme_minimal() +
        labs(y = "Position", x = "Chromosome", title = "Plant Genome Variants")
    })
  })
  
  observeEvent(input$btn_phylo, {
    req(input$file2)
    
    # Read the Newick file for phylogenetic tree
    phylo_tree <- read.tree(input$file2$datapath)
    
    # Generate the phylogenetic tree plot
    output$phyloTree <- renderPlot({
      plot(phylo_tree, main = "Phylogenetic Tree")
    })
  })
}

# Run the application
shinyApp(ui = ui, server = server, options = list(height = 1000, width=1000))
