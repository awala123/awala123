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
library(plotly)
library("ggbio")
library("circlize")
library("igraph")
library("ggraph")
library(FactoMineR)
library(factoextra)
library(ape)
library(ggtree)


# Define UI for application using shinydashboard
ui <- dashboardPage(
  dashboardHeader(title = "Plant Genome Visualization"),
  dashboardSidebar(
    sidebarMenu(
      menuItem("Genome Visualization", tabName = "genome_viz", icon = icon("leaf")),
      menuItem("SNP Density Plot", tabName = "snp_density", icon = icon("chart-bar")),
      menuItem("Manhattan Plot", tabName = "manhattan_plot", icon = icon("map")),
      menuItem("Genomic Variation Heatmap", tabName = "heatmap_plot", icon = icon("fire")),
      menuItem("Phylogenetic Analysis", tabName = "phylo_analysis", icon = icon("tree")),
      menuItem("Phylogenetic Bootstrap", tabName = "bootstrap_analysis", icon = icon("tree")),
      menuItem("Circular Genome", tabName = "circular_genome", icon = icon("circle-notch")),
      menuItem("Co-expression Network", tabName = "coexpression_network", icon = icon("project-diagram")),
      menuItem("Functional Annotations", tabName = "functional_annotations", icon = icon("dna")),
      menuItem("PCA Analysis", tabName = "pca_analysis", icon = icon("chart-line"))
    )
  ),
  dashboardBody(
    tabItems(
      # Tab for Genome Visualization
      tabItem(tabName = "genome_viz",
              fluidRow(
                selectInput("server1", "Server:",
                            c("Breedbase"="https://musabase.org", "BrAPI test server" = "https://test-server.brapi.org")),
                #fileInput('file1', 'Choose VCF File', accept = c('.vcf')),
                actionButton("btn_viz", "Process and Visualize")
              ),
              fluidRow(
                DTOutput('variantTable'),
                plotOutput("genomePlot")
              )
      ),
      #Tab for SNP Plots
      tabItem(tabName = "snp_density",
              fluidRow(
                #fileInput('file3', 'Choose VCF File', accept = c('.vcf')),
                actionButton("btn_density", "Generate SNP Density Plot")
              ),
              fluidRow(
                plotOutput("densityPlot")
              )
      ),
      # Tab for Manhattan plot
      tabItem(tabName = "manhattan_plot",
              fluidRow(
                #fileInput('file4', 'Choose VCF File', accept = c('.vcf')),
                actionButton("btn_manhattan", "Generate Manhattan Plot")
              ),
              fluidRow(
                plotOutput("manhattanPlot")
              )
      ),
      
      #Tab for Heatmap plot
      tabItem(tabName = "heatmap_plot",
              fluidRow(
                #fileInput('file6', 'Choose VCF File', accept = c('.vcf')),
                numericInput("bin_size","Set bin size:", value=5000, min=1000, max=50000, step=1000),
                textInput("chromosome", "Chr", value = "chr02", width = NULL, placeholder = NULL),
                actionButton("btn_heatmap", "Generate Heatmap")
              ),
              fluidRow(
                plotOutput("heatmapPlot", height = "600px")
              )
      ),
      # Tab for Circular Genome Visualization
      tabItem(tabName = "circular_genome",
              fluidRow(
                fileInput('file5', 'Choose VCF File', accept = c('.vcf')),
                sliderInput("zoom_region", "Zoom into Genomic Region", min = 1, max = 5000000, value = c(1, 100000)),
                actionButton("btn_circular", "Generate Circular Genome Plot")
              ),
              fluidRow(
                plotOutput("circularGenomePlot", height = "600px")
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
      ),
      #Tab for Phylogenetic Bootstrap
      tabItem(tabName = "bootstrap_analysis",
              fluidRow(
                fileInput('file10', 'Upload Phylogenetic Tree (Newick Format)', accept = c('.nwk', '.newick')),
                sliderInput("bootstrap_replicates", "Number of Bootstrap Replicates:", min = 100, max = 1000, value = 500),
                actionButton("btn_bootstrap", "Compute Bootstrap Confidence")
              ),
              fluidRow(
                plotOutput("bootstrapTreePlot", height = "600px")
              )
      ),
      #Tab for coexpression_network
      tabItem(tabName = "coexpression_network",
              fluidRow(
                fileInput('file7', 'Upload Expression Data (CSV)', accept = c('.csv')),
                actionButton("btn_network", "Generate Network")
              ),
              fluidRow(
                plotOutput("networkPlot", height = "600px")
              )
      ),
      
      #Tab for Functional Annotation
      tabItem(tabName = "functional_annotations",
              fluidRow(
                fileInput('file8', 'Upload Gene Annotation File (GFF/GTF)', accept = c('.gff', '.gtf')),
                textInput("gene_search", "Search for Gene:", ""),
                actionButton("btn_annotations", "Load Annotations")
              ),
              fluidRow(
                plotOutput("annotationPlot", height = "600px")
              )
      ),
      
      #Tab for PCA Analysis
      tabItem(tabName = "pca_analysis",
              fluidRow(
                fileInput('file9', 'Upload Genomic Data (CSV)', accept = c('.csv')),
                selectInput("pc_x", "Select X-axis PCA Component:", choices = c("PC1", "PC2", "PC3")),
                selectInput("pc_y", "Select Y-axis PCA Component:", choices = c("PC2", "PC3", "PC4")),
                actionButton("btn_pca", "Run PCA")
              ),
              fluidRow(
                plotlyOutput("pcaPlot", height = "600px")
              )
      )
    )
  )
)

# Server logic to process data and generate visualizations
server <- function(input, output) {
  
  observeEvent(input$btn_viz, {
    req(input$server1)
    
    # Read the VCF file
    #vcf <- readVcf(input$file1$datapath, genome = "plant_genome")

    #set brapi
    brapi_url  <- paste0(input$server1, "/brapi/v2")
    call_url  <- paste(brapi_url, "variantsets", "811p14","calls","?page=0&pageSize=500",sep="/")

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
    #vcf_df <- as.data.frame(info(vcf))
   
    # Output the variant table
    output$variantTable <- renderDT({
      dataMat
    }, options = list(pageLength = 20, scrollX = TRUE))
    
    ##### plot
    variant_url  <- paste(brapi_url, "variantsets", "811p14","variants","?page=0&pageSize=50",sep="/")
    
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
      ggplot(flatten_results2, aes(y = flatten_results2$start, x = flatten_results2$referenceName, color = flatten_results2$referenceBases)) +
        geom_point() +
        theme_minimal() +
        labs(y = "Position", x = "Chromosome", title = "Plant Genome Variants",color = "Ref") 
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
  
  observeEvent(input$btn_bootstrap, {
    req(input$file10)
    
    # Read the phylogenetic tree
    phylo_tree <- read.tree(input$file10$datapath)
    
    # Compute bootstrap support
    bootstrapped_tree <- boot.phylo(phylo_tree, phylo_tree$tip.label, B = input$bootstrap_replicates)
    
    # Add bootstrap values to tree
    phylo_tree$node.label <- bootstrapped_tree
    
    # Generate bootstrap phylogenetic tree visualization
    output$bootstrapTreePlot <- renderPlot({
      ggtree(phylo_tree, aes(color = as.numeric(node.label))) +
        geom_tiplab() +
        theme_minimal() +
        scale_color_gradient(low = "red", high = "blue") +
        labs(title = "Phylogenetic Tree with Bootstrap Confidence Intervals")
    })
  })
  
  brapi_url  <- paste("https://musabase.org", "brapi/v2", sep="/")
  variant_url  <- paste(brapi_url, "variantsets", "811p14","variants","?page=0&pageSize=500", sep="/")
  
  #make request
  req <- httr2::request(utils::URLencode(variant_url))
  req <- httr2::req_method(req, "GET")
  req <- httr2::req_headers(req, "Accept-Encoding" = "gzip, deflate")
  
  #handle repsonse
  response <- httr2::req_perform(req)
  variant_vcf <- jsonlite::fromJSON(httr2::resp_body_string(response), flatten = TRUE)$result$data
  
  observeEvent(input$btn_density, {
    #req(input$file3)
    
    # Read the VCF file
    #vcf <- readVcf(input$file3$datapath, genome = "plant_genome")
    
    # Extract position data and summarize SNP density
    #snp_positions <- as.numeric(info(vcf)$POS)
    snp_positions <- variant_vcf$start
    
    # Create a density plot
    output$densityPlot <- renderPlot({
      ggplot(data.frame(Position = snp_positions), aes(x = Position)) +
        geom_density(fill = "blue", alpha = 0.4) +
        theme_minimal() +
        labs(x = "Genomic Position", y = "Density", title = "SNP Density Plot")
    })
  })
  
  observeEvent(input$btn_manhattan, {
    #req(input$file4)

    # Generate Manhattan plot
    output$manhattanPlot <- renderPlot({
      ggplot(variant_vcf, aes(x = variant_vcf$start, y = -log10( runif(33, 0, 1) ), color = variant_vcf$referenceName)) +
        geom_point() +
        theme_minimal() +
        labs(x = "Genomic Position", y = "-log10(P-value)", title = "Genome-wide Manhattan Plot")
    })
  })
 
  observeEvent(input$btn_circular, {
    req(input$file5)
    
    # Read the VCF file
    vcf <- readVcf(input$file5$datapath, genome = "plant_genome")
    print(vcf)
    # Convert VCF data into GRanges for visualization
    gr <- as(vcf, "GRanges")
    
    # Apply zoom filter
    gr_filtered <- gr[seqnames(gr) %in% paste0("chr", 1:12) & 
                        start(gr) >= input$zoom_region[1] & start(gr) <= input$zoom_region[2]]
    
    # Generate circular genome visualization
    output$circularGenomePlot <- renderPlot({
      autoplot(gr_filtered, layout = "circular", aes(fill = allele_frequency(gr))) +
        theme_minimal() +
        labs(title = "Circular Genome Visualization with SNPs and Annotations")
    })
  })
 
  observeEvent(input$btn_heatmap, {
    #req(input$file6)
    
    # Read the VCF file
    #vcf <- readVcf(input$file6$datapath, genome = "plant_genome")
    #print(variant_vcf)
    vcf_df = data.frame(variant_vcf$start,variant_vcf$referenceName,runif(33,0,1))
    colnames(vcf_df) = c("POS","CHROM","VAR")
    print(vcf_df)
    
    # # Extract relevant genomic data
    # vcf_df <- data.frame(POS = as.numeric(info(vcf)$POS),
    #                      CHROM = as.factor(info(vcf)$CHROM),
    #                      VAR = runif(nrow(info(vcf)), 0, 1))  # Simulated variation scores
    
    # Filter based on selected chromosome
    #vcf_filtered <- vcf_df[vcf_df$CHROM == input$chromosome, ]
    vcf_filtered = vcf_df
    print(vcf_filtered)
    # Aggregate into bins for heatmap visualization
    #print(input$bin_size)
    vcf_binned <- vcf_filtered %>%
      mutate(Binned_POS = floor(vcf_filtered$POS / input$bin_size) * input$bin_size) %>%
      group_by(Binned_POS) %>%
      summarize(Mean_VAR = mean(vcf_filtered$VAR), nrow = 2)
    
    # Generate Heatmap
    output$heatmapPlot <- renderPlot({
      heatmap_matrix <- matrix(vcf_filtered$VAR,nrow=6)
      print(vcf_binned)
      print(heatmap_matrix)
      heatmap(heatmap_matrix,
              xlab = "Genomic Variation",
              # col = colorRamp2(c(min(vcf_filtered$VAR), max(vcf_filtered$VAR)), c("blue", "red")),
              main = paste("Genomic Variation Heatmap -", input$chromosome),
              cluster_columns = FALSE)
    })
  }) 
  
  observeEvent(input$btn_network, {
    req(input$file7)
    
    # Read the expression data
    expression_data <- read.csv(input$file7$datapath, row.names = 1)
    
    # Compute correlation matrix
    correlation_matrix <- cor(expression_data, method = "pearson")
    
    # Extract edges above the user-defined threshold
    edge_list <- which(correlation_matrix > input$cor_threshold, arr.ind = TRUE)
    edges <- data.frame(from = rownames(correlation_matrix)[edge_list[,1]],
                        to = colnames(correlation_matrix)[edge_list[,2]],
                        weight = correlation_matrix[edge_list])
    
    # Create igraph object
    graph <- graph_from_data_frame(edges, directed = FALSE)
    
    # Color nodes based on expression level
    node_colors <- rowMeans(expression_data)  
    V(graph)$color <- colorRampPalette(c("blue", "red"))(length(unique(node_colors)))[rank(node_colors)]
    
    # Generate interactive co-expression network plot
    output$networkPlot <- renderPlotly({
      g <- ggraph(graph, layout = input$layout_type) +
        geom_edge_link(aes(edge_alpha = weight), color = "gray") +
        geom_node_point(aes(color = color, size = degree(graph)), show.legend = FALSE) +
        geom_node_text(aes(label = name), repel = TRUE, size = 4) +
        theme_void() +
        labs(title = "Co-expression Network Visualization")
      
      ggplotly(g)
    })
  })
  
  observeEvent(input$btn_annotations, {
    req(input$file8)
    
    # Read the Gene Annotation File
    annotation_data <- import(input$file8$datapath)
    
    # Convert into a GRanges object
    gr_annotations <- as(annotation_data, "GRanges")
    
    # Apply gene search filter
    if (input$gene_search != "") {
      gr_annotations <- gr_annotations[grep(input$gene_search, elementMetadata(gr_annotations)$gene_name, ignore.case = TRUE)]
    }
    
    # Apply zoom filter
    gr_filtered <- gr_annotations[start(gr_annotations) >= input$zoom_region[1] & start(gr_annotations) <= input$zoom_region[2]]
    
    # Generate interactive functional annotation plot
    output$annotationPlot <- renderPlotly({
      g <- autoplot(gr_filtered, aes(fill = feature), layout = "linear") +
        theme_minimal() +
        labs(title = "Functional Annotations Overlay")
      
      ggplotly(g, tooltip = "feature")
    })
  })
  
  observeEvent(input$btn_pca, {
    req(input$file9)
    
    # Read genomic dataset
    genomic_data <- read.csv(input$file9$datapath, row.names = 1)
    
    # Run PCA
    pca_results <- PCA(genomic_data, graph = FALSE)
    
    # Extract PCA coordinates
    pca_df <- data.frame(Sample = rownames(pca_results$ind$coord),
                         PC1 = pca_results$ind$coord[,1],
                         PC2 = pca_results$ind$coord[,2],
                         PC3 = pca_results$ind$coord[,3],
                         PC4 = pca_results$ind$coord[,4])
    
    # Generate interactive PCA plot
    output$pcaPlot <- renderPlotly({
      g <- ggplot(pca_df, aes_string(x = input$pc_x, y = input$pc_y, color = "Sample")) +
        geom_point(size = 4) +
        theme_minimal() +
        labs(title = "PCA Analysis of Genomic Data", x = input$pc_x, y = input$pc_y)
      
      ggplotly(g)
    })
  })
}

# Run the application
shinyApp(ui = ui, server = server, options = list(height = 900, width=1500))
