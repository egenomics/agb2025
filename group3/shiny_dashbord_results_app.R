
###### AGB RESULTS DELIVERY #########


# Load required libraries
library(shiny)
library(shinydashboard)
library(ggplot2)
library(phyloseq)
library(tidyverse)
library(vegan)
library(readr)
library(DT)
library(ape)
library(bslib)
library(plotly)
library(reshape2)
library(Rtsne)
library(umap)

# Define UI
ui <- dashboardPage(
  skin = "blue", 
  
  # Header from dashboard style but with navbar styling 
  dashboardHeader(
    title = div(icon("microscope"), "AGB Microbiome Results"), 
    titleWidth = 300
  ), 
  
  # Sidebar with dashboard style 
  dashboardSidebar(
    width = 300, 
    sidebarMenu(
      menuItem(" Home", tabName = "home", icon = icon("home")), 
      menuItem(" Data Upload", tabName = "data", icon = icon("upload")),
      menuItem(" Taxonomy", tabName = "taxonomy", icon = icon("bacteria")),
      menuItem(" Alpha Diversity", tabName = "alpha", icon = icon("chart-line")),
      menuItem(" Beta Diversity", tabName = "beta", icon = icon("project-diagram")),
      # menuItem("️ Functional Pathways", tabName = "function", icon = icon("dna")),
      menuItem(" Metadata Explorer", tabName = "metadata", icon = icon("table")),
      menuItem(" Reports", tabName = "reports", icon = icon("file-export"))
    )
  ),
  
  # Body with dashboard style 
  dashboardBody(
    # Use BS theme for consistent styling 
    tags$head(
      tags$link(rel = "stylesheet", type = "text/css", href = "https://cdnjs.cloudflare.com/ajax/libs/font-awesome/5.15.4/css/all.min.css"), 
      tags$style(HTML("
                      .box.box-solid.box-primary>.box-header {
                        background: #3c8dbc;
                        color: #fff;
                      }
                      .nav-tabs-custom>.nav-tabs>li.active {
                        border-top-color: #3c8dbc;
                      }
                      .content-wrapper, .right-side {
                        background-color: #f4f6f9;
                      }
                      .box {
                        border-radius: 3px;
                        box-shadow: 0 2px 5px rgba(0,0,0,0.1);
                      }
                      .well {
                        background-color: #fff;
                        border: 1px solid #e3e3e3;
                        border-radius: 3px;
                        box-shadow: 0 1px 3px rgba(0,0,0,0.05);
                      }
                      "))
    ), 
    
    tabItems(
      # Home tab 
      tabItem(tabName = "home", 
              fluidRow(
                box(
                  title = "Welcome to AGB Results", 
                  status = "primary", 
                  solidHeader = TRUE, 
                  width = 12, 
                  div(
                    h3("Comprehensive Microbiome Analysis Dashboard"), 
                    hr(), 
                    fluidRow(
                      column(
                        width = 6, 
                        h4("🔍 Features"), 
                        tags$ul(
                          tags$li(strong("Taxonomic Profiles:"), "Explore the microbial composition at various taxonomic levels"),
                          tags$li(strong("Diversity Analysis:"), "Assess alpha and beta diversity metrics"), 
                          tags$li(strong("Metadata Exploration:"), "Correlate microbial features with sample metadata"), 
                          tags$li(strong("Report Generation:"), "Generate customized reports for your findings")
                        )
                      ), 
                      column(
                        width = 6, 
                        h4("🚀 Getting Started"), 
                        p("Begin by uploading your data files in the Data Upload tab"),
                        tags$ul(
                          tags$li("Upload ASV/OTU table, taxonomy assignments, and sample metadata"),
                          tags$li("Explore taxonomic composition in the Taxonomy tab"), 
                          tags$li("Analyze diversity metrics in the Alpha and Beta Diversity tabs"), 
                          tags$li("Generate reports with key findings")
                        )
                      )
                    )
                  )
                )
              )
    ), 
    
    # Data Upload Tab 
    tabItem(tabName = "data",
            fluidRow(
              box(
                title = "Upload Microbiome Data", 
                status = "primary", 
                solidHeader = TRUE, 
                width = 12, 
                fluidRow(
                  column(
                    width = 6,
                    fileInput("asv_table", "Upload ASV/OTU Table (CSV or TSV format):",
                              accept = c(".csv", ".tsv", ".txt")), 
                    fileInput("taxonomy_file", "Upload Taxonomy File (CSV or TSV format):",
                              accept = c(".csv", ".tsv", ".txt")),
                    fileInput("metadata_file", "Upload Sample Metadata (CSV or TSV format):",
                              accept = c(".csv", ".tsv", ".txt")),
                    actionButton("load_data", "Load Data", 
                                 icon = icon("upload"),
                                 class = "btn-success"),
                    actionButton("generate_example", "Generate Example Data", 
                                 icon = icon("database"),
                                 class = "btn-info")
                  ),
                  column(
                    width = 6,
                    checkboxInput("load_example", "Load Example Data Instead", value = FALSE),
                    conditionalPanel(
                      condition = "input.load_example == true",
                      selectInput("example_dataset", "Select Example Dataset:",
                                  choices = c("Human Gut Microbiome", "Soil Microbiome"))
                    ),
                    hr(),
                    h4("Data Preview:"),
                    dataTableOutput("metadata_preview")
                  )
                )
              )
            ), 
            fluidRow(
              box(
                title = "Data Summary", 
                status = "info", 
                solidHeader = TRUE,
                width = 12,
                verbatimTextOutput("data_summary")
              )
            )
            ),
    
    # Taxonomy Tab 
    tabItem(tabName = "taxonomy",
            fluidRow(
              box(
                title = "Taxonomy Settings", 
                status = "primary", 
                solidHeader = TRUE,
                width = 3,
                selectInput("tax_level", "Taxonomic Level:",
                            choices = c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"),
                            selected = "Phylum"),
                selectInput("group_var", "Group By:",
                            choices = c("None"), selected = "None"),
                sliderInput("min_abundance", "Minimum Abundance (%):", 
                            min = 0, max = 5, value = 1, step = 0.1),
                checkboxInput("sort_by_abundance", "Sort by Abundance", value = TRUE),
                checkboxInput("show_legend", "Show Legend", value = TRUE),
                radioButtons("plot_type", "Plot Type:",
                             choices = c("Stacked Bar" = "bar", "Heatmap" = "heatmap"),
                             selected = "bar"),
                hr(),
                downloadButton("download_tax_plot", "Download Plot", class = "btn-primary")
              ),
              box(
                title = "Taxonomic Composition", 
                status = "info", 
                solidHeader = TRUE,
                width = 9,
                plotOutput("taxonomy_plot", height = 600)
              )
            ),
            fluidRow(
              box(
                title = "Top Taxa Table", 
                status = "info", 
                solidHeader = TRUE,
                width = 12,
                dataTableOutput("taxa_table")
              )
            )
    ),
    
    # Alpha Diversity Tab
    tabItem(tabName = "alpha",
            fluidRow(
              box(
                title = "Alpha Diversity Settings", 
                status = "primary", 
                solidHeader = TRUE,
                width = 3,
                selectInput("alpha_metric", "Alpha Diversity Metric:",
                            choices = c("Shannon", "Observed", "Simpson", "InvSimpson", "Faith PD", "Chao1"),
                            selected = "Shannon"),
                selectInput("alpha_group", "Group Samples By:",
                            choices = c("None"), selected = "None"),
                conditionalPanel(
                  condition = "input.alpha_group != 'None'",
                  selectInput("stat_test", "Statistical Test:",
                              choices = c("None", "ANOVA", "Kruskal-Wallis", "T-test", "Wilcoxon"),
                              selected = "None")
                ),
                sliderInput("rarefaction_depth", "Rarefaction Depth:", 
                            min = 1000, max = 10000, value = 5000, step = 1000),
                checkboxInput("add_boxplot", "Add Box Plot", value = TRUE),
                checkboxInput("add_points", "Show Individual Points", value = TRUE),
                hr(),
                downloadButton("download_alpha_plot", "Download Plot", class = "btn-primary")
              ),
              box(
                title = "Alpha Diversity Plot", 
                status = "info", 
                solidHeader = TRUE,
                width = 9,
                plotlyOutput("alpha_plot", height = 500),
                conditionalPanel(
                  condition = "input.stat_test != 'None'",
                  h4("Statistical Test Results:"),
                  verbatimTextOutput("alpha_stats")
                )
              )
            ),
            fluidRow(
              box(
                title = "Rarefaction Curve", 
                status = "info", 
                solidHeader = TRUE,
                width = 12,
                plotlyOutput("rarefaction_plot", height = 400)
              )
            )
    ),
    
    # Beta Diversity Tab
    tabItem(tabName = "beta",
            fluidRow(
              box(
                title = "Beta Diversity Settings", 
                status = "primary", 
                solidHeader = TRUE,
                width = 3,
                selectInput("beta_metric", "Beta Diversity Metric:",
                            choices = c("Bray-Curtis", "Jaccard", "UniFrac", "Weighted UniFrac"),
                            selected = "Bray-Curtis"),
                selectInput("ordination", "Ordination Method:",
                            choices = c("PCoA", "NMDS", "t-SNE", "UMAP"),
                            selected = "PCoA"),
                selectInput("beta_group", "Color By:",
                            choices = c("None"), selected = "None"),
                selectInput("shape_by", "Shape By:",
                            choices = c("None"), selected = "None"),
                conditionalPanel(
                  condition = "input.beta_group != 'None'",
                  selectInput("permanova", "Run PERMANOVA:",
                              choices = c("No", "Yes"),
                              selected = "No")
                ),
                conditionalPanel(
                  condition ="input.beta_group != 'None'",
                  radioButtons("plot_type", "Visualization Type:",
                               choices = c("Ellipses" = "ellipse", "Boxplots" = "boxplot"), 
                               selected = "ellipse")
                ),
                conditionalPanel(
                  condition = "input.plot_type == 'ellipse' && input.beta_group != 'None'", 
                  checkboxInput("add_ellipse", "Add Confidence Ellipses", value = TRUE)
                ),
                checkboxInput("add_labels", "Add Sample Labels", value = FALSE),
                hr(),
                downloadButton("download_beta_plot", "Download Plot", class = "btn-primary")
              ),
              box(
                title = "Beta Diversity Plot", 
                status = "info", 
                solidHeader = TRUE,
                width = 9,
                plotlyOutput("beta_plot", height = 500),
                conditionalPanel(
                  condition = "input.permanova == 'Yes'",
                  h4("PERMANOVA Results:"),
                  verbatimTextOutput("permanova_results")
                )
              )
            )
    ),
    
    # Metadata Explorer Tab
    tabItem(tabName = "metadata",
            fluidRow(
              box(
                title = "Metadata Exploration", 
                status = "primary", 
                solidHeader = TRUE,
                width = 3,
                selectInput("meta_x", "X-axis Variable:",
                            choices = c("None"), selected = "None"),
                selectInput("meta_y", "Y-axis Variable:",
                            choices = c("None"), selected = "None"),
                selectInput("meta_color", "Color By:",
                            choices = c("None"), selected = "None"),
                radioButtons("meta_plot_type", "Plot Type:",
                             choices = c("Scatter" = "scatter", 
                                         "Box Plot" = "box", 
                                         "Bar Plot" = "bar"),
                             selected = "scatter"),
                hr(),
                downloadButton("download_meta_plot", "Download Plot", class = "btn-primary")
              ),
              box(
                title = "Metadata Visualization", 
                status = "info", 
                solidHeader = TRUE,
                width = 9,
                plotOutput("metadata_plot", height = 500)
              )
            ),
            fluidRow(
              box(
                title = "Correlation Analysis", 
                status = "info", 
                solidHeader = TRUE,
                width = 12,
                plotOutput("correlation_plot", height = 400)
              )
            )
    ),
    
    # Reports Tab
    tabItem(tabName = "reports",
            fluidRow(
              box(
                title = "Generate Reports", 
                status = "primary", 
                solidHeader = TRUE,
                width = 12,
                fluidRow(
                  column(
                    width = 6,
                    h4("Report Options"),
                    checkboxGroupInput("report_sections", "Include Sections:",
                                       choices = c("Taxonomy Analysis" = "taxonomy",
                                                   "Alpha Diversity" = "alpha",
                                                   "Beta Diversity" = "beta",
                                                   "Functional Analysis" = "function",
                                                   "Metadata Exploration" = "metadata"),
                                       selected = c("taxonomy", "alpha", "beta")),
                    radioButtons("report_format", "Report Format:",
                                 choices = c("PDF" = "pdf", 
                                             "HTML" = "html",
                                             "Word Document" = "docx"),
                                 selected = "html"),
                    textInput("report_title", "Report Title:", 
                              value = "Microbiome Analysis Report"),
                    textAreaInput("report_comments", "Additional Comments:", 
                                  rows = 3),
                    hr(),
                    downloadButton("generate_report", "Generate Report", class = "btn-success")
                  ),
                  column(
                    width = 6,
                    h4("Report Preview"),
                    htmlOutput("report_preview")
                  )
                )
              )
            )
    )
    )
    )
    
  )

  
# Server logic with mock data

server <- function(input, output, session) {
  
  # Reactive values to store data 
  phyloseq_obj <- reactiveVal(NULL)
  metadata <- reactiveVal(NULL)
  
  # Function to generate example data similar to Moving Pictures tutorial
  generate_example_data <- function() {
    # Create example OTU table 
    n_samples <- 30
    n_otus <- 200
    
    # Sample names 
    sample_names <- paste0("Sample_", 1:n_samples)
    
    # OTU ids 
    otu_ids <- paste0("OTU_", 1:n_otus)
    
    # Create an abundance matrix with some realistic patterns 
    set.seed(123)
    abundance_matrix <- matrix(0, nrow = n_otus, ncol = n_samples)
    
    # Add some structure to the data 
    for (i in 1:n_otus) {
      if (i <= 50) { # First 50 OTUs more abundant in first 10 samples
        abundance_matrix[i, 1:10] <- rpois(10, lambda = sample(10:50, 1))
        abundance_matrix[i, 11:30] <- rpois(20, lambda = sample(1:5, 1))
      } else if (i <= 100) {  # Next 50 OTUs more abundant in next 10 samples
        abundance_matrix[i, 11:20] <- rpois(10, lambda = sample(10:50, 1))
        abundance_matrix[i, c(1:10, 21:30)] <- rpois(20, lambda = sample(1:5, 1))
      } else if (i <= 150) {  # Next 50 OTUs more abundant in last 10 samples
        abundance_matrix[i, 21:30] <- rpois(10, lambda = sample(10:50, 1))
        abundance_matrix[i, 1:20] <- rpois(20, lambda = sample(1:5, 1))
      } else {  # Remaining OTUs randomly distributed
        abundance_matrix[i, ] <- rpois(n_samples, lambda = sample(1:10, 1))
      }
    }
    
    colnames(abundance_matrix) <- sample_names
    rownames(abundance_matrix) <- otu_ids
    
    # Create taxonomy table 
    tax_levels <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")
    
    # Define some common bacterial taxa for realism 
    kingdoms <- c("Bacteria", "Archaea")
    phyla <- c("Proteobacteria", "Firmicutes", "Bacteroidetes", "Actinobacteria", "Cyanobacteria")
    classes <- c("Alphaproteobacteria", "Bacilli", "Bacteroidia", "Actinobacteria", "Clostridia")
    orders <- c("Rhizobiales", "Lactobacillales", "Bacteroidales", "Bifidobacteriales", "Clostridiales")
    families <- c("Rhizobiaceae", "Lactobacillaceae", "Bacteroidaceae", "Bifidobacteriaceae", "Clostridiaceae")
    genera <- c("Rhizobium", "Lactobacillus", "Bacteroides", "Bifidobacterium", "Clostridium")
    
    # Create taxonomy matrix 
    taxonomy_matrix <- matrix("", nrow = n_otus, ncol = length(tax_levels))
    
    for (i in 1:n_otus) {
      phylum_idx <- sample(1:length(phyla), 1)
      class_idx <- sample(1:length(classes), 1)
      order_idx <- sample(1:length(orders), 1)
      family_idx <- sample(1:length(families), 1)
      genus_idx <- sample(1:length(genera), 1)
      
      taxonomy_matrix[i, ] <- c(
        sample(kingdoms, 1), 
        phyla[phylum_idx], 
        classes[class_idx],
        orders[order_idx], 
        families[family_idx], 
        genera[genus_idx],
        paste0(genera[genus_idx], " sp.", sample(1:20, 1))
      )
    }
    
    colnames(taxonomy_matrix) <- tax_levels
    rownames(taxonomy_matrix) <- otu_ids
    
    # Create sample metadata 
    body_sites <- c("gut", "skin", "oral", "nasal")
    subjects <- c("subject1", "subject2", "subject3", "subject4", "subject5")
    
    # Distribute samples evenly among body sites to mimic Moving Pictures tutorial
    site_assignments <- rep(body_sites, length.out = n_samples)
    
    sample_metadata <- data.frame(
      SampleID = sample_names, 
      BodySite = site_assignments, 
      Subject = sample(subjects, n_samples, replace = TRUE), 
      DaySinceExperimentStart = sample(0:30, n_samples, replace = TRUE),
      AntibioticUsage = sample(c("Yes", "No"), n_samples, replace = TRUE, prob = c(0.3, 0.7))
    )
    
    rownames(sample_metadata) <- sample_metadata$SampleID
    
    # Convert to phyloseq format 
    otu_table <- otu_table(abundance_matrix, taxa_are_rows = TRUE)
    tax_table <- tax_table(taxonomy_matrix)
    sample_data <- sample_data(sample_metadata)
    
    # Create phyloseq object 
    ps <- phyloseq(otu_table, tax_table, sample_data)
    
    return(list(phyloseq = ps, metadata = sample_metadata))
  }
  
  # Handle generate example data button
  observeEvent(input$generate_example, {
    # Generate example data using the existing function
    example_data <- generate_example_data()
    
    # Set the reactive values 
    phyloseq_obj(example_data$phyloseq)
    metadata(example_data$metadata)
    
    # Update select inputs with available metadata columns 
    meta_columns <- colnames(example_data$metadata)[-1]
    updateSelectInput(session, "group_var", 
                      choices = c("None", meta_columns), 
                      selected = "BodySite")
    updateSelectInput(session, "alpha_group", 
                      choices = c("None", meta_columns), 
                      selected = "BodySite")
    updateSelectInput(session, "beta_group", 
                      choices = c("None", meta_columns), 
                      selected = "BodySite")
    updateSelectInput(session, "shape_by", 
                      choices = c("None", meta_columns), 
                      selected = "None")
    updateSelectInput(session, "meta_x",
                      choices = c("None", meta_columns),
                      selected = meta_columns[1])
    updateSelectInput(session, "meta_y", 
                      choices = c("None", meta_columns),
                      selected = meta_columns[2])
    updateSelectInput(session, "meta_color", 
                      choices = c("None", meta_columns), 
                      selected = "BodySite")
    showNotification("Example data loaded automatically on startup!", type = "message")
  })
  
  # Handle data loading 
  observeEvent(input$load_data, {
    if (input$load_example) {
      # Load example data based on selected dataset
      if(exists("input$example_dataset") && input$example_dataset == "Soil Microbiome") {
        # Generate different example dataset if soil is selected 
        example_data <- generate_example_data()
      } else {
        # Default to human gut / moving pictures style data 
        example_data <- generate_example_data()
      }
      
      phyloseq_obj(example_data$phyloseq)
      metadata(example_data$metadata)
      
      # Update selectInputs with available metadata columns 
      meta_columns <- colnames(example_data$metadata)[-1] # Exclude SampleID
      updateSelectInput(session, "group_var", 
                        choices = c("None", meta_columns), 
                        selected = "BodySite")
      updateSelectInput(session, "alpha_group", 
                        choices = c("None", meta_columns), 
                        selected = "BodySite")
      updateSelectInput(session, "beta_group", 
                        choices = c("None", meta_columns), 
                        selected = "BodySite")
      updateSelectInput(session, "shape_by", 
                        choices = c("None", meta_columns), 
                        selected = "None")
      updateSelectInput(session, "meta_x",
                        choices = c("None", meta_columns),
                        selected = meta_columns[1])
      updateSelectInput(session, "meta_y", 
                        choices = c("None", meta_columns),
                        selected = meta_columns[2])
      updateSelectInput(session, "meta_color", 
                        choices = c("None", meta_columns), 
                        selected = "BodySite")
      showNotification("Example data loaded successfully!", type = "message")
    } else {
      # Code to load real data would go here 
      # For now, placeholder with an error message 
      showNotification("Real data loading not implemented yet. Please use the example data option.", 
                       type = "error", duration = 10)
    }
  })
  
  # Generate data summary 
  output$data_summary <- renderPrint({
    ps <- phyloseq_obj()
    if (is.null(ps)) {
      return("No data loaded. Please upload data files or load example data.")
    }
    
    cat("Microbiome Dataset Summary:\n")
    cat("Number of samples:", nsamples(ps), "\n")
    cat("Number of taxa:", ntaxa(ps), "\n")
    cat("Total reads:", sum(sample_sums(ps)), "\n")
    cat("Mean reads per sample:", mean(sample_sums(ps)), "\n")
    cat("Median reads per sample:", median(sample_sums(ps)), "\n")
    cat("\nTop phyla by abundance:\n")
    top_phyla <- tax_glom(ps, taxrank = "Phylum") %>%
      transform_sample_counts(function(x) x / sum(x) * 100) %>%
      psmelt() %>%
      group_by(Phylum) %>%
      summarize(MeanRelativeAbundance = mean(Abundance)) %>%
      arrange(desc(MeanRelativeAbundance)) %>%
      head(5)
    print(top_phyla)
    cat("\nSample distribution by body site:\n")
    sample_counts <- table(sample_data(ps)$BodySite)
    print(sample_counts)
  })
  
  # Preview metadata 
  output$metadata_preview <- renderDataTable({
    md <- metadata()
    if (is.null(md)) {
      return(NULL)
    }
    datatable(md, options = list(pageLength = 5, scrollX = TRUE))
  })
  
  # Taxonomy plot 
  output$taxonomy_plot <- renderPlot({
    ps <- phyloseq_obj()
    if (is.null(ps)) {
      return(NULL)
    }
    
    # Get selected taxonomic level 
    tax_level <- input$tax_level
    
    # Agglomerate at the selected taxonomic level 
    ps_glom <- tax_glom(ps, taxrank = tax_level)
    
    # Transform to relavtive abundance
    ps_rel <- transform_sample_counts(ps_glom, function(x) x / sum(x) * 100)
    
    # Melt to long format for ggplot 
    ps_melt <- psmelt(ps_rel)
    
    # Filter low abundance taxa 
    ps_melt <- ps_melt %>%
      group_by(get(tax_level)) %>%
      mutate(MeanAbundance = mean(Abundance)) %>%
      ungroup()
    
    low_abundance_taxa <- ps_melt %>%
      filter(MeanAbundance < input$min_abundance) %>%
      pull(!!sym(tax_level)) %>%
      unique()
    
    if (length(low_abundance_taxa) > 0) {
      ps_melt <- ps_melt %>%
        mutate(!!tax_level := ifelse(get(tax_level) %in% low_abundance_taxa, "Other", get(tax_level)))
    }
    
    # Create plot based on plot type
    if (input$plot_type == "bar") {
      # Group by the selected variable if not "None"
      if (input$group_var != "None") {
        group_var <- input$group_var
        
        # Reorder samples by group and abundance if requested 
        if (input$sort_by_abundance) {
          sample_order <- ps_melt %>%
            group_by(Sample, !!sym(group_var)) %>%
            summarize(TotalAbundance = sum(Abundance), .groups = "drop") %>%
            arrange(!!sym(group_var), desc(TotalAbundance)) %>%
            pull(Sample)
          ps_melt$Sample <- factor(ps_melt$Sample, levels = sample_order)
        }
        
        # Create the grouped stacked bar plot 
        p <- ggplot(ps_melt, aes(x = Sample, y = Abundance, fill = get(tax_level))) +
          geom_bar(stat = "identity") +
          labs(x = "Sample", y = "Relative Abundance (%)", fill = "tax_level") +
          theme_minimal() +
          theme(
            axis.text.x = element_text(angle = 45, hjust = 1, size = 8), 
            legend.position = if (input$show_legend) "right" else "none"
          ) +
          facet_grid(. ~ get(group_var), scales = "free_x", space = "free_x")
      } else {
        # Create the non-grouped stacked bar plot. Reorder samples by abundance if requested
        if (input$sort_by_abundance) {
          sample_order <- ps_melt %>%
            group_by(Sample) %>%
            summarize(TotalAbundance = sum(Abundance), .groups = "drop") %>%
            arrange(desc(TotalAbundance)) %>%
            pull(Sample)
          ps_melt$Sample <- factor(ps_melt$Sample, levels = sample_order)
        }
        
        p <- ggplot(ps_melt, aes(x = Sample, y = Abundance, fill = get(tax_level))) +
          geom_bar(stat = "identity") +
          labs(x = "Sample", y = "Relative Abundance (%)", fill = tax_level) +
          theme_minimal() +
          theme(
            axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
            legend.position = if(input$show_legend) "right" else "none"
          )
      }
      
      # Use a colorful palette 
      if (length(unique(ps_melt[[tax_level]])) <= 9) {
        p <- p + scale_fill_brewer(palette = "Set1")
      } else {
        p <- p + scale_fill_viridis_d()
      }
    } else if (input$plot_type == "heatmap") {
      # Create heatmap for taxonomic abundance. Summarize data for heatmap
      if (input$group_var != "None") {
        # Group samples for heatmap
        heat_data <- ps_melt %>%
          group_by(get(tax_level), !!sym(input$group_var)) %>%
          summarize(MeanAbundance = mean(Abundance), .groups = "drop") %>%
          spread(key = !!sym(input$group_var), value = MeanAbundance, fill = 0)
        
        # Convert back to long format for ggplot 
        heat_data_long <- gather(heat_data, key = "Group", value = "Abundance", -1)
        colnames(heat_data_long)[1] <- tax_level
        
        # Order taxa by overall abundance 
        taxa_order <- ps_melt %>%
          group_by(!!sym(tax_level)) %>%
          summarize(MeanAbundance = mean(Abundance), .groups = "drop") %>%
          arrange(desc(MeanAbundance)) %>%
          pull(!!sym(tax_level))
        heat_data_long[[tax_level]] <- factor(heat_data_long[[tax_level]], levels = taxa_order)
        
        # Create heatmap
        p <- ggplot(heat_data_long, aes(x = Group, y =!!sym(tax_level), fill = Abundance)) +
          geom_tile() +
          scale_fill_viridis_c(option = "plasma") +
          labs(x = input$group_var, y = tax_level, fill = "Mean\nAbundance (%)") +
          theme_minimal() +
          theme(
            axis.text.x = element_text(angle = 45, hjust = 1),
            legend.position = "right"
          )
      } else {
        # Sample-level heatmap. Order taxa and samples
        taxa_order <- ps_melt %>%
          group_by(!!sym(tax_level)) %>%
          summarize(MeanAbundance = mean(Abundance), .groups = "drop") %>%
          arrange(desc(MeanAbundance)) %>%
          head(20) %>% # Limit to top taxa for readability 
          pull(!!sym(tax_level))
        
        # Filter for top taxa 
        ps_melt_filtered <- ps_melt %>%
          filter(!!sym(tax_level) %in% taxa_order)
        
        ps_melt_filtered[[tax_level]] <- factor(ps_melt_filtered[[tax_level]], levels = taxa_order)
        
        # Create heatmap
        p <- ggplot(ps_melt_filtered, aes(x = Sample, y = !!sym(tax_level), fill = Abundance)) +
          geom_tile() +
          scale_fill_viridis_c(option = "plasma") +
          labs(x = "Sample", y = tax_level, fill = "Abundance (%)") +
          theme_minimal() +
          theme(
            axis.text.x = element_text(angle = 90, hjust = 1, size = 7),
            legend.position = "right"
          )
      }
    }
    return(p)
  })
  
  # Top taxa table
  output$taxa_table <- renderDataTable({
    ps <- phyloseq_obj()
    if (is.null(ps)) {
      return(NULL)
    }
    # Get selected taxonomic level 
    tax_level <- input$tax_level
    
    # Agglomerate and calculate relative abundance 
    ps_glom <- tax_glom(ps, taxrank = tax_level)
    ps_rel <- transform_sample_counts(ps_glom, function(x) x / sum(x) * 100)
    ps_melt <- psmelt(ps_rel)
    
    # Summarize by selected taxonomic level 
    taxa_summary <- ps_melt %>%
      group_by(!!sym(tax_level)) %>%
      summarize(
        MeanRelativeAbundance = mean(Abundance),
        MaxRelativeAbundance = max(Abundance), 
        Prevalence = sum(Abundance > 0) / length(unique(Sample)) * 100, 
        .groups = "drop"
      ) %>%
      arrange(desc(MeanRelativeAbundance))
    
    datatable(taxa_summary, 
              options = list(pageLength = 10, scrollX = TRUE),
              rownames = FALSE) %>%
      formatRound(columns = c("MeanRelativeAbundance", "MaxRelativeAbundance", "Prevalence"), digits = 2)
  })
  
  # Alpha diversity plot interactive
  output$alpha_plot <- renderPlotly({
    ps <-  phyloseq_obj()
    if (is.null(ps)) {
      return(NULL)
    }
    
    # Calculate alpha diversity 
    alpha_metric <- input$alpha_metric 
    if (!alpha_metric %in% c("Faith PD")) { # Standandard metrics 
      alpha_div <- estimate_richness(ps, measures = alpha_metric)
      alpha_div$SampleID <- rownames(alpha_div)
      
      # Merge with metadata 
      md <- as(sample_data(ps), "data.frame")
      alpha_data <- merge(alpha_div, md, by.x = "SampleID", by.y = "row.names")
      
      # Set up plotting 
      if (input$alpha_group != "None") {
        group_var <- input$alpha_group
        
        # Create the alpha diversity plot
        p <- ggplot(alpha_data, aes(x = !!sym(group_var), y = !!sym(alpha_metric), 
                                    color = !!sym(group_var),
                                    text = paste("Sample:", SampleID, 
                                                 "<br>", alpha_metric, ":", round(!!sym(alpha_metric), 2),
                                                 "<br>Group:", !!sym(group_var)))) +
          labs(x = group_var, y = paste(alpha_metric, "Diversity")) +
          theme_minimal() +
          theme(legend.position = "none")
        
        # Add points if requested 
        if(input$add_points) {
          p <- p + geom_jitter(width = 0.2, alpha = 0.7, size = 3)
        }
        
        # Add boxplot if requested 
        if (input$add_boxplot) {
          p <- p + geom_boxplot(alpha = 0.5, outlier.shape = NA)
        }
        
        # Add statistical test if requested 
        if (input$stat_test != "None" && length(unique(alpha_data[[group_var]])) > 1) {
          test_result <- NULL
          p_value <- NA
          
          #### 
          if (input$stat_test == "Krustal-Wallis") {
            # Run Kruskal-Wallis test
            test_formula <- as.formula(paste(alpha_metric, "~", group_var))
            test_result <- kruskal.test(test_formula, data = alpha_data)
            p_value <- test_result$p.value
          } else if (input$stat_test == "ANOVA") {
            # Run ANOVA 
            test_formula <- as.formula(paste(alpha_metric, "~", group_var))
            test_result <- aov(test_formula, data = alpha_data)
            p_value <- summary(test_result)[[1]]$"Pr(>F)"[1]
          } else if (input$stat_test %in% c("T-test", "Wilcoxon") && 
                     length(unique(alpha_data[[group_var]])) == 2) {
            # Run t-test or Wilcoxon test for two groups
            group_levels <- unique(alpha_data[[group_var]])
            group1_data <- alpha_data[alpha_data[[group_var]] == group_levels[1], alpha_metric]
            group2_data <- alpha_data[alpha_data[[group_var]] == group_levels[2], alpha_metric]
            
            if (input$stat_test == "T-test") {
              test_result <- t.test(group1_data, group2_data)
            } else {
              test_result <- wilcox.test(group1_data, group2_data)
            }
            p_value <- test_result$p.value
          }
          
          # Add p-value annotation if test was run 
          if (!is.na(p_value)) {
            p_value_formatted <- format(p_value, digits =3)
            if (p_value < 0.001) p_value_formatted <- "p < 0.001"
            
            p <- p + labs(subtitle = paste(input$stat_test, "p-value:", p_value_formatted))
            
            # Store test result for detailes output 
            output$alpha_stats <- renderPrint({
              print(test_result)
            })
          }
        }
      } else {
        # Simple histogram if no grouping 
        p <- ggplot(alpha_data, aes(x = !!sym(alpha_metric),
                                    text = paste("Value:", round(!!sym(alpha_metric), 2), 
                                                 "<br>Sample:", SampleID))) +
          geom_histogram(fill = "steelblue", color ="black", bins = 15) +
          labs(x = paste(alpha_metric, "Diversity"), y = "Count") +
          theme_minimal()
      }
    } else {
      # Placeholder for Faith's PD (would require phylogenetic tree)
      p <- ggplot() +
        annotate("text", x = 0.5, y = 0.5, 
                 label = "Faith's PD requires a phylogenetic tree.\nNot available in example data.") +
        theme_void()
    } ##############
    # Convert ggplot to plotly
    p_interactive <- ggplotly(p, tooltip = "text") %>%
      layout(hoverlabel = list(bgcolor ="white", font = list(size = 12)), 
             dragmode = "zoom") %>%
      config(displayModeBar = TRUE, scrollZoom = TRUE)
    return(p_interactive)
  })
  
  # Rarefaction curve 
  output$rarefaction_plot <- renderPlotly({
    ps <- phyloseq_obj() 
    if (is.null(ps)) {
      return(NULL)
    }
    
    # Get OTU table and ensure samples are rows 
    otu <- as(otu_table(ps), "matrix")
    if (taxa_are_rows(ps)) {
      otu <- t(otu) 
    }
    
    # Check if we have data to work with
    if (nrow(otu) == 0 || ncol(otu) == 0) {
      return(NULL)
    }
    
    # Calculate rarefaction 
    raremax <- min(rowSums(otu))
    step_size <- max(floor(raremax/20), 1) # Ensure step size is at least 1
    
    # Get rarefaction data using vegan tidy approach 
    rarecurve_data <- vegan::rarecurve(otu, step = step_size, tidy = TRUE)
    
    # Rename columns for clarity 
    colnames(rarecurve_data) <- c("SampleName", "Reads", "Species")
    
    # Add grouping information if requested 
    if (input$alpha_group != "None") {
      # Get metadata 
      md <- as(sample_data(ps), "data.frame")
      md$SampleName <- rownames(md)
      
      # Merge with metadata 
      rarecurve_data$SampleName <- as.character(rarecurve_data$SampleName)
      rarecurve_data <- merge(rarecurve_data, md, by = "SampleName")
      
      # Check if grouping variable exists in the metadata
      if (!input$alpha_group %in% colnames(rarecurve_data)) {
        # Create plot with grouping
        p <- ggplot(rarecurve_data, aes(x = Reads, y = Species, group = SampleName,
                                        color = .data[[input$alpha_group]], 
                                        text = paste("Sample:", SampleName, 
                                                     "<br>Reads:", Reads, 
                                                     "<br>ASVs/OTUs:", Species))) +
          geom_line(linewidth = 0.5, alpha = 0.7) +
          labs(x = "Sequencing Depth", y = "Observed ASVs/OTUs",
               title = "Rarefaction Curve") +
          theme_minimal() +
          theme(legend.position = "right")
      } else {
        # Fallback without grouping 
        p <- ggplot(rarecurve_data, aes(x = Reads, y = Species, group = SampleName, 
                                        color = SampleName, 
                                        text = paste("Sample:", SampleName, 
                                                     "<br>Reads:", Reads, 
                                                     "<br>ASVs/OTUs:", Species))) +
          geom_line(linewidth = 0.5, alpha = 0.7) +
          labs(x = "Sequencing Depth", y = "Observed ASVs/OTUs",
               title = "Rarefaction Curve") +
          theme_minimal()
      }
    } else {
      # Create plot without grouping 
      p <- ggplot(rarecurve_data, aes(x = Reads, y = Species, group = SampleName,
                                      color = SampleName, 
                                      text = paste("Sample:", SampleName,
                                                   "<br>Reads:", Reads, 
                                                   "<br>ASVs/OTUs:", Species))) +
        geom_line(linewidth = 0.5, alpha = 0.7) +
        labs(x = "Sequencing Depth", y = "Observed ASVs/OTUs",
             title = "Rarefaction Curve") +
        theme_minimal() +
        scale_colour_viridis_d(option = "plasma")
    }
    
    # Convert to plotly for interactivity 
    p_interactive <- ggplotly(p, tooltip = "text")
    
    # Improve layout 
    p_interactive <- p_interactive %>%
      layout(hoverlabel = list(bgcolor = "white"), 
             legend = list(title = list(text = input$alpha_group))) %>%
      config(displayModeBar = TRUE, scrollZoom = TRUE)
    
    return(p_interactive)
  })
  
  # Beta diversity plot
  output$beta_plot <- renderPlotly({
    ps <- phyloseq_obj()
    if (is.null(ps)) {
      return(NULL)
    }
    # Select distance metric
    dist_method <- switch(input$beta_metric,
                          "Bray-Curtis" = "bray",
                          "Jaccard" = "jaccard",
                          "UniFrac" = "unifrac",
                          "Weighted UniFrac" = "wunifrac",
                          "bray") # Default to Bray-Curtis
    # Skip UniFrac methods if no tree is available
    if (dist_method %in% c("unifrac", "wunifrac") && is.null(phy_tree(ps, errorIfNull = FALSE))) {
       p <- ggplot() +
               annotate("text", x = 0.5, y = 0.5,
                        label = "UniFrac requires a phylogenetic tree.\nNot available in example data.") +
               theme_void()
       return(ggplotly(p))
    }
    # Calculate distance matrix
    dist_matrix <- phyloseq::distance(ps, method = dist_method)
    # Perform ordination
    ord_method <- switch(input$ordination,
                         "PCoA" = "PCoA",
                         "NMDS" = "NMDS",
                         "t-SNE" = "t-SNE",
                         "UMAP" = "UMAP",
                         "PCoA") # Default to PCoA
    # Extract metadata
    metadata_df <- as(sample_data(ps), "data.frame")
    
    # Handle different ordination methods 
    if (ord_method == "PCoA") {
      ord <- ordinate(ps, method = ord_method, distance = dist_matrix)
      ord_data <- plot_ordination(ps, ord, justDF = TRUE)
      axis_labels <- c(paste0("PCo1 [", round(ord$values$Relative_eig[1] * 100, 1), "%"),
                       paste0("PCo2 [", round(ord$values$Relative_eig[2] * 100, 1), "%"))
    } else if (ord_method == "NMDS") {
      ord <- ordinate(ps, method = ord_method, distance = dist_matrix)
      ord_data <- plot_ordination(ps, ord, justDF = TRUE)
      axis_labels <- c("NMDS1", "NMDS2")
    } else if (ord_method == "t-SNE") {
      tsne_result <- Rtsne::Rtsne(dist_matrix, is_distance = TRUE, perplexity = min(30, nrow(metadata_df) - 1))
      ord_data <- data.frame(
        Axis.1 = tsne_result$Y[, 1], 
        Axis.2 = tsne_result$Y[, 2], 
        row.names = rownames(metadata_df)
      )
      # Add sample data 
      ord_data <- cbind(ord_data, metadata_df)
      axis_labels <- c("t-SNE1", "t-SNE2")
    } else if (ord_method == "UMAP") {
      umap_config <- umap::umap.defaults
      umap_config$n_neighbors <- min(15, nrow(metadata_df) - 1)
      umap_result <- umap::umap(dist_matrix, config = umap_config)
      ord_data <- data.frame(
        Axis.1 <- umap_result$layout[, 1], 
        Axis.2 <- umap_result$layout[, 2], 
        row.names = rownames(metadata_df)
      )
      ord_data <- cbind(ord_data, metadata_df)
      axis_labels <- c("UMAP1", "UMAP2")
    }
    
    # Check if display boxplot or ordination plot
    if (input$beta_group != "None" && input$plot_type == "boxplot") {
      # Create distance to centroid data for boxplots 
      disp <- vegan::betadisper(dist_matrix, group = metadata_df[[input$beta_group]])
      
      # Extract distance to centroid data 
      boxplot_data <- data.frame(
        Group = metadata_df[[input$beta_group]], 
        Distance = dist$distances, 
        Sample = rownames(metadata_df)
      )
      
      # Create boxplot 
      p <- ggplot(boxplot_data, aes(x = Group, y = Distance, fill = Group)) +
        geom_boxplot(alpha = 0.7, outlier.shape = NA) +
        geom_jitter(aes(text = Sample), width = 0.2, alpha = 0.7) +
        theme_minimal() +
        labs(
          x = input$beta_group, 
          y = paste("Distance to centroid -", input$beta_metric), 
          title = paste("Beta dispersion using", input$beta_metric, "distances")
        )
      
      # Run PERMANOVA if requested 
      if (input$permanova == "Yes") {
        # Run PERMANOVA test
        perm_formula <- as.formula(paste("dist_matrix ~", input$beta_group))
        permanova_result <- vegan::adonis2(perm_formula, data = metadata_df)
        
        # Display results 
        output$permanova_results <- renderPrint({
          cat("PERMANOVA Results for", input$beta_metric, "distance by", input$beta_group, "\n\n")
          print(permanova_result)
          # Calculate and display variance explained 
          r2 <- permanova_result$R2[1]
          cat("\nVariance explained:", round(r2 * 100, 2), "%\n")
          
          # Add BETADISPER results (test of homogeneity of dispersions)
          cat("\nBETADISPER Results (Homogeneity of Dispersions Test):\n")
          anova_result <- anova(disp)
          print(anova_result)
        })
        
        # Add R² and p-value to plot title 
        r2 <- permanova_result$R2[1]
        p_val <- permanova_result$`Pr(>F)`[1]
        p_val_text <- ifelse(p_val < 0.001, "p < 0.001", paste("p =", round(p_val, 3)))
        
        # Add betadisper p-value 
        betadisp_pval <- anova(disp)$`Pr(>F)`[1]
        betadisp_text <- ifelse(betadisp_pval < 0.001, "p < 0.001", paste("p =", round(betadisp_pval, 3)))
        
        p <- p + labs(subtitle = paste0("PERMANOVA: R² =", round(r2, 3), ", ", p_val_text, 
                                        "\nBETADISPER: ", betadisp_text))
      }
      # Make the boxplot interactive 
      int <- ggplotly(p, tooltip = c("x", "y", "text"))
      return(int)
    } else {
      # Regular ordination plot. Set up plot aesthetics 
      aes_list <- list(x = sym("Axis.1"), y = sym("Axis.2"))
      
      # Add color aesthetic if grouping is provided 
      if (input$beta_group != "None") {
        aes_list$color <- sym(input$beta_group)
      }
      # Add shape aesthetic if requested 
      if (input$shape_by != "None" && input$shape_by != input$beta_group) {
        aes_list$shape <- sym(input$shape_by)
      }
      # Create tooltip text
      if (input$beta_group != "None" && input$shape_by != "None" && input$shape_by != input$beta_group) {
        ord_data$tooltip <- paste("Sample:", rownames(ord_data), 
                                  "\n", input$beta_group, ":", ord_data[[input$beta_group]],
                                  "\n", input$shape_by, ":", ord_data[[input$shape_by]])
      } else if (input$beta_group != "None") {
        ord_data$tooltip <- paste("Sample:", rownames(ord_data), 
                                  "\n", input$beta_group, ":", ord_data[[input$beta_group]])
      } else {
        ord_data$tooltip <- paste("Sample:", rownames(ord_data))
      }
      
      aes_list$text <- sym("tooltip")
      
      # Create the base plot
      p <- ggplot(ord_data, do.call(aes, aes_list)) +
        geom_point(size = 3, alpha = 0.7) +
        theme_minimal() +
        labs(
          x = axis_labels[1],
          y = axis_labels[2],
          title = paste(input$beta_metric, "Distance -", ord_method, "Ordination")
        )
      
      # Add ellipses if requested and grouping is provided
      if (input$add_ellipse && input$beta_group != "None" && input$plot_type == "ellipse") {
        p <- p + stat_ellipse(aes(color = !!sym(input$beta_group)), type = "norm", level = 0.95)
      }
      
      # Add sample labels if requested
      if (input$add_labels) {
        p <- p + geom_text(aes(label = rownames(ord_data)), size = 3, vjust = -1, check_overlap = TRUE)
      }
      
      # Run PERMANOVA if requested
      if (input$beta_group != "None" && input$permanova == "Yes") {
        # Run PERMANOVA test
        perm_formula <- as.formula(paste("dist_matrix ~", input$beta_group))
        permanova_result <- vegan::adonis2(perm_formula, data = metadata_df)
        
        # Display results
        output$permanova_results <- renderPrint({
          cat("PERMANOVA Results for", input$beta_metric, "distance by", input$beta_group, "\n\n")
          print(permanova_result)
          # Calculate and display variance explained
          r2 <- permanova_result$R2[1]
          cat("\nVariance explained:", round(r2 * 100, 2), "%\n")
        })
        
        # Add R² and p-value to plot title
        r2 <- permanova_result$R2[1]
        p_val <- permanova_result$`Pr(>F)`[1]
        p_val_text <- ifelse(p_val < 0.001, "p < 0.001", paste("p =", round(p_val, 3)))
        p <- p + labs(subtitle = paste0("PERMANOVA: R² = ", round(r2, 3), ", ", p_val_text))
      }
      # Convert to plotly
      plt <- ggplotly(p, tooltip = "text")
      
      # Customize plotly layout
      plt <- plt %>% layout(
        hoverlabel = list(bgcolor = "white", font = list(size = 12)),
        legend = list(title = list(text = input$beta_group))
      )
      return(plt)
    }
  })
  
  # Download handler for beta diversity plot
  output$download_beta_plot <- downloadHandler(
    filename = function() {
      paste("beta_diversity_", input$beta_metric, "_", input$ordination, "_", Sys.Date(), ".png", sep = "")
    },
    content = function(file) {
      # Get the plotly object
      p <- plotly_build(output$beta_plot)
      
      # Save as a PNG file
      png(file, width = 1000, height = 800, res = 100)
      print(ggplotly(p))
      dev.off()
    }
  )
  
  
  ###### AGB RESULTS DELIVERY #########
  
  
  # Load required libraries
  library(shiny)
  library(shinydashboard)
  library(ggplot2)
  library(phyloseq)
  library(tidyverse)
  library(vegan)
  library(readr)
  library(DT)
  library(ape)
  library(bslib)
  library(plotly)
  library(reshape2)
  library(Rtsne)
  library(umap)
  
  # Define UI
  ui <- dashboardPage(
    skin = "blue", 
    
    # Header from dashboard style but with navbar styling 
    dashboardHeader(
      title = div(icon("microscope"), "AGB Microbiome Results"), 
      titleWidth = 300
    ), 
    
    # Sidebar with dashboard style 
    dashboardSidebar(
      width = 300, 
      sidebarMenu(
        menuItem(" Home", tabName = "home", icon = icon("home")), 
        menuItem(" Data Upload", tabName = "data", icon = icon("upload")),
        menuItem(" Taxonomy", tabName = "taxonomy", icon = icon("bacteria")),
        menuItem(" Alpha Diversity", tabName = "alpha", icon = icon("chart-line")),
        menuItem(" Beta Diversity", tabName = "beta", icon = icon("project-diagram")),
        # menuItem("️ Functional Pathways", tabName = "function", icon = icon("dna")),
        menuItem(" Metadata Explorer", tabName = "metadata", icon = icon("table")),
        menuItem(" Reports", tabName = "reports", icon = icon("file-export"))
      )
    ),
    
    # Body with dashboard style 
    dashboardBody(
      # Use BS theme for consistent styling 
      tags$head(
        tags$link(rel = "stylesheet", type = "text/css", href = "https://cdnjs.cloudflare.com/ajax/libs/font-awesome/5.15.4/css/all.min.css"), 
        tags$style(HTML("
                      .box.box-solid.box-primary>.box-header {
                        background: #3c8dbc;
                        color: #fff;
                      }
                      .nav-tabs-custom>.nav-tabs>li.active {
                        border-top-color: #3c8dbc;
                      }
                      .content-wrapper, .right-side {
                        background-color: #f4f6f9;
                      }
                      .box {
                        border-radius: 3px;
                        box-shadow: 0 2px 5px rgba(0,0,0,0.1);
                      }
                      .well {
                        background-color: #fff;
                        border: 1px solid #e3e3e3;
                        border-radius: 3px;
                        box-shadow: 0 1px 3px rgba(0,0,0,0.05);
                      }
                      "))
      ), 
      
      tabItems(
        # Home tab 
        tabItem(tabName = "home", 
                fluidRow(
                  box(
                    title = "Welcome to AGB Results", 
                    status = "primary", 
                    solidHeader = TRUE, 
                    width = 12, 
                    div(
                      h3("Comprehensive Microbiome Analysis Dashboard"), 
                      hr(), 
                      fluidRow(
                        column(
                          width = 6, 
                          h4("🔍 Features"), 
                          tags$ul(
                            tags$li(strong("Taxonomic Profiles:"), "Explore the microbial composition at various taxonomic levels"),
                            tags$li(strong("Diversity Analysis:"), "Assess alpha and beta diversity metrics"), 
                            tags$li(strong("Metadata Exploration:"), "Correlate microbial features with sample metadata"), 
                            tags$li(strong("Report Generation:"), "Generate customized reports for your findings")
                          )
                        ), 
                        column(
                          width = 6, 
                          h4("🚀 Getting Started"), 
                          p("Begin by uploading your data files in the Data Upload tab"),
                          tags$ul(
                            tags$li("Upload ASV/OTU table, taxonomy assignments, and sample metadata"),
                            tags$li("Explore taxonomic composition in the Taxonomy tab"), 
                            tags$li("Analyze diversity metrics in the Alpha and Beta Diversity tabs"), 
                            tags$li("Generate reports with key findings")
                          )
                        )
                      )
                    )
                  )
                )
        ), 
        
        # Data Upload Tab 
        tabItem(tabName = "data",
                fluidRow(
                  box(
                    title = "Upload Microbiome Data", 
                    status = "primary", 
                    solidHeader = TRUE, 
                    width = 12, 
                    fluidRow(
                      column(
                        width = 6,
                        fileInput("asv_table", "Upload ASV/OTU Table (CSV or TSV format):",
                                  accept = c(".csv", ".tsv", ".txt")), 
                        fileInput("taxonomy_file", "Upload Taxonomy File (CSV or TSV format):",
                                  accept = c(".csv", ".tsv", ".txt")),
                        fileInput("metadata_file", "Upload Sample Metadata (CSV or TSV format):",
                                  accept = c(".csv", ".tsv", ".txt")),
                        actionButton("load_data", "Load Data", 
                                     icon = icon("upload"),
                                     class = "btn-success"),
                        actionButton("generate_example", "Generate Example Data", 
                                     icon = icon("database"),
                                     class = "btn-info")
                      ),
                      column(
                        width = 6,
                        checkboxInput("load_example", "Load Example Data Instead", value = FALSE),
                        conditionalPanel(
                          condition = "input.load_example == true",
                          selectInput("example_dataset", "Select Example Dataset:",
                                      choices = c("Human Gut Microbiome", "Soil Microbiome"))
                        ),
                        hr(),
                        h4("Data Preview:"),
                        dataTableOutput("metadata_preview")
                      )
                    )
                  )
                ), 
                fluidRow(
                  box(
                    title = "Data Summary", 
                    status = "info", 
                    solidHeader = TRUE,
                    width = 12,
                    verbatimTextOutput("data_summary")
                  )
                )
        ),
        
        # Taxonomy Tab 
        tabItem(tabName = "taxonomy",
                fluidRow(
                  box(
                    title = "Taxonomy Settings", 
                    status = "primary", 
                    solidHeader = TRUE,
                    width = 3,
                    selectInput("tax_level", "Taxonomic Level:",
                                choices = c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"),
                                selected = "Phylum"),
                    selectInput("group_var", "Group By:",
                                choices = c("None"), selected = "None"),
                    sliderInput("min_abundance", "Minimum Abundance (%):", 
                                min = 0, max = 5, value = 1, step = 0.1),
                    checkboxInput("sort_by_abundance", "Sort by Abundance", value = TRUE),
                    checkboxInput("show_legend", "Show Legend", value = TRUE),
                    radioButtons("plot_type", "Plot Type:",
                                 choices = c("Stacked Bar" = "bar", "Heatmap" = "heatmap"),
                                 selected = "bar"),
                    hr(),
                    downloadButton("download_tax_plot", "Download Plot", class = "btn-primary")
                  ),
                  box(
                    title = "Taxonomic Composition", 
                    status = "info", 
                    solidHeader = TRUE,
                    width = 9,
                    plotOutput("taxonomy_plot", height = 600)
                  )
                ),
                fluidRow(
                  box(
                    title = "Top Taxa Table", 
                    status = "info", 
                    solidHeader = TRUE,
                    width = 12,
                    dataTableOutput("taxa_table")
                  )
                )
        ),
        
        # Alpha Diversity Tab
        tabItem(tabName = "alpha",
                fluidRow(
                  box(
                    title = "Alpha Diversity Settings", 
                    status = "primary", 
                    solidHeader = TRUE,
                    width = 3,
                    selectInput("alpha_metric", "Alpha Diversity Metric:",
                                choices = c("Shannon", "Observed", "Simpson", "InvSimpson", "Faith PD", "Chao1"),
                                selected = "Shannon"),
                    selectInput("alpha_group", "Group Samples By:",
                                choices = c("None"), selected = "None"),
                    conditionalPanel(
                      condition = "input.alpha_group != 'None'",
                      selectInput("stat_test", "Statistical Test:",
                                  choices = c("None", "ANOVA", "Kruskal-Wallis", "T-test", "Wilcoxon"),
                                  selected = "None")
                    ),
                    sliderInput("rarefaction_depth", "Rarefaction Depth:", 
                                min = 1000, max = 10000, value = 5000, step = 1000),
                    checkboxInput("add_boxplot", "Add Box Plot", value = TRUE),
                    checkboxInput("add_points", "Show Individual Points", value = TRUE),
                    hr(),
                    downloadButton("download_alpha_plot", "Download Plot", class = "btn-primary")
                  ),
                  box(
                    title = "Alpha Diversity Plot", 
                    status = "info", 
                    solidHeader = TRUE,
                    width = 9,
                    plotlyOutput("alpha_plot", height = 500),
                    conditionalPanel(
                      condition = "input.stat_test != 'None'",
                      h4("Statistical Test Results:"),
                      verbatimTextOutput("alpha_stats")
                    )
                  )
                ),
                fluidRow(
                  box(
                    title = "Rarefaction Curve", 
                    status = "info", 
                    solidHeader = TRUE,
                    width = 12,
                    plotlyOutput("rarefaction_plot", height = 400)
                  )
                )
        ),
        
        # Beta Diversity Tab
        tabItem(tabName = "beta",
                fluidRow(
                  box(
                    title = "Beta Diversity Settings", 
                    status = "primary", 
                    solidHeader = TRUE,
                    width = 3,
                    selectInput("beta_metric", "Beta Diversity Metric:",
                                choices = c("Bray-Curtis", "Jaccard", "UniFrac", "Weighted UniFrac"),
                                selected = "Bray-Curtis"),
                    selectInput("ordination", "Ordination Method:",
                                choices = c("PCoA", "NMDS", "t-SNE", "UMAP"),
                                selected = "PCoA"),
                    selectInput("beta_group", "Color By:",
                                choices = c("None"), selected = "None"),
                    selectInput("shape_by", "Shape By:",
                                choices = c("None"), selected = "None"),
                    conditionalPanel(
                      condition = "input.beta_group != 'None'",
                      selectInput("permanova", "Run PERMANOVA:",
                                  choices = c("No", "Yes"),
                                  selected = "No")
                    ),
                    conditionalPanel(
                      condition ="input.beta_group != 'None'",
                      radioButtons("plot_type", "Visualization Type:",
                                   choices = c("Ellipses" = "ellipse", "Boxplots" = "boxplot"), 
                                   selected = "ellipse")
                    ),
                    conditionalPanel(
                      condition = "input.plot_type == 'ellipse' && input.beta_group != 'None'", 
                      checkboxInput("add_ellipse", "Add Confidence Ellipses", value = TRUE)
                    ),
                    checkboxInput("add_labels", "Add Sample Labels", value = FALSE),
                    hr(),
                    downloadButton("download_beta_plot", "Download Plot", class = "btn-primary")
                  ),
                  box(
                    title = "Beta Diversity Plot", 
                    status = "info", 
                    solidHeader = TRUE,
                    width = 9,
                    plotlyOutput("beta_plot", height = 500),
                    conditionalPanel(
                      condition = "input.permanova == 'Yes'",
                      h4("PERMANOVA Results:"),
                      verbatimTextOutput("permanova_results")
                    )
                  )
                )
        ),
        
        # Metadata Explorer Tab
        tabItem(tabName = "metadata",
                fluidRow(
                  box(
                    title = "Metadata Exploration", 
                    status = "primary", 
                    solidHeader = TRUE,
                    width = 3,
                    selectInput("meta_x", "X-axis Variable:",
                                choices = c("None"), selected = "None"),
                    selectInput("meta_y", "Y-axis Variable:",
                                choices = c("None"), selected = "None"),
                    selectInput("meta_color", "Color By:",
                                choices = c("None"), selected = "None"),
                    radioButtons("meta_plot_type", "Plot Type:",
                                 choices = c("Scatter" = "scatter", 
                                             "Box Plot" = "box", 
                                             "Bar Plot" = "bar"),
                                 selected = "scatter"),
                    hr(),
                    downloadButton("download_meta_plot", "Download Plot", class = "btn-primary")
                  ),
                  box(
                    title = "Metadata Visualization", 
                    status = "info", 
                    solidHeader = TRUE,
                    width = 9,
                    plotOutput("metadata_plot", height = 500)
                  )
                ),
                fluidRow(
                  box(
                    title = "Correlation Analysis", 
                    status = "info", 
                    solidHeader = TRUE,
                    width = 12,
                    plotOutput("correlation_plot", height = 400)
                  )
                )
        ),
        
        # Reports Tab
        tabItem(tabName = "reports",
                fluidRow(
                  box(
                    title = "Generate Reports", 
                    status = "primary", 
                    solidHeader = TRUE,
                    width = 12,
                    fluidRow(
                      column(
                        width = 6,
                        h4("Report Options"),
                        checkboxGroupInput("report_sections", "Include Sections:",
                                           choices = c("Taxonomy Analysis" = "taxonomy",
                                                       "Alpha Diversity" = "alpha",
                                                       "Beta Diversity" = "beta",
                                                       "Functional Analysis" = "function",
                                                       "Metadata Exploration" = "metadata"),
                                           selected = c("taxonomy", "alpha", "beta")),
                        radioButtons("report_format", "Report Format:",
                                     choices = c("PDF" = "pdf", 
                                                 "HTML" = "html",
                                                 "Word Document" = "docx"),
                                     selected = "html"),
                        textInput("report_title", "Report Title:", 
                                  value = "Microbiome Analysis Report"),
                        textAreaInput("report_comments", "Additional Comments:", 
                                      rows = 3),
                        hr(),
                        downloadButton("generate_report", "Generate Report", class = "btn-success")
                      ),
                      column(
                        width = 6,
                        h4("Report Preview"),
                        htmlOutput("report_preview")
                      )
                    )
                  )
                )
        )
      )
    )
    
  )
  
  
  # Server logic with mock data
  
  server <- function(input, output, session) {
    
    # Reactive values to store data 
    phyloseq_obj <- reactiveVal(NULL)
    metadata <- reactiveVal(NULL)
    
    # Function to generate example data similar to Moving Pictures tutorial
    generate_example_data <- function() {
      # Create example OTU table 
      n_samples <- 30
      n_otus <- 200
      
      # Sample names 
      sample_names <- paste0("Sample_", 1:n_samples)
      
      # OTU ids 
      otu_ids <- paste0("OTU_", 1:n_otus)
      
      # Create an abundance matrix with some realistic patterns 
      set.seed(123)
      abundance_matrix <- matrix(0, nrow = n_otus, ncol = n_samples)
      
      # Add some structure to the data 
      for (i in 1:n_otus) {
        if (i <= 50) { # First 50 OTUs more abundant in first 10 samples
          abundance_matrix[i, 1:10] <- rpois(10, lambda = sample(10:50, 1))
          abundance_matrix[i, 11:30] <- rpois(20, lambda = sample(1:5, 1))
        } else if (i <= 100) {  # Next 50 OTUs more abundant in next 10 samples
          abundance_matrix[i, 11:20] <- rpois(10, lambda = sample(10:50, 1))
          abundance_matrix[i, c(1:10, 21:30)] <- rpois(20, lambda = sample(1:5, 1))
        } else if (i <= 150) {  # Next 50 OTUs more abundant in last 10 samples
          abundance_matrix[i, 21:30] <- rpois(10, lambda = sample(10:50, 1))
          abundance_matrix[i, 1:20] <- rpois(20, lambda = sample(1:5, 1))
        } else {  # Remaining OTUs randomly distributed
          abundance_matrix[i, ] <- rpois(n_samples, lambda = sample(1:10, 1))
        }
      }
      
      colnames(abundance_matrix) <- sample_names
      rownames(abundance_matrix) <- otu_ids
      
      # Create taxonomy table 
      tax_levels <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")
      
      # Define some common bacterial taxa for realism 
      kingdoms <- c("Bacteria", "Archaea")
      phyla <- c("Proteobacteria", "Firmicutes", "Bacteroidetes", "Actinobacteria", "Cyanobacteria")
      classes <- c("Alphaproteobacteria", "Bacilli", "Bacteroidia", "Actinobacteria", "Clostridia")
      orders <- c("Rhizobiales", "Lactobacillales", "Bacteroidales", "Bifidobacteriales", "Clostridiales")
      families <- c("Rhizobiaceae", "Lactobacillaceae", "Bacteroidaceae", "Bifidobacteriaceae", "Clostridiaceae")
      genera <- c("Rhizobium", "Lactobacillus", "Bacteroides", "Bifidobacterium", "Clostridium")
      
      # Create taxonomy matrix 
      taxonomy_matrix <- matrix("", nrow = n_otus, ncol = length(tax_levels))
      
      for (i in 1:n_otus) {
        phylum_idx <- sample(1:length(phyla), 1)
        class_idx <- sample(1:length(classes), 1)
        order_idx <- sample(1:length(orders), 1)
        family_idx <- sample(1:length(families), 1)
        genus_idx <- sample(1:length(genera), 1)
        
        taxonomy_matrix[i, ] <- c(
          sample(kingdoms, 1), 
          phyla[phylum_idx], 
          classes[class_idx],
          orders[order_idx], 
          families[family_idx], 
          genera[genus_idx],
          paste0(genera[genus_idx], " sp.", sample(1:20, 1))
        )
      }
      
      colnames(taxonomy_matrix) <- tax_levels
      rownames(taxonomy_matrix) <- otu_ids
      
      # Create sample metadata 
      body_sites <- c("gut", "skin", "oral", "nasal")
      subjects <- c("subject1", "subject2", "subject3", "subject4", "subject5")
      
      # Distribute samples evenly among body sites to mimic Moving Pictures tutorial
      site_assignments <- rep(body_sites, length.out = n_samples)
      
      sample_metadata <- data.frame(
        SampleID = sample_names, 
        BodySite = site_assignments, 
        Subject = sample(subjects, n_samples, replace = TRUE), 
        DaySinceExperimentStart = sample(0:30, n_samples, replace = TRUE),
        AntibioticUsage = sample(c("Yes", "No"), n_samples, replace = TRUE, prob = c(0.3, 0.7))
      )
      
      rownames(sample_metadata) <- sample_metadata$SampleID
      
      # Convert to phyloseq format 
      otu_table <- otu_table(abundance_matrix, taxa_are_rows = TRUE)
      tax_table <- tax_table(taxonomy_matrix)
      sample_data <- sample_data(sample_metadata)
      
      # Create phyloseq object 
      ps <- phyloseq(otu_table, tax_table, sample_data)
      
      return(list(phyloseq = ps, metadata = sample_metadata))
    }
    
    # Handle generate example data button
    observeEvent(input$generate_example, {
      # Generate example data using the existing function
      example_data <- generate_example_data()
      
      # Set the reactive values 
      phyloseq_obj(example_data$phyloseq)
      metadata(example_data$metadata)
      
      # Update select inputs with available metadata columns 
      meta_columns <- colnames(example_data$metadata)[-1]
      updateSelectInput(session, "group_var", 
                        choices = c("None", meta_columns), 
                        selected = "BodySite")
      updateSelectInput(session, "alpha_group", 
                        choices = c("None", meta_columns), 
                        selected = "BodySite")
      updateSelectInput(session, "beta_group", 
                        choices = c("None", meta_columns), 
                        selected = "BodySite")
      updateSelectInput(session, "shape_by", 
                        choices = c("None", meta_columns), 
                        selected = "None")
      updateSelectInput(session, "meta_x",
                        choices = c("None", meta_columns),
                        selected = meta_columns[1])
      updateSelectInput(session, "meta_y", 
                        choices = c("None", meta_columns),
                        selected = meta_columns[2])
      updateSelectInput(session, "meta_color", 
                        choices = c("None", meta_columns), 
                        selected = "BodySite")
      showNotification("Example data loaded automatically on startup!", type = "message")
    })
    
    # Handle data loading 
    observeEvent(input$load_data, {
      if (input$load_example) {
        # Load example data based on selected dataset
        if(exists("input$example_dataset") && input$example_dataset == "Soil Microbiome") {
          # Generate different example dataset if soil is selected 
          example_data <- generate_example_data()
        } else {
          # Default to human gut / moving pictures style data 
          example_data <- generate_example_data()
        }
        
        phyloseq_obj(example_data$phyloseq)
        metadata(example_data$metadata)
        
        # Update selectInputs with available metadata columns 
        meta_columns <- colnames(example_data$metadata)[-1] # Exclude SampleID
        updateSelectInput(session, "group_var", 
                          choices = c("None", meta_columns), 
                          selected = "BodySite")
        updateSelectInput(session, "alpha_group", 
                          choices = c("None", meta_columns), 
                          selected = "BodySite")
        updateSelectInput(session, "beta_group", 
                          choices = c("None", meta_columns), 
                          selected = "BodySite")
        updateSelectInput(session, "shape_by", 
                          choices = c("None", meta_columns), 
                          selected = "None")
        updateSelectInput(session, "meta_x",
                          choices = c("None", meta_columns),
                          selected = meta_columns[1])
        updateSelectInput(session, "meta_y", 
                          choices = c("None", meta_columns),
                          selected = meta_columns[2])
        updateSelectInput(session, "meta_color", 
                          choices = c("None", meta_columns), 
                          selected = "BodySite")
        showNotification("Example data loaded successfully!", type = "message")
      } else {
        # Code to load real data would go here 
        # For now, placeholder with an error message 
        showNotification("Real data loading not implemented yet. Please use the example data option.", 
                         type = "error", duration = 10)
      }
    })
    
    # Generate data summary 
    output$data_summary <- renderPrint({
      ps <- phyloseq_obj()
      if (is.null(ps)) {
        return("No data loaded. Please upload data files or load example data.")
      }
      
      cat("Microbiome Dataset Summary:\n")
      cat("Number of samples:", nsamples(ps), "\n")
      cat("Number of taxa:", ntaxa(ps), "\n")
      cat("Total reads:", sum(sample_sums(ps)), "\n")
      cat("Mean reads per sample:", mean(sample_sums(ps)), "\n")
      cat("Median reads per sample:", median(sample_sums(ps)), "\n")
      cat("\nTop phyla by abundance:\n")
      top_phyla <- tax_glom(ps, taxrank = "Phylum") %>%
        transform_sample_counts(function(x) x / sum(x) * 100) %>%
        psmelt() %>%
        group_by(Phylum) %>%
        summarize(MeanRelativeAbundance = mean(Abundance)) %>%
        arrange(desc(MeanRelativeAbundance)) %>%
        head(5)
      print(top_phyla)
      cat("\nSample distribution by body site:\n")
      sample_counts <- table(sample_data(ps)$BodySite)
      print(sample_counts)
    })
    
    # Preview metadata 
    output$metadata_preview <- renderDataTable({
      md <- metadata()
      if (is.null(md)) {
        return(NULL)
      }
      datatable(md, options = list(pageLength = 5, scrollX = TRUE))
    })
    
    # Taxonomy plot 
    output$taxonomy_plot <- renderPlot({
      ps <- phyloseq_obj()
      if (is.null(ps)) {
        return(NULL)
      }
      
      # Get selected taxonomic level 
      tax_level <- input$tax_level
      
      # Agglomerate at the selected taxonomic level 
      ps_glom <- tax_glom(ps, taxrank = tax_level)
      
      # Transform to relavtive abundance
      ps_rel <- transform_sample_counts(ps_glom, function(x) x / sum(x) * 100)
      
      # Melt to long format for ggplot 
      ps_melt <- psmelt(ps_rel)
      
      # Filter low abundance taxa 
      ps_melt <- ps_melt %>%
        group_by(get(tax_level)) %>%
        mutate(MeanAbundance = mean(Abundance)) %>%
        ungroup()
      
      low_abundance_taxa <- ps_melt %>%
        filter(MeanAbundance < input$min_abundance) %>%
        pull(!!sym(tax_level)) %>%
        unique()
      
      if (length(low_abundance_taxa) > 0) {
        ps_melt <- ps_melt %>%
          mutate(!!tax_level := ifelse(get(tax_level) %in% low_abundance_taxa, "Other", get(tax_level)))
      }
      
      # Create plot based on plot type
      if (input$plot_type == "bar") {
        # Group by the selected variable if not "None"
        if (input$group_var != "None") {
          group_var <- input$group_var
          
          # Reorder samples by group and abundance if requested 
          if (input$sort_by_abundance) {
            sample_order <- ps_melt %>%
              group_by(Sample, !!sym(group_var)) %>%
              summarize(TotalAbundance = sum(Abundance), .groups = "drop") %>%
              arrange(!!sym(group_var), desc(TotalAbundance)) %>%
              pull(Sample)
            ps_melt$Sample <- factor(ps_melt$Sample, levels = sample_order)
          }
          
          # Create the grouped stacked bar plot 
          p <- ggplot(ps_melt, aes(x = Sample, y = Abundance, fill = get(tax_level))) +
            geom_bar(stat = "identity") +
            labs(x = "Sample", y = "Relative Abundance (%)", fill = "tax_level") +
            theme_minimal() +
            theme(
              axis.text.x = element_text(angle = 45, hjust = 1, size = 8), 
              legend.position = if (input$show_legend) "right" else "none"
            ) +
            facet_grid(. ~ get(group_var), scales = "free_x", space = "free_x")
        } else {
          # Create the non-grouped stacked bar plot. Reorder samples by abundance if requested
          if (input$sort_by_abundance) {
            sample_order <- ps_melt %>%
              group_by(Sample) %>%
              summarize(TotalAbundance = sum(Abundance), .groups = "drop") %>%
              arrange(desc(TotalAbundance)) %>%
              pull(Sample)
            ps_melt$Sample <- factor(ps_melt$Sample, levels = sample_order)
          }
          
          p <- ggplot(ps_melt, aes(x = Sample, y = Abundance, fill = get(tax_level))) +
            geom_bar(stat = "identity") +
            labs(x = "Sample", y = "Relative Abundance (%)", fill = tax_level) +
            theme_minimal() +
            theme(
              axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
              legend.position = if(input$show_legend) "right" else "none"
            )
        }
        
        # Use a colorful palette 
        if (length(unique(ps_melt[[tax_level]])) <= 9) {
          p <- p + scale_fill_brewer(palette = "Set1")
        } else {
          p <- p + scale_fill_viridis_d()
        }
      } else if (input$plot_type == "heatmap") {
        # Create heatmap for taxonomic abundance. Summarize data for heatmap
        if (input$group_var != "None") {
          # Group samples for heatmap
          heat_data <- ps_melt %>%
            group_by(get(tax_level), !!sym(input$group_var)) %>%
            summarize(MeanAbundance = mean(Abundance), .groups = "drop") %>%
            spread(key = !!sym(input$group_var), value = MeanAbundance, fill = 0)
          
          # Convert back to long format for ggplot 
          heat_data_long <- gather(heat_data, key = "Group", value = "Abundance", -1)
          colnames(heat_data_long)[1] <- tax_level
          
          # Order taxa by overall abundance 
          taxa_order <- ps_melt %>%
            group_by(!!sym(tax_level)) %>%
            summarize(MeanAbundance = mean(Abundance), .groups = "drop") %>%
            arrange(desc(MeanAbundance)) %>%
            pull(!!sym(tax_level))
          heat_data_long[[tax_level]] <- factor(heat_data_long[[tax_level]], levels = taxa_order)
          
          # Create heatmap
          p <- ggplot(heat_data_long, aes(x = Group, y =!!sym(tax_level), fill = Abundance)) +
            geom_tile() +
            scale_fill_viridis_c(option = "plasma") +
            labs(x = input$group_var, y = tax_level, fill = "Mean\nAbundance (%)") +
            theme_minimal() +
            theme(
              axis.text.x = element_text(angle = 45, hjust = 1),
              legend.position = "right"
            )
        } else {
          # Sample-level heatmap. Order taxa and samples
          taxa_order <- ps_melt %>%
            group_by(!!sym(tax_level)) %>%
            summarize(MeanAbundance = mean(Abundance), .groups = "drop") %>%
            arrange(desc(MeanAbundance)) %>%
            head(20) %>% # Limit to top taxa for readability 
            pull(!!sym(tax_level))
          
          # Filter for top taxa 
          ps_melt_filtered <- ps_melt %>%
            filter(!!sym(tax_level) %in% taxa_order)
          
          ps_melt_filtered[[tax_level]] <- factor(ps_melt_filtered[[tax_level]], levels = taxa_order)
          
          # Create heatmap
          p <- ggplot(ps_melt_filtered, aes(x = Sample, y = !!sym(tax_level), fill = Abundance)) +
            geom_tile() +
            scale_fill_viridis_c(option = "plasma") +
            labs(x = "Sample", y = tax_level, fill = "Abundance (%)") +
            theme_minimal() +
            theme(
              axis.text.x = element_text(angle = 90, hjust = 1, size = 7),
              legend.position = "right"
            )
        }
      }
      return(p)
    })
    
    # Top taxa table
    output$taxa_table <- renderDataTable({
      ps <- phyloseq_obj()
      if (is.null(ps)) {
        return(NULL)
      }
      # Get selected taxonomic level 
      tax_level <- input$tax_level
      
      # Agglomerate and calculate relative abundance 
      ps_glom <- tax_glom(ps, taxrank = tax_level)
      ps_rel <- transform_sample_counts(ps_glom, function(x) x / sum(x) * 100)
      ps_melt <- psmelt(ps_rel)
      
      # Summarize by selected taxonomic level 
      taxa_summary <- ps_melt %>%
        group_by(!!sym(tax_level)) %>%
        summarize(
          MeanRelativeAbundance = mean(Abundance),
          MaxRelativeAbundance = max(Abundance), 
          Prevalence = sum(Abundance > 0) / length(unique(Sample)) * 100, 
          .groups = "drop"
        ) %>%
        arrange(desc(MeanRelativeAbundance))
      
      datatable(taxa_summary, 
                options = list(pageLength = 10, scrollX = TRUE),
                rownames = FALSE) %>%
        formatRound(columns = c("MeanRelativeAbundance", "MaxRelativeAbundance", "Prevalence"), digits = 2)
    })
    
  
    # Alpha diversity plot interactive
    output$alpha_plot <- renderPlotly({
      ps <-  phyloseq_obj()
      if (is.null(ps)) {
        return(NULL)
      }
      
      # Calculate alpha diversity 
      alpha_metric <- input$alpha_metric 
      if (!alpha_metric %in% c("Faith PD")) { # Standandard metrics 
        alpha_div <- estimate_richness(ps, measures = alpha_metric)
        alpha_div$SampleID <- rownames(alpha_div)
        
        # Merge with metadata 
        md <- as(sample_data(ps), "data.frame")
        alpha_data <- merge(alpha_div, md, by.x = "SampleID", by.y = "row.names")
        
        # Set up plotting 
        if (input$alpha_group != "None") {
          group_var <- input$alpha_group
          
          # Create the alpha diversity plot
          p <- ggplot(alpha_data, aes(x = !!sym(group_var), y = !!sym(alpha_metric), 
                                      color = !!sym(group_var),
                                      text = paste("Sample:", SampleID, 
                                                   "<br>", alpha_metric, ":", round(!!sym(alpha_metric), 2),
                                                   "<br>Group:", !!sym(group_var)))) +
            labs(x = group_var, y = paste(alpha_metric, "Diversity")) +
            theme_minimal() +
            theme(legend.position = "none")
          
          # Add points if requested 
          if(input$add_points) {
            p <- p + geom_jitter(width = 0.2, alpha = 0.7, size = 3)
          }
          
          # Add boxplot if requested 
          if (input$add_boxplot) {
            p <- p + geom_boxplot(alpha = 0.5, outlier.shape = NA)
          }
          
          # Add statistical test if requested 
          if (input$stat_test != "None" && length(unique(alpha_data[[group_var]])) > 1) {
            test_result <- NULL
            p_value <- NA
            
            #### 
            if (input$stat_test == "Krustal-Wallis") {
              # Run Kruskal-Wallis test
              test_formula <- as.formula(paste(alpha_metric, "~", group_var))
              test_result <- kruskal.test(test_formula, data = alpha_data)
              p_value <- test_result$p.value
            } else if (input$stat_test == "ANOVA") {
              # Run ANOVA 
              test_formula <- as.formula(paste(alpha_metric, "~", group_var))
              test_result <- aov(test_formula, data = alpha_data)
              p_value <- summary(test_result)[[1]]$"Pr(>F)"[1]
            } else if (input$stat_test %in% c("T-test", "Wilcoxon") && 
                       length(unique(alpha_data[[group_var]])) == 2) {
              # Run t-test or Wilcoxon test for two groups
              group_levels <- unique(alpha_data[[group_var]])
              group1_data <- alpha_data[alpha_data[[group_var]] == group_levels[1], alpha_metric]
              group2_data <- alpha_data[alpha_data[[group_var]] == group_levels[2], alpha_metric]
              
              if (input$stat_test == "T-test") {
                test_result <- t.test(group1_data, group2_data)
              } else {
                test_result <- wilcox.test(group1_data, group2_data)
              }
              p_value <- test_result$p.value
            }
            
            # Add p-value annotation if test was run 
            if (!is.na(p_value)) {
              p_value_formatted <- format(p_value, digits =3)
              if (p_value < 0.001) p_value_formatted <- "p < 0.001"
              
              p <- p + labs(subtitle = paste(input$stat_test, "p-value:", p_value_formatted))
              
              # Store test result for detailes output 
              output$alpha_stats <- renderPrint({
                print(test_result)
              })
            }
          }
        } else {
          # Simple histogram if no grouping 
          p <- ggplot(alpha_data, aes(x = !!sym(alpha_metric),
                                      text = paste("Value:", round(!!sym(alpha_metric), 2), 
                                                   "<br>Sample:", SampleID))) +
            geom_histogram(fill = "steelblue", color ="black", bins = 15) +
            labs(x = paste(alpha_metric, "Diversity"), y = "Count") +
            theme_minimal()
        }
      } else {
        # Placeholder for Faith's PD (would require phylogenetic tree)
        p <- ggplot() +
          annotate("text", x = 0.5, y = 0.5, 
                   label = "Faith's PD requires a phylogenetic tree.\nNot available in example data.") +
          theme_void()
      } ##############
      # Convert ggplot to plotly
      p_interactive <- ggplotly(p, tooltip = "text") %>%
        layout(hoverlabel = list(bgcolor ="white", font = list(size = 12)), 
               dragmode = "zoom") %>%
        config(displayModeBar = TRUE, scrollZoom = TRUE)
      return(p_interactive)
    })
    
    # Rarefaction curve 
    output$rarefaction_plot <- renderPlotly({
      ps <- phyloseq_obj() 
      if (is.null(ps)) {
        return(NULL)
      }
      
      # Get OTU table and ensure samples are rows 
      otu <- as(otu_table(ps), "matrix")
      if (taxa_are_rows(ps)) {
        otu <- t(otu) 
      }
      
      # Check if we have data to work with
      if (nrow(otu) == 0 || ncol(otu) == 0) {
        return(NULL)
      }
      
      # Calculate rarefaction 
      raremax <- min(rowSums(otu))
      step_size <- max(floor(raremax/20), 1) # Ensure step size is at least 1
      
      # Get rarefaction data using vegan tidy approach 
      rarecurve_data <- vegan::rarecurve(otu, step = step_size, tidy = TRUE)
      
      # Rename columns for clarity 
      colnames(rarecurve_data) <- c("SampleName", "Reads", "Species")
      
      # Add grouping information if requested 
      if (input$alpha_group != "None") {
        # Get metadata 
        md <- as(sample_data(ps), "data.frame")
        md$SampleName <- rownames(md)
        
        # Merge with metadata 
        rarecurve_data$SampleName <- as.character(rarecurve_data$SampleName)
        rarecurve_data <- merge(rarecurve_data, md, by = "SampleName")
        
        # Check if grouping variable exists in the metadata
        if (!input$alpha_group %in% colnames(rarecurve_data)) {
          # Create plot with grouping
          p <- ggplot(rarecurve_data, aes(x = Reads, y = Species, group = SampleName,
                                          color = .data[[input$alpha_group]], 
                                          text = paste("Sample:", SampleName, 
                                                       "<br>Reads:", Reads, 
                                                       "<br>ASVs/OTUs:", Species))) +
            geom_line(linewidth = 0.5, alpha = 0.7) +
            labs(x = "Sequencing Depth", y = "Observed ASVs/OTUs",
                 title = "Rarefaction Curve") +
            theme_minimal() +
            theme(legend.position = "right")
        } else {
          # Fallback without grouping 
          p <- ggplot(rarecurve_data, aes(x = Reads, y = Species, group = SampleName, 
                                          color = SampleName, 
                                          text = paste("Sample:", SampleName, 
                                                       "<br>Reads:", Reads, 
                                                       "<br>ASVs/OTUs:", Species))) +
            geom_line(linewidth = 0.5, alpha = 0.7) +
            labs(x = "Sequencing Depth", y = "Observed ASVs/OTUs",
                 title = "Rarefaction Curve") +
            theme_minimal()
        }
      } else {
        # Create plot without grouping 
        p <- ggplot(rarecurve_data, aes(x = Reads, y = Species, group = SampleName,
                                        color = SampleName, 
                                        text = paste("Sample:", SampleName,
                                                     "<br>Reads:", Reads, 
                                                     "<br>ASVs/OTUs:", Species))) +
          geom_line(linewidth = 0.5, alpha = 0.7) +
          labs(x = "Sequencing Depth", y = "Observed ASVs/OTUs",
               title = "Rarefaction Curve") +
          theme_minimal() +
          scale_colour_viridis_d(option = "plasma")
      }
      
      # Convert to plotly for interactivity 
      p_interactive <- ggplotly(p, tooltip = "text")
      
      # Improve layout 
      p_interactive <- p_interactive %>%
        layout(hoverlabel = list(bgcolor = "white"), 
               legend = list(title = list(text = input$alpha_group))) %>%
        config(displayModeBar = TRUE, scrollZoom = TRUE)
      
      return(p_interactive)
    })
    
    # Beta diversity plot
    output$beta_plot <- renderPlotly({
      ps <- phyloseq_obj()
      if (is.null(ps)) {
        return(NULL)
      }
      # Select distance metric
      dist_method <- switch(input$beta_metric,
                            "Bray-Curtis" = "bray",
                            "Jaccard" = "jaccard",
                            "UniFrac" = "unifrac",
                            "Weighted UniFrac" = "wunifrac",
                            "bray") # Default to Bray-Curtis
      # Skip UniFrac methods if no tree is available
      if (dist_method %in% c("unifrac", "wunifrac") && is.null(phy_tree(ps, errorIfNull = FALSE))) {
        p <- ggplot() +
          annotate("text", x = 0.5, y = 0.5,
                   label = "UniFrac requires a phylogenetic tree.\nNot available in example data.") +
          theme_void()
        return(ggplotly(p))
      }
      # Calculate distance matrix
      dist_matrix <- phyloseq::distance(ps, method = dist_method)
      # Perform ordination
      ord_method <- switch(input$ordination,
                           "PCoA" = "PCoA",
                           "NMDS" = "NMDS",
                           "t-SNE" = "t-SNE",
                           "UMAP" = "UMAP",
                           "PCoA") # Default to PCoA
      # Extract metadata
      metadata_df <- as(sample_data(ps), "data.frame")
      
      # Handle different ordination methods 
      if (ord_method == "PCoA") {
        ord <- ordinate(ps, method = ord_method, distance = dist_matrix)
        ord_data <- plot_ordination(ps, ord, justDF = TRUE)
        axis_labels <- c(paste0("PCo1 [", round(ord$values$Relative_eig[1] * 100, 1), "%"),
                         paste0("PCo2 [", round(ord$values$Relative_eig[2] * 100, 1), "%"))
      } else if (ord_method == "NMDS") {
        ord <- ordinate(ps, method = ord_method, distance = dist_matrix)
        ord_data <- plot_ordination(ps, ord, justDF = TRUE)
        axis_labels <- c("NMDS1", "NMDS2")
      } else if (ord_method == "t-SNE") {
        tsne_result <- Rtsne::Rtsne(dist_matrix, is_distance = TRUE, perplexity = min(30, nrow(metadata_df) - 1))
        ord_data <- data.frame(
          Axis.1 = tsne_result$Y[, 1], 
          Axis.2 = tsne_result$Y[, 2], 
          row.names = rownames(metadata_df)
        )
        # Add sample data 
        ord_data <- cbind(ord_data, metadata_df)
        axis_labels <- c("t-SNE1", "t-SNE2")
      } else if (ord_method == "UMAP") {
        umap_config <- umap::umap.defaults
        umap_config$n_neighbors <- min(15, nrow(metadata_df) - 1)
        umap_result <- umap::umap(dist_matrix, config = umap_config)
        ord_data <- data.frame(
          Axis.1 <- umap_result$layout[, 1], 
          Axis.2 <- umap_result$layout[, 2], 
          row.names = rownames(metadata_df)
        )
        ord_data <- cbind(ord_data, metadata_df)
        axis_labels <- c("UMAP1", "UMAP2")
      }
      
      # Check if display boxplot or ordination plot
      if (input$beta_group != "None" && input$plot_type == "boxplot") {
        # Create distance to centroid data for boxplots 
        disp <- vegan::betadisper(dist_matrix, group = metadata_df[[input$beta_group]])
        
        # Extract distance to centroid data 
        boxplot_data <- data.frame(
          Group = metadata_df[[input$beta_group]], 
          Distance = dist$distances, 
          Sample = rownames(metadata_df)
        )
        
        # Create boxplot 
        p <- ggplot(boxplot_data, aes(x = Group, y = Distance, fill = Group)) +
          geom_boxplot(alpha = 0.7, outlier.shape = NA) +
          geom_jitter(aes(text = Sample), width = 0.2, alpha = 0.7) +
          theme_minimal() +
          labs(
            x = input$beta_group, 
            y = paste("Distance to centroid -", input$beta_metric), 
            title = paste("Beta dispersion using", input$beta_metric, "distances")
          )
        
        # Run PERMANOVA if requested 
        if (input$permanova == "Yes") {
          # Run PERMANOVA test
          perm_formula <- as.formula(paste("dist_matrix ~", input$beta_group))
          permanova_result <- vegan::adonis2(perm_formula, data = metadata_df)
          
          # Display results 
          output$permanova_results <- renderPrint({
            cat("PERMANOVA Results for", input$beta_metric, "distance by", input$beta_group, "\n\n")
            print(permanova_result)
            # Calculate and display variance explained 
            r2 <- permanova_result$R2[1]
            cat("\nVariance explained:", round(r2 * 100, 2), "%\n")
            
            # Add BETADISPER results (test of homogeneity of dispersions)
            cat("\nBETADISPER Results (Homogeneity of Dispersions Test):\n")
            anova_result <- anova(disp)
            print(anova_result)
          })
          
          # Add R² and p-value to plot title 
          r2 <- permanova_result$R2[1]
          p_val <- permanova_result$`Pr(>F)`[1]
          p_val_text <- ifelse(p_val < 0.001, "p < 0.001", paste("p =", round(p_val, 3)))
          
          # Add betadisper p-value 
          betadisp_pval <- anova(disp)$`Pr(>F)`[1]
          betadisp_text <- ifelse(betadisp_pval < 0.001, "p < 0.001", paste("p =", round(betadisp_pval, 3)))
          
          p <- p + labs(subtitle = paste0("PERMANOVA: R² =", round(r2, 3), ", ", p_val_text, 
                                          "\nBETADISPER: ", betadisp_text))
        }
        # Make the boxplot interactive 
        int <- ggplotly(p, tooltip = c("x", "y", "text"))
        return(int)
      } else {
        # Regular ordination plot. Set up plot aesthetics 
        aes_list <- list(x = sym("Axis.1"), y = sym("Axis.2"))
        
        # Add color aesthetic if grouping is provided 
        if (input$beta_group != "None") {
          aes_list$color <- sym(input$beta_group)
        }
        # Add shape aesthetic if requested 
        if (input$shape_by != "None" && input$shape_by != input$beta_group) {
          aes_list$shape <- sym(input$shape_by)
        }
        # Create tooltip text
        if (input$beta_group != "None" && input$shape_by != "None" && input$shape_by != input$beta_group) {
          ord_data$tooltip <- paste("Sample:", rownames(ord_data), 
                                    "\n", input$beta_group, ":", ord_data[[input$beta_group]],
                                    "\n", input$shape_by, ":", ord_data[[input$shape_by]])
        } else if (input$beta_group != "None") {
          ord_data$tooltip <- paste("Sample:", rownames(ord_data), 
                                    "\n", input$beta_group, ":", ord_data[[input$beta_group]])
        } else {
          ord_data$tooltip <- paste("Sample:", rownames(ord_data))
        }
        
        aes_list$text <- sym("tooltip")
        
        # Create the base plot
        p <- ggplot(ord_data, do.call(aes, aes_list)) +
          geom_point(size = 3, alpha = 0.7) +
          theme_minimal() +
          labs(
            x = axis_labels[1],
            y = axis_labels[2],
            title = paste(input$beta_metric, "Distance -", ord_method, "Ordination")
          )
        
        # Add ellipses if requested and grouping is provided
        if (input$add_ellipse && input$beta_group != "None" && input$plot_type == "ellipse") {
          p <- p + stat_ellipse(aes(color = !!sym(input$beta_group)), type = "norm", level = 0.95)
        }
        
        # Add sample labels if requested
        if (input$add_labels) {
          p <- p + geom_text(aes(label = rownames(ord_data)), size = 3, vjust = -1, check_overlap = TRUE)
        }
        
        # Run PERMANOVA if requested
        if (input$beta_group != "None" && input$permanova == "Yes") {
          # Run PERMANOVA test
          perm_formula <- as.formula(paste("dist_matrix ~", input$beta_group))
          permanova_result <- vegan::adonis2(perm_formula, data = metadata_df)
          
          # Display results
          output$permanova_results <- renderPrint({
            cat("PERMANOVA Results for", input$beta_metric, "distance by", input$beta_group, "\n\n")
            print(permanova_result)
            # Calculate and display variance explained
            r2 <- permanova_result$R2[1]
            cat("\nVariance explained:", round(r2 * 100, 2), "%\n")
          })
          
          # Add R² and p-value to plot title
          r2 <- permanova_result$R2[1]
          p_val <- permanova_result$`Pr(>F)`[1]
          p_val_text <- ifelse(p_val < 0.001, "p < 0.001", paste("p =", round(p_val, 3)))
          p <- p + labs(subtitle = paste0("PERMANOVA: R² = ", round(r2, 3), ", ", p_val_text))
        }
        # Convert to plotly
        plt <- ggplotly(p, tooltip = "text")
        
        # Customize plotly layout
        plt <- plt %>% layout(
          hoverlabel = list(bgcolor = "white", font = list(size = 12)),
          legend = list(title = list(text = input$beta_group))
        )
        return(plt)
      }
    })
    
    # Download handler for beta diversity plot
    output$download_beta_plot <- downloadHandler(
      filename = function() {
        paste("beta_diversity_", input$beta_metric, "_", input$ordination, "_", Sys.Date(), ".png", sep = "")
      },
      content = function(file) {
        # Get the plotly object
        p <- plotly_build(output$beta_plot)
        
        # Save as a PNG file
        png(file, width = 1000, height = 800, res = 100)
        print(ggplotly(p))
        dev.off()
      }
    )
  }
  
  # Metadata plot
  output$metadata_plot <- renderPlot({
    ps <- phyloseq_obj()
    if (is.null(ps)) {
      return(NULL)
    }
    
    meta_df <- as(sample_data(ps), "data.frame")
    
    # Extract plot type
    plot_type <- input$meta_plot_type
    x_var <- input$meta_x
    y_var <- input$meta_y
    color_var <- input$meta_color
    
    # Skip plotting if X or Y is not selected and required
    if (x_var == "None" || (plot_type %in% c("scatter", "box") && y_var == "None")) {
      return(NULL)
    }
    
    p <- NULL
    
    # Scatter Plot
    if (plot_type == "scatter") {
      p <- ggplot(meta_df, aes_string(x = x_var, y = y_var)) +
        geom_point(aes_string(color = if (color_var != "None") color_var else NULL), size = 3, alpha = 0.7) +
        theme_minimal() +
        labs(x = x_var, y = y_var, color = if (color_var != "None") color_var else NULL)
      
      # Box Plot
    } else if (plot_type == "box") {
      p <- ggplot(meta_df, aes_string(x = x_var, y = y_var)) +
        geom_boxplot(aes_string(fill = if (color_var != "None") color_var else NULL), alpha = 0.7) +
        theme_minimal() +
        labs(x = x_var, y = y_var, fill = if (color_var != "None") color_var else NULL)
      
      # Bar Plot
    } else if (plot_type == "bar") {
      p <- ggplot(meta_df, aes_string(x = x_var)) +
        geom_bar(aes_string(fill = if (color_var != "None") color_var else NULL)) +
        theme_minimal() +
        labs(x = x_var, y = "Count", fill = if (color_var != "None") color_var else NULL)
    }
    
    return(p)
  })
  
  output$correlation_plot <- renderPlot({
    md <- metadata()
    if (is.null(md)) return(NULL)
    
    # Select only numeric columns for correlation
    numeric_md <- md %>%
      select(where(is.numeric))
    
    # Check if at least 2 numeric variables exist
    if (ncol(numeric_md) < 2) {
      plot.new()
      text(0.5, 0.5, "Not enough numeric metadata for correlation analysis.", cex = 1.2)
      return()
    }
    
    # Compute correlation matrix
    cor_mat <- cor(numeric_md, use = "pairwise.complete.obs", method = "spearman")
    
    # Melt for ggplot
    cor_df <- reshape2::melt(cor_mat)
    colnames(cor_df) <- c("Var1", "Var2", "Correlation")
    
    # Plot heatmap
    ggplot(cor_df, aes(x = Var1, y = Var2, fill = Correlation)) +
      geom_tile(color = "white") +
      scale_fill_gradient2(low = "blue", high = "red", mid = "white", 
                           midpoint = 0, limit = c(-1,1), space = "Lab", 
                           name = "Spearman\nCorrelation") +
      theme_minimal() +
      theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) +
      coord_fixed()
  })
  
  
  
  
  
  
  
  
  #   function(input, output) {
  #   output$taxaBarPlot <- renderPlot({
  #     mock_data <- data.frame(
  #       Sample = rep(paste0("Sample", 1:5), each = 4),
  #       Genus = rep(c("Bacteroides", "Lactobacillus", "Firmicutes", "Proteobacteria"), 5),
  #       Abundance = runif(20, 5, 50)
  #     )
  #     ggplot(mock_data, aes(x = Sample, y = Abundance, fill = Genus)) +
  #       geom_bar(stat = "identity", position = "stack") +
  #       labs(title = "Mock Genus-Level Abundance", y = "% Abundance") +
  #       theme_minimal(base_size = 15)
  #   })
  # 
  #   output$diversityPlot <- renderPlot({
  #     data <- data.frame(
  #       Group = rep(c("Control", "Treatment"), each = 30),
  #       Shannon = c(rnorm(30, 3.5, 0.5), rnorm(30, 2.8, 0.4))
  #     )
  #     ggplot(data[data$Group == input$diversityGroup, ], aes(x = Group, y = Shannon, fill = Group)) +
  #       geom_violin() +
  #       geom_boxplot(width = 0.1, fill = "white", outlier.shape = NA) +
  #       theme_minimal(base_size = 15) +
  #       scale_fill_brewer(palette = "Pastel2") +
  #       labs(title = paste("Shannon Diversity -", input$diversityGroup))
  #   })
  # 
  #   output$functionPlot <- renderPlot({
  #     pathways <- data.frame(
  #       Pathway = c("Metabolism", "Immune", "Signal", "Transport"),
  #       Score = runif(4, 0, 100)
  #     )
  #     ggplot(pathways, aes(x = reorder(Pathway, Score), y = Score, fill = Pathway)) +
  #       geom_col() +
  #       coord_flip() +
  #       theme_minimal(base_size = 15) +
  #       scale_fill_brewer(palette = "Set3") +
  #       labs(title = "Mock Pathway Activity", x = "", y = "Score")
  #   })
  # 
  #   output$metaBoxPlot <- renderPlot({
  #     data <- data.frame(
  #       Sex = rep(c("Male", "Female"), each = 40),
  #       Value = c(rnorm(40, 2.5, 0.3), rnorm(40, 2.8, 0.4))
  #     )
  #     ggplot(data, aes(x = Sex, y = Value, fill = Sex)) +
  #       geom_boxplot() +
  #       theme_minimal(base_size = 15) +
  #       scale_fill_manual(values = c("steelblue", "salmon")) +
  #       labs(title = "Mock Metadata Grouping by Sex", y = "Example Value")
  #   })
  # 
  #   output$downloadReport <- downloadHandler(
  #     filename = function() {
  #       "AGB_mock_report.txt"
  #     },
  #     content = function(file) {
  #       writeLines("This is a placeholder report.\nData and insights will appear here once available.", file)
  #     }
  #   )
  # }
  
  # Run app
  shinyApp(ui, server)
  
  
  
  
  
  
}



  
  
  
  
  
  
#   function(input, output) {
#   output$taxaBarPlot <- renderPlot({
#     mock_data <- data.frame(
#       Sample = rep(paste0("Sample", 1:5), each = 4),
#       Genus = rep(c("Bacteroides", "Lactobacillus", "Firmicutes", "Proteobacteria"), 5),
#       Abundance = runif(20, 5, 50)
#     )
#     ggplot(mock_data, aes(x = Sample, y = Abundance, fill = Genus)) +
#       geom_bar(stat = "identity", position = "stack") +
#       labs(title = "Mock Genus-Level Abundance", y = "% Abundance") +
#       theme_minimal(base_size = 15)
#   })
# 
#   output$diversityPlot <- renderPlot({
#     data <- data.frame(
#       Group = rep(c("Control", "Treatment"), each = 30),
#       Shannon = c(rnorm(30, 3.5, 0.5), rnorm(30, 2.8, 0.4))
#     )
#     ggplot(data[data$Group == input$diversityGroup, ], aes(x = Group, y = Shannon, fill = Group)) +
#       geom_violin() +
#       geom_boxplot(width = 0.1, fill = "white", outlier.shape = NA) +
#       theme_minimal(base_size = 15) +
#       scale_fill_brewer(palette = "Pastel2") +
#       labs(title = paste("Shannon Diversity -", input$diversityGroup))
#   })
# 
#   output$functionPlot <- renderPlot({
#     pathways <- data.frame(
#       Pathway = c("Metabolism", "Immune", "Signal", "Transport"),
#       Score = runif(4, 0, 100)
#     )
#     ggplot(pathways, aes(x = reorder(Pathway, Score), y = Score, fill = Pathway)) +
#       geom_col() +
#       coord_flip() +
#       theme_minimal(base_size = 15) +
#       scale_fill_brewer(palette = "Set3") +
#       labs(title = "Mock Pathway Activity", x = "", y = "Score")
#   })
# 
#   output$metaBoxPlot <- renderPlot({
#     data <- data.frame(
#       Sex = rep(c("Male", "Female"), each = 40),
#       Value = c(rnorm(40, 2.5, 0.3), rnorm(40, 2.8, 0.4))
#     )
#     ggplot(data, aes(x = Sex, y = Value, fill = Sex)) +
#       geom_boxplot() +
#       theme_minimal(base_size = 15) +
#       scale_fill_manual(values = c("steelblue", "salmon")) +
#       labs(title = "Mock Metadata Grouping by Sex", y = "Example Value")
#   })
# 
#   output$downloadReport <- downloadHandler(
#     filename = function() {
#       "AGB_mock_report.txt"
#     },
#     content = function(file) {
#       writeLines("This is a placeholder report.\nData and insights will appear here once available.", file)
#     }
#   )
# }

# Run app
shinyApp(ui, server)





