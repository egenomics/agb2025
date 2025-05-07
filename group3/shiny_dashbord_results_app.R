
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
library(RColorBrewer)
# library(leaflet)

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
      menuItem(" Data Upload & Overview", tabName = "data_overview", icon = icon("upload")),
      menuItem(" Taxonomic Composition", tabName = "taxonomy", icon = icon("bacteria")),
      menuItem(" Diversity Analysis", tabName = "diversity", icon = icon("chart-line")),
      menuItem(" Clinical Factors", tabName = "clinical", icon = icon("notes-medical")),
      menuItem(" Lifestyle Impact", tabName = "lifestyle", icon = icon("heartbeat")),
      menuItem(" Medical Interventions", tabName = "interventions", icon = icon("prescription-bottle-medical")), 
      menuItem(" Multi-factor Analysis", tabName = "multifactor", icon = icon("sitemap")), 
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
    
    # TAB 1: Data Upload & Overview
    tabItem(tabName = "data_overview",
            fluidRow(
              box(
                title = "Upload Microbiome Data", 
                status = "primary", 
                solidHeader = TRUE, 
                width = 8, 
                fluidRow(
                  column(
                    width = 12,
                    fileInput("sample_metadata", "Upload Sample Metadata (CSV or TSV format):",
                              accept = c(".csv", ".tsv", ".txt")),
                    fileInput("run_metadata", "Upload Run Metadata (CSV or TSV format):",
                              accept = c(".csv", ".tsv", ".txt")),
                    fileInput("taxonomy_file", "Upload Taxonomy Data (TSV format):",
                              accept = c(".tsv", ".txt")),
                    actionButton("load_data", "Load Data", 
                                 icon = icon("upload"),
                                 class = "btn-success"),
                    actionButton("generate_example", "Generate Example Data", 
                                 icon = icon("database"),
                                 class = "btn-info")
                  )
                )
              ),
              column(
                width = 4, 
                valueBoxOutput("total_samples_box", width = 12), 
                valueBoxOutput("total_species_box", width = 12), 
                valueBoxOutput("total_runs_box", width = 12)
              )
            ), 
            
            # Metadata Previews
            fluidRow(
              box(
                title = "Metadata Tables", 
                status = "primary",
                solidHeader = TRUE,
                width = 12, 
                tabBox(
                  width = 12, 
                  tabPanel("Sample Metadata",
                           dataTableOutput("sample_metadata_preview")), 
                  tabPanel("Run Metadata",
                           dataTableOutput("run_metadata_preview")), 
                  tabPanel("Taxonomy Data Preview", 
                           dataTableOutput("taxonomy_data_preview"))
                )
              ) 
            ), 
            
            # Data Overview
            fluidRow(
              box(
                title = "Data Summary",
                status = "primary", 
                solidHeader = TRUE, 
                width = 12,
                uiOutput("data_summary_dynamic")
              )
            )
    ),
    #-------------------------------------------------------------
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
  data_store <- reactiveValues(
    sample_metadata = NULL,
    run_metadata = NULL, 
    taxonomy_data = NULL
  )
  
  # Generate example data when button is clicked 
  observeEvent(input$generate_example, {
    # Generate sample metadata 
    sample_metadata <- data.frame(
      Run_ID = paste0("RUN", sprintf("%03d", 1:50)),
      Sample_ID = paste0("SAMPLE", sprintf("%03d", 1:50)),
      Institution = sample(c("Hospital A", "Hospital B", "Hospital C"), 50, replace = TRUE),
      Department = sample(c("Gastroenterology", "Immunology", "Microbiology"), 50, replace = TRUE),
      Collection_Date = as.character(sample(seq(as.Date('2023-01-01'), as.Date('2025-04-30'), by="day"), 50)),
      Collection_Storage_Temperature = sample(c("-80", "-20", "4", "Room Temperature"), 50, replace = TRUE),
      Analyst_Processor_Name = sample(c("John Doe", "Jane Smith", "Alex Johnson"), 50, replace = TRUE),
      Gender = sample(c("Male", "Female", "Other"), 50, replace = TRUE),
      Age = sample(18:80, 50, replace = TRUE),
      Geographical_origin = sample(c("North America", "Europe", "Asia", "Africa", "South America"), 50, replace = TRUE),
      Ongoing_conditions = sample(c("None", "Diabetes", "Hypertension", "IBS", "Crohn's Disease"), 50, replace = TRUE),
      Appendix_removed = sample(c("Yes", "No"), 50, replace = TRUE),
      Neurological_Disorders = sample(c("None", "Alzheimer's", "Parkinson's", "Multiple Sclerosis"), 50, replace = TRUE),
      Allergies = sample(c("None", "Pollen", "Food", "Medicine", "Multiple"), 50, replace = TRUE),
      Bowel_movement_frequency = sample(c("Daily", "2-3 times a week", "4-6 times a week", "Multiple times daily"), 50, replace = TRUE),
      Bowel_movement_quality = sample(c("Normal", "Constipated", "Loose", "Variable"), 50, replace = TRUE),
      Antibiotic_intake = sample(c("None in past year", "Within past month", "Within past 6 months", "Currently taking"), 50, replace = TRUE),
      Medications = sample(c("None", "Antidepressants", "Blood pressure medication", "Multiple medications"), 50, replace = TRUE),
      Cancer = sample(c("No", "Yes, currently", "Yes, in remission"), 50, replace = TRUE),
      Cancer_treatment = sample(c("None", "Chemotherapy", "Radiation", "Surgery", "Combination"), 50, replace = TRUE),
      BMI = round(rnorm(50, mean = 25, sd = 4), 1),
      Exercise_frequency = sample(c("Never", "1-2 times a week", "3-4 times a week", "5+ times a week"), 50, replace = TRUE),
      Exercise_intensity = sample(c("Low", "Moderate", "High", "Variable"), 50, replace = TRUE),
      Smoking_status = sample(c("Never smoked", "Former smoker", "Current smoker"), 50, replace = TRUE),
      Cigarettes_per_day = sample(c(NA, "1-5", "6-10", "11-20", "20+"), 50, replace = TRUE),
      Stopped_smoking = sample(c(NA, "Less than 1 year ago", "1-5 years ago", "5+ years ago"), 50, replace = TRUE),
      Alcohol_consumption = sample(c("Never", "Occasionally", "Weekly", "Daily"), 50, replace = TRUE),
      Alcohol_frequency = sample(c(NA, "1-2 drinks", "3-4 drinks", "5+ drinks"), 50, replace = TRUE),
      Drinks_per_day = sample(c(NA, "1", "2", "3", "4+"), 50, replace = TRUE)
    )
    
    # Generate run metadata 
    run_metadata <- data.frame(
      Timestamp = as.character(sample(seq(as.Date('2023-01-01'), as.Date('2025-04-30'), by="day"), 30)),
      Run_ID = paste0("RUN", sprintf("%03d", 1:30)),
      Sequencing_Date = as.character(sample(seq(as.Date('2023-01-15'), as.Date('2025-05-01'), by="day"), 30)),
      Sequencing_Platform = sample(c("Illumina", "Ion Torrent", "Oxford Nanopore", "PacBio"), 30, replace = TRUE),
      Sequencing_Type = sample(c("16S rRNA", "Shotgun", "ITS", "Amplicon"), 30, replace = TRUE),
      Expected_read_length = sample(c("150bp", "250bp", "300bp", "400bp"), 30, replace = TRUE),
      Sequencing_depth_target = paste0(sample(c("10", "20", "30", "50", "100"), 30, replace = TRUE), "M reads"),
      Library_preparation_kit = sample(c("Nextera XT", "TruSeq Nano", "NEBNext", "Swift Biosciences"), 30, replace = TRUE),
      Technician_name = sample(c("Maria Garcia", "Thomas Johnson", "Sophia Lee", "Michael Brown"), 30, replace = TRUE)
    )
    
    # Generate taonomy data 
    species_names <- c(
      "Bacteroides fragilis", "Escherichia coli", "Lactobacillus acidophilus", "Bifidobacterium longum", 
      "Prevotella copri", "Faecalibacterium prausnitzii", "Akkermansia muciniphila", "Ruminococcus bromii",
      "Clostridium difficile", "Enterococcus faecalis", "Blautia obeum", "Roseburia intestinalis",
      "Streptococcus thermophilus", "Eubacterium rectale", "Methanobrevibacter smithii"
    )
    
    taxonomy_data <- data.frame(
      Sample_ID = rep(paste0("SAMPLE", sprintf("%03d", 1:50)), each = length(species_names)),
      Species = rep(species_names, 50),
      Abundance = round(runif(50 * length(species_names), min = 0, max = 10000))
    )
    
    # Store the data in reactive values 
    data_store$sample_metadata <- sample_metadata
    data_store$run_metadata <- run_metadata
    data_store$taxonomy_data <- taxonomy_data 
    
    # Show notification 
    showNotification("Example data generated successfully!", type = "message")
  })
  
  # Load data when button is clicked 
  observeEvent(input$load_data, {
    # Check if all required files are uploaded 
    req(input$sample_metadata, input$run_metadata, input$taxonomy_file)
    
    # Read sample metadata 
    sample_metadata <- read.csv(input$sample_metadata$datapath)
    
    # Read run metadata 
    run_metadata <- read.csv(input$run_metadata$datapath)
    
    # Read taxonomy data 
    taxonomy_data <- read.delim(input$taxonomy_file$datapath, sep = "\t")
    
    # Store the data in reactive values 
    data_store$sample_metadata <- sample_metadata
    data_store$run_metadata <- run_metadata
    data_store$taxonomy_data <- taxonomy_data
    
    # Show notification 
    showNotification("Data loaded successfully!", type = "message")
  })
  
  # Value boxes
  output$total_samples_box <- renderValueBox({
    req(data_store$sample_metadata)
    valueBox(
      nrow(data_store$sample_metadata), 
      "Total Samples", 
      icon = icon("users"), 
      color = "blue"
    )
  })
  
  output$total_species_box <- renderValueBox({
    req(data_store$taxonomy_data)
    valueBox(
      length(unique(data_store$taxonomy_data$Species)), 
      "Unique Species", 
      icon = icon("bacteria"), 
      color = "green"
    )
  })
  
  output$total_runs_box <- renderValueBox({
    req(data_store$run_metadata)
    valueBox(
      nrow(data_store$run_metadata), 
      "Total Runs", 
      icon = icon("server"),
      color = "purple"
    )
  })
  
  # Dynamic summary content 
  output$data_summary_dynamic <- renderUI({
    req(data_store$sample_metadata, data_store$run_metadata, data_store$taxonomy_data)
    
    # General Overview 
    total_samples <- nrow(data_store$sample_metadata)
    total_runs <- nrow(data_store$run_metadata)
    total_species <- length(unique(data_store$taxonomy_data$Species))
    sample_date_range <- paste(min(as.Date(data_store$sample_metadata$Collection_Date)), "to", 
                               max(as.Date(data_store$sample_metadata$Collection_Date)))
    seq_date_range <- paste(min(as.Date(data_store$sample_metadata$Sequencing_Date)), "to", 
                            max(as.Date(data_store$sample_metadata$Sequencing_Date)))
    
    # Sequencing Summary
    platform_counts <- table(data_store$run_metadata$Sequencing_Platform)
    platforms <- paste(names(platform_counts), paste0("(", platform_counts, ")"), collapse = ", ")
    seq_types <- paste(unique(data_store$run_metadata$Sequencing_Type), collapse = ", ")
    lib_kits <- paste(unique(data_store$run_metadata$Library_preparation_kit), collapse = ", ")
    read_lengths <- paste(unique(data_store$run_metadata$Expected_read_length), collapse = ", ")
    depth_targets <- paste(unique(data_store$run_metadata$Sequencing_depth_target), collapse = ", ")
    
    # Demographics 
    gender_counts <- table(data_store$sample_metadata$Gender)
    gender_dist <- paste(names(gender_counts), paste0("(", gender_counts, ")"), collapse = ", ")
    age_range <- paste(min(data_store$sample_metadata$Age), "-", max(data_store$sample_metadata$Age), 
                       "(mean:", round(mean(data_store$sample_metadata$Age), 1), ")")
    bmi_stats <- paste("mean:", round(mean(data_store$sample_metadata$BMI), 1),
                       "± SD:", round(sd(data_store$sample_metadata$BMI), 1))
    
    geo_counts <- sort(table(data_store$sample_metadata$Geographical_origin), decreasing = TRUE)
    geo_origins <- paste(names(geo_counts), paste0("(", geo_counts, ")"), collapse = ", ")
    
    condition_counts <- sort(table(data_store$sample_metadata$Ongoing_conditions), decreasing = TRUE)
    top_conditions <- paste(names(condition_counts)[1:min(3, length(condition_counts))],
                            paste0("(", condition_counts[1:min(3, length(condition_counts))], ")"), 
                            collapse = ", ")
    
    # Health & Lifestyle 
    neuro_counts <- table(data_store$sample_metadata$Neurological_Disorders)
    neuro_disorders <- paste(names(neuro_counts), paste0("(", neuro_counts, ")"), collapse = ", ")
    
    allergy_counts <- table(data_store$sample_metadata$Allergies)
    allergies <- paste(names(allergy_counts), paste0("(", allergy_counts, ")"), collapse = ", ")
    
    smoking_counts <- table(data_store$sample_metadata$Smoking)
    smoking <- paste(names(smoking_counts), paste0("(", smoking_counts, ")"), collapse = ", ")
    
    alcohol_counts <- table(data_store$sample_metadata$Alcohol_Consumption)
    alcohol <- paste(names(alcohol_counts), paste0("(", alcohol_counts, ")"), collapse = ", ")
    
    exercise_counts <- table(data_store$sample_metadata$Exercise_frequency)
    exercise <- paste(names(exercise_counts), paste0("(", exercise_counts, ")"), collapse = ", ")
    
    # Taxonomy Overview
    # Calculate mean abundance for each species across all samples
    species_abundance <- aggregate(Abundance ~ Species, data = data_store$taxonomy_data, mean)
    species_abundance <- species_abundance[order(species_abundance$Abundance, decreasing = TRUE), ]
    top_species <- paste(species_abundance$Species[1:min(5, nrow(species_abundance))], collapse = ", ")
    
    # Count samples where each species is present
    species_presence <- table(data_store$taxonomy_data$Species[data_store$taxonomy_data$Abundance > 0])
    common_species <- sum(species_presence > (total_samples * 0.5))
    
    # Calculate species richness per sample
    species_richness <- aggregate(Species ~ Sample_ID, 
                                  data = data_store$taxonomy_data[data_store$taxonomy_data$Abundance > 0, ], 
                                  function(x) length(unique(x)))
    mean_richness <- round(mean(species_richness$Species), 1)
    
    # Return an HTML structure
    HTML(paste0('
      <div style="padding: 15px;">
        <h4>📊 <strong>General Overview</strong></h4>
        <ul>
          <li><strong>Total Samples:</strong> ', total_samples, '</li>
          <li><strong>Total Runs:</strong> ', total_runs, '</li>
          <li><strong>Total Unique Species:</strong> ', total_species, '</li>
          <li><strong>Sample Collection Date Range:</strong> ', sample_date_range, '</li>
          <li><strong>Sequencing Date Range:</strong> ', seq_date_range, '</li>
        </ul>
        
        <h4>🧬 <strong>Sequencing Summary</strong></h4>
        <ul>
          <li><strong>Platforms used:</strong> ', platforms, '</li>
          <li><strong>Types:</strong> ', seq_types, '</li>
          <li><strong>Library Kits:</strong> ', lib_kits, '</li>
          <li><strong>Read lengths:</strong> ', read_lengths, '</li>
          <li><strong>Depth targets:</strong> ', depth_targets, '</li>
        </ul>
        
        <h4>👥 <strong>Participant Demographics</strong></h4>
        <ul>
          <li><strong>Gender:</strong> ', gender_dist, '</li>
          <li><strong>Age range:</strong> ', age_range, '</li>
          <li><strong>BMI:</strong> ', bmi_stats, '</li>
          <li><strong>Geographical Origins:</strong> ', geo_origins, '</li>
          <li><strong>Ongoing conditions (top 3):</strong> ', top_conditions, '</li>
        </ul>
        
        <h4>🧠 <strong>Health & Lifestyle</strong></h4>
        <ul>
          <li><strong>Neurological Disorders:</strong> ', neuro_disorders, '</li>
          <li><strong>Allergies:</strong> ', allergies, '</li>
          <li><strong>Smoking:</strong> ', smoking, '</li>
          <li><strong>Alcohol Consumption:</strong> ', alcohol, '</li>
          <li><strong>Exercise Frequency:</strong> ', exercise, '</li>
        </ul>
        
        <h4>🦠 <strong>Taxonomy Overview</strong></h4>
        <ul>
          <li><strong>Most Abundant Species:</strong> ', top_species, '</li>
          <li><strong>Sparsity:</strong> species detected in >50% of samples: ', common_species, '</li>
          <li><strong>Mean Species Richness per Sample:</strong> ', mean_richness, '</li>
        </ul>
      </div>
    '))
  })
  
  # Render metadata table previews 
  output$sample_metadata_preview <- renderDataTable({
    req(data_store$sample_metadata)
    datatable(data_store$sample_metadata, 
              options = list(pageLength = 10, scrollX = TRUE))
  })
  
  output$run_metadata_preview <- renderDataTable({
    req(data_store$run_metadata)
    datatable(data_store$run_metadata, 
              options = list(pageLength = 10, scrollX = TRUE))
  })
  
  output$taxonomy_data_preview <- renderDataTable({
    req(data_store$taxonomy_data)
    datatable(data_store$taxonomy_data, 
              options = list(pageLength = 10, scrollX = TRUE))
  })
  
  
  # -----------------------------------------------
  
  
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
        
        # Get the number of unique groups for the palette 
        n_groups <- length(unique(alpha_data[[group_var]]))
        palette_colors <- get_palette(n_groups)
        
        # Create the alpha diversity plot
        p <- ggplot(alpha_data, aes(x = !!sym(group_var), y = !!sym(alpha_metric), 
                                    color = !!sym(group_var),
                                    text = paste("Sample:", SampleID, 
                                                 "<br>", alpha_metric, ":", round(!!sym(alpha_metric), 2),
                                                 "<br>Group:", !!sym(group_var)))) +
          labs(x = group_var, y = paste(alpha_metric, "Diversity")) +
          scale_color_manual(values = palette_colors) +
          scale_fill_manual(values = palette_colors) +
          theme_minimal() +
          theme(legend.position = "right")
        
        # Add points if requested 
        if(input$add_points) {
          p <- p + geom_jitter(width = 0.2, alpha = 0.7, size = 3)
        }
        
        # Add boxplot if requested 
        if (input$add_boxplot) {
          p <- p + geom_boxplot(alpha = 0.2, fill = "white", outlier.shape = NA)
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
        Distance = disp$distances, 
        Sample = rownames(metadata_df)
      )
      # Get the number of groups for the palette 
      n_groups <- length(unique(boxplot_data$Group))
      palette_colors <- get_palette(n_groups)
      
      # Create boxplot 
      p <- ggplot(boxplot_data, aes(x = Group, y = Distance, fill = Group)) +
        geom_boxplot(alpha = 0.2, fill = "white", outlier.shape = NA) +
        geom_jitter(aes(text = Sample), width = 0.2, alpha = 0.7) +
        scale_color_manual(values = palette_colors) +
        scale_fill_manual(values = palette_colors) +
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





