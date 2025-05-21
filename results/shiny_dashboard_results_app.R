
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
library(ComplexHeatmap)
library(circlize)
library(ggbeeswarm)
library(dplyr)
library(tidyr)
library(stringr)
library(viridis)
library(FSA)


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
      menuItem(" Clinical & Lifestyle Correlations", tabName = "multifactor", icon = icon("notes-medical")),
      menuItem(" Medical Interventions", tabName = "interventions", icon = icon("prescription-bottle-medical")), 
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
                      fileInput("control_sample_metadata", "Upload Control Sample Metadata (CSV or TSV format):", 
                                accept = c(".csv", ".tsv", ".txt")),
                      fileInput("control_run_metadata", "Upload Control Run Metadata (CSV or TSV format):", 
                                accept = c(".csv", ".tsv", ".txt")),
                      fileInput("control_taxonomy", "Upload Control Taxonomy Data (TSV format):", 
                                accept = c(".tsv", ".txt")),
                      actionButton("load_data", "Load Data", 
                                   icon = icon("upload"),
                                   class = "btn-success"),
                      actionButton("generate_example", "Generate Example Data", 
                                   icon = icon("database"),
                                   class = "btn-info")
                    )
                  ),
                  br(), 
                  fluidRow(
                    box(
                      title = "Sample Selection",
                      status = "warning",
                      solidHeader = TRUE, 
                      width = 12,
                      fluidRow(
                        column(
                          width = 4,
                          selectizeInput("selected_sample_id", "Select Sample to Analyze:",
                                         choices = NULL,
                                         options = list(
                                           placeholder = "Choose a sample...", 
                                           onInitialize = I('function() { this.setValue(""); }')
                                         ))
                        ),
                        column(
                          width = 8, 
                          htmlOutput("selected_sample_info"),
                          htmlOutput("comparison_info")
                        )
                      )
                    )
                  )
                ),
                column(
                  width = 4, 
                  valueBoxOutput("total_samples_box", width = 12), 
                  valueBoxOutput("total_species_box", width = 12), 
                  valueBoxOutput("total_runs_box", width = 12), 
                  valueBoxOutput("total_control_samples_box", width = 12)
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
                             fluidRow(
                               column(width = 12, 
                                      radioButtons("sample_metadata_filter", "Filter by Condition:",
                                                   choices = c("All", "Control", "Sample"),
                                                   selected = "All", 
                                                   inline = TRUE))
                             ),
                             dataTableOutput("sample_metadata_preview")), 
                    tabPanel("Run Metadata",
                             fluidRow(
                               column(width = 12, 
                                      radioButtons("run_metadata_filter", "Filter by Condition:",
                                                   choices = c("All", "Control", "Sample"),
                                                   selected = "All", 
                                                   inline = TRUE))
                             ),
                             dataTableOutput("run_metadata_preview")), 
                    tabPanel("Taxonomy Data Preview", 
                             fluidRow(
                               column(width = 12, 
                                      radioButtons("taxonomy_data_filter", "Filter by Condition:",
                                                   choices = c("All", "Control", "Sample"),
                                                   selected = "All", 
                                                   inline = TRUE))
                             ),
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
      
      # TAB 2: Taxonomy Composition Tab 
      tabItem(tabName = "taxonomy",
              fluidRow(
                box(
                  title = "Taxonomy Composition Analysis", 
                  status = "primary", 
                  solidHeader = TRUE,
                  width = 12)
              ), 
              # Analysis controls sidebar 
              fluidRow(
                column(width = 3, 
                       box(
                         title = "Analysis Controls", 
                         status = "primary", 
                         solidHeader = TRUE, 
                         width = 12, 
                         
                         # Taxonomic level selection
                         selectInput("tax_level", "Taxonomic Level:", 
                                     choices = c("Phylum", "Class", "Order", "Family", "Genus", "Species"), 
                                     selected = "Species"), 
                         
                         # Filtering options 
                         radioButtons("sample_group", "Sample Group:", 
                                      choices = c("All Samples", "Patients Only", "Controls Only"), 
                                      selected = "All Samples"),
                         
                         hr(), 
                         
                         # Metadata grouping options 
                         selectInput("metadata_group", "Group by Metadata:", 
                                     choices = c("None", "Age", "Gender", "BMI"),
                                     selected = "None"), 
                         
                         checkboxInput("compare_groups", "Compare Patients vs Controls", value = TRUE), 
                         
                         hr(), 
                         
                         # Statistical testing options 
                         selectInput("stat_test", "Statistical Test:", 
                                     choices = c("None", "PERMANOVA", "DESeq2"), 
                                     selected = "None"), 
                         ), 
                       # Sample selection box 
                       box(
                         title = "Sample Selection", 
                         status = "warning", 
                         solidHeader = TRUE, 
                         width = 12, 
                         fluidRow(
                           column(
                             width = 12, 
                             selectizeInput("selected_sample_id", "Select Sample to Analyze:", 
                                            choices = NULL, 
                                            options = list(
                                              placeholder = "Choose a sample...",
                                              onInitialize = I('function() { this.setValue(""); }')
                                            ))
                           )
                         )
                       )
                ),
                # Main content area 
                column(width = 9,
                       # Quick stat boxes
                       fluidRow(
                         valueBoxOutput("total_taxa_box", width = 3), 
                         valueBoxOutput("dominant_taxa_box", width = 3),
                         valueBoxOutput("unique_taxa_box", width = 3)
                       ),
                       # Main visualizations 
                       box(
                         title = "Taxonomic Visualizations",
                         status ="info", 
                         solidHeader = TRUE, 
                         width = 12, 
                         height = "750px", 
                         tabBox(
                           width = 12, 
                           height = "700px", 
                          
                           # Tab 1: Relative Abundance bar plot
                           tabPanel("Relative Abundance", 
                                    fluidRow(
                                      # Left side controls 
                                      column(width = 3, 
                                             wellPanel(
                                               sliderInput("abundance_threshold", "Minimum Abundance (%):", 
                                                           min = 0, max = 5, value = 0.5, step = 0.1), 
                                               hr(), 
                                               checkboxInput("sort_abundance", "Sort by Abundance", value = TRUE),
                                               selectInput("abundance_color_scheme", "Color Scheme:", 
                                                           choices = c("Default", "Viridis", "Set1", "Set2", "Paired"), 
                                                           selected = "Default")
                                             )),
                                      # Plot area 
                                      column(width = 9, 
                                             div(style = "height: 600px; overflow-y: auto;", 
                                                 plotOutput("relative_abundance_plot", height = "550px")
                                             )
                                      )
                                    ),
                                    fluidRow(
                                      column(width = 12, 
                                             align = "right", 
                                             downloadButton("download_abundance_plot", "Download Plot",
                                                            class = "btn-sm btn-info")
                                      )
                                    )
                           ), 
                           ######### LEFT SIDE CONTROLS OR NOT? ########
                           # Tab 2: Heatmap Visualization
                           tabPanel("Heatmap", 
                                    fluidRow(
                                      column(width = 12, 
                                             # Color palette options 
                                             fluidRow(
                                               column(width = 6, 
                                                      offset = 3, 
                                                      selectInput("heatmap_color", "Color Palette:", 
                                                                  choices = c("Viridis", "Magma", "Plasma", "Inferno", "Blues", "RdBu"),
                                                                  selected = "Viridis"), 
                                                      )
                                             ),
                                             div(style = "height: 600px; overflow-y: auto;", 
                                                 plotOutput("heatmap_plot", height = "550px")
                                             )
                                      )
                                    ), 
                                    fluidRow(
                                      column(width = 12, 
                                             align = "right", 
                                             downloadButton("download_heatmap", "Download Plot", 
                                                            class = "btn-sm btn-info")
                                      )
                                    )
                           ),
                           # Tab 3: Clustering Visualization
                           tabPanel("Clustering",
                                    fluidRow(
                                      column(width = 12, 
                                             # Clustering methods and options 
                                             fluidRow(
                                               column(width = 6, 
                                                      offset = 3, 
                                                      selectInput("cluster_method", "Cluster Method:", 
                                                         choices = c("PCA", "t-SNE", "UMAP", "NMDS"), 
                                                         selected = "PCA"), 
                                                      checkboxInput("show_ellipses", "Show Group Ellipses", value = TRUE), 
                                                      checkboxInput("show_labels", "Show Sample Labels", value = FALSE)
                                                      )
                                             ), 
                                             div(style = "height: 600px; overflow-y: auto;",
                                                 plotOutput("clustering_plot", height = "550px")
                                             )
                                      )
                                    ), 
                                    fluidRow(
                                      column(width = 12, 
                                             align = "right", 
                                             downloadButton("download_clustering_plot", "Download Plot", 
                                                            class = "btn-sm btn-info")
                                      )
                                    )
                           )
                       )
                     )
                     )
              )
      ),
                  
      # TAB 3: Diversity Analysis Tab
      tabItem(tabName = "diversity",
              fluidRow(
                
                # LEFT PANEL: Only Diversity Type
                box(
                  title = "Diversity Settings",
                  status = "primary",
                  solidHeader = TRUE,
                  width = 3,
                  collapsible = TRUE,
                  
                  radioButtons("diversity_type", "Diversity Type:",
                               choices = c("Alpha Diversity" = "alpha", 
                                           "Beta Diversity" = "beta"),
                               selected = "alpha"),
                  checkboxInput("subdivide_samples", "Split samples by condition", value = FALSE),
                  
                  # Only show if checked
                  conditionalPanel(
                    condition = "input.subdivide_samples == true",
                    selectInput("condition_column", "Clinical condition to group by:",
                                choices = c("Ongoing_conditions", "Neurological_Disorders", "Allergies", "Cancer"),
                                selected = "Ongoing_conditions")
                  )
                ),
                
                # MIDDLE PANEL: Plot options
                box(
                  title = "Visualization Options",
                  status = "primary",
                  solidHeader = TRUE,
                  width = 3,
                  collapsible = TRUE,
                  
                  conditionalPanel(
                    condition = "input.diversity_type == 'alpha'",
                    
                    selectInput("alpha_metric", "Alpha Metric:",
                                choices = c("Observed OTUs", "Shannon", "Simpson"),
                                selected = "Shannon"),
                    
                    uiOutput("alpha_plot_type_ui"),
                    uiOutput("otu_warning_ui"),
                    
                    selectInput("alpha_color_palette", "Color Palette:",
                                choices = c("Set1", "Dark2", "Pastel1", "Paired", "Viridis"),
                                selected = "Set1")
                  ),
                  
                  conditionalPanel(
                    condition = "input.diversity_type == 'beta'",
                    
                    selectInput("distance_metric", "Distance Metric:",
                                choices = c("Bray-Curtis", "Euclidean", "Jaccard"),
                                selected = "Bray-Curtis"),
                    
                    selectInput("ordination_method", "Ordination Method:",
                                choices = c("PCoA", "NMDS"),
                                selected = "PCoA"),
                    
                    checkboxInput("show_group_ellipses", "Show Group Ellipses", value = TRUE),
                    checkboxInput("show_clustering_dendrogram", "Show Clustering Dendrogram", value = TRUE)
                  ),
                  
                  checkboxInput("show_diversity_legend", "Show Legend", value = TRUE),
                  hr(),
                  downloadButton("download_diversity_plot", "Download Plot", class = "btn-primary")
                ),
                
                # RIGHT PANEL: Plot Output
                # RIGHT PANEL: Plot Output
                box(
                  title = "Diversity Analysis",
                  status = "info",
                  solidHeader = TRUE,
                  width = 6,
                  
                  # ALPHA DIVERSITY TABS
                  conditionalPanel(
                    condition = "input.diversity_type == 'alpha'",
                    tabsetPanel(
                      tabPanel("Alpha Diversity Plot", plotlyOutput("diversity_main_plot", height = "600px")),
                      tabPanel("Alpha Statistical Summary", verbatimTextOutput("alpha_stats_table"))
                    )
                  ),
                  
                  # BETA DIVERSITY TABS
                  conditionalPanel(
                    condition = "input.diversity_type == 'beta'",
                    tabsetPanel(
                      tabPanel("Beta Ordination Plot", plotlyOutput("ordination_plot", height = "600px")),
                      tabPanel("Beta Statistical Summary", verbatimTextOutput("beta_stats_table")),
                      tabPanel("Clustering Dendrogram", plotOutput("clustering_plot", height = "600px"))
                    )
                  )
                )
                
              )
      ),
      
      
      
      # Other possible tabs 
      tabItem(tabName = "clinical", 
              h3("Clinical Factors - Under Development")
              ),
      
      # TAB 4: Lifestyle Tab 
      tabItem(tabName = "multifactor",
              fluidRow(
                box(
                  title = "Parallel Coordinates Plot",
                  status = "primary",
                  solidHeader = TRUE,
                  width = 12,
                  helpText("Visualize the relationship between microbial taxa and selected metadata."),
                  uiOutput("parallel_var_select"),
                  plotlyOutput("parallel_plot", height = "600px")
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
    taxonomy_data = NULL,
    control_sample_metadata = NULL, 
    control_run_metadata = NULL, 
    control_taxonomy_data = NULL, 
    error_message = NULL, 
    success_message = NULL
  )
  
  # Generate example data when button is clicked 
  observeEvent(input$generate_example, {
    
    data_store$error_message <- NULL
    
    # Generate sample metadata 
    sample_metadata <- data.frame(
      Run_ID = paste0("RUN", sprintf("%03d", 1:50)),
      Sample_ID = paste0("SAMPLE", sprintf("%03d", 1:50)),
      Collection_Date = as.character(sample(seq(as.Date('2023-01-01'), as.Date('2025-04-30'), by="day"), 50)),
      Collection_Storage_Temperature = sample(c("-80", "-20", "4", "Room Temperature"), 50, replace = TRUE),
      Analyst_Processor_Name = sample(c("John Doe", "Jane Smith", "Alex Johnson"), 50, replace = TRUE),
      Gender = sample(c("Male", "Female", "Other"), 50, replace = TRUE),
      Age = sample(18:80, 50, replace = TRUE),
      Ongoing_conditions = sample(c("None", "Diabetes", "Hypertension", "IBS", "Crohn's Disease"), 50, replace = TRUE),
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
    
    # Generate taxonomy data 
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
    
    # Generate CONTROL sample metadata with some different characteristics 
    control_sample_metadata <- data.frame(
      Run_ID = paste0("CTRL_RUN", sprintf("%03d", 1:30)),
      Sample_ID = paste0("CTRL_SAMPLE", sprintf("%03d", 1:30)),
      Collection_Date = as.character(sample(seq(as.Date('2023-01-01'), as.Date('2025-04-30'), by="day"), 30)),
      Collection_Storage_Temperature = sample(c("-80", "-20", "4", "Room Temperature"), 30, replace = TRUE),
      Analyst_Processor_Name = sample(c("John Doe", "Jane Smith", "Alex Johnson"), 30, replace = TRUE),
      Gender = sample(c("Male", "Female", "Other"), 30, replace = TRUE),
      Age = sample(18:80, 30, replace = TRUE),
      Ongoing_conditions = "None",  # Controls have no conditions
      Neurological_Disorders = "None",  # Controls have no neurological disorders
      Allergies = sample(c("None", "Pollen", "Food", "Medicine"), 30, replace = TRUE),
      Bowel_movement_frequency = sample(c("Daily", "2-3 times a week", "4-6 times a week"), 30, replace = TRUE),
      Bowel_movement_quality = sample(c("Normal", "Variable"), 30, replace = TRUE, prob = c(0.8, 0.2)),
      Antibiotic_intake = sample(c("None in past year", "Within past 6 months"), 30, replace = TRUE, prob = c(0.9, 0.1)),
      Medications = "None",  # Controls take no medications
      Cancer = "No",  # Controls have no cancer
      Cancer_treatment = "None",  # Controls have no cancer treatment
      BMI = round(rnorm(30, mean = 22, sd = 2), 1),  # Controls have healthier BMI
      Exercise_frequency = sample(c("3-4 times a week", "5+ times a week"), 30, replace = TRUE),  # Controls exercise more
      Exercise_intensity = sample(c("Moderate", "High"), 30, replace = TRUE),
      Smoking_status = sample(c("Never smoked", "Former smoker"), 30, replace = TRUE, prob = c(0.9, 0.1)),
      Cigarettes_per_day = NA,
      Stopped_smoking = NA,
      Alcohol_consumption = sample(c("Never", "Occasionally"), 30, replace = TRUE),
      Alcohol_frequency = NA,
      Drinks_per_day = NA
    )
    
    # Generate CONTROL run metadata 
    control_run_metadata <- data.frame(
      Timestamp = as.character(sample(seq(as.Date('2023-01-01'), as.Date('2025-04-30'), by="day"), 20)),
      Run_ID = paste0("CTRL_RUN", sprintf("%03d", 1:20)),
      Sequencing_Date = as.character(sample(seq(as.Date('2023-01-15'), as.Date('2025-05-01'), by="day"), 20)),
      Sequencing_Platform = sample(c("Illumina", "Ion Torrent"), 20, replace = TRUE),  # Controls use only two platforms
      Sequencing_Type = sample(c("16S rRNA", "Shotgun"), 20, replace = TRUE),  # Controls use only two types
      Expected_read_length = sample(c("250bp", "300bp"), 20, replace = TRUE),
      Sequencing_depth_target = paste0(sample(c("20", "30", "50"), 20, replace = TRUE), "M reads"),
      Library_preparation_kit = sample(c("Nextera XT", "TruSeq Nano"), 20, replace = TRUE),
      Technician_name = sample(c("Maria Garcia", "Thomas Johnson"), 20, replace = TRUE)
    )
    
    # Generate CONTROL taxonomy data with different abundances 
    control_taxonomy_data <- data.frame(
      Sample_ID = rep(paste0("CTRL_SAMPLE", sprintf("%03d", 1:30)), each = length(species_names)),
      Species = rep(species_names, 30),
      Abundance = round(runif(30 * length(species_names), min = 0, max = 10000))
    )
    
    # Modify control taxonomy to have more of the "healthy" microbiome species
    healthy_indices <- which(control_taxonomy_data$Species %in% c(
      "Faecalibacterium prausnitzii", "Akkermansia muciniphila", 
      "Bifidobacterium longum", "Lactobacillus acidophilus"
    ))
    
    # Increase abundance of healthy species
    control_taxonomy_data$Abundance[healthy_indices] <- 
      control_taxonomy_data$Abundance[healthy_indices] * sample(seq(1.5, 3, 0.1), length(healthy_indices), replace = TRUE)
    
    # Decrease abundance of potentially problematic species
    problematic_indices <- which(control_taxonomy_data$Species %in% c(
      "Clostridium difficile", "Bacteroides fragilis"
    ))
    control_taxonomy_data$Abundance[problematic_indices] <- 
      control_taxonomy_data$Abundance[problematic_indices] * sample(seq(0.1, 0.4, 0.05), length(problematic_indices), replace = TRUE)
    
    # Ensure abundances don't exceed our max range
    control_taxonomy_data$Abundance <- pmin(control_taxonomy_data$Abundance, 10000)
    
    # Store the data in reactive values 
    data_store$sample_metadata <- sample_metadata
    data_store$run_metadata <- run_metadata
    data_store$taxonomy_data <- taxonomy_data 
    data_store$control_sample_metadata <- control_sample_metadata
    data_store$control_run_metadata <- control_run_metadata
    data_store$control_taxonomy_data <- control_taxonomy_data
    data_store$success_message <- "Example data generated successfully! Ready for analysis."
    
    # Show notification 
    showNotification("Example data generated successfully!", type = "message")
  })
  
  # Load data when button is clicked 
  observeEvent(input$load_data, {
    # Clear previous messages 
    data_store$error_message = NULL
    data_store$success_message = NULL
    
    ##### HANDLE POTENTIAL ERRORS ####
    tryCatch({
      # Check if all required files are uploaded 
      if (is.null(input$sample_metadata) || is.null(input$run_metadata) || is.null(input$taxonomy_file) ||
          is.null(input$control_sample_metadata) || is.null(input$control_run_metadata) || is.null(input$control_taxonomy)) {
        data_store$error_message <- "Please upload all required files"
        return(NULL)
      }
      
      # Read sample metadata 
      sample_metadata <- read.csv(input$sample_metadata$datapath, header = TRUE)
      
      # Read run metadata 
      run_metadata <- read.csv(input$run_metadata$datapath, header = TRUE)
      
      # Read taxonomy data 
      taxonomy_data <- read.delim(input$taxonomy_file$datapath, sep = "\t", header = TRUE)
      
      # Read control sample metadata 
      control_sample_metadata <- read.csv(input$control_sample_metadata$datapath, header = TRUE)
      
      # Read control run metadata 
      control_run_metadata <- read.csv(input$control_run_metadata$datapath, header = TRUE)
      
      # Read control taxonomy data 
      control_taxonomy_data <- read.delim(input$control_taxonomy$datapath, sep = "\t", header = TRUE)
      
      # Store the data in reactive values 
      data_store$sample_metadata <- sample_metadata
      data_store$run_metadata <- run_metadata
      data_store$taxonomy_data <- taxonomy_data
      data_store$control_sample_metadata <- control_sample_metadata
      data_store$control_run_metadata <- control_run_metadata
      data_store$control_taxonomy_data <- control_taxonomy_data
      
      # Show notification 
      showNotification("Data loaded successfully!", type = "message")
    }, error = function(e) {
      # Handle any error that occur during file reading
      data_store$error_message <- paste("Error loading data:", e$message)
      showNotification(paste("Error:", e$message), type = "error")
    })
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
  
  # Add a control samples value box
  output$total_control_samples_box <- renderValueBox({
    req(data_store$control_sample_metadata)
    valueBox(
      nrow(data_store$control_sample_metadata), 
      "Control Samples", 
      icon = icon("vial"), 
      color = "yellow"
    )
  })
  
  ###############################################################
  ##### FUTURE: Handle missing values in Dynamic summary ########
  ###############################################################
  
  
  # Dynamic summary content 
  output$data_summary_dynamic <- renderUI({
    req(data_store$sample_metadata, data_store$run_metadata, data_store$taxonomy_data)
    
    # General Overview 
    total_samples <- nrow(data_store$sample_metadata)
    total_runs <- nrow(data_store$run_metadata)
    total_species <- length(unique(data_store$taxonomy_data$Species))
    sample_date_range <- paste(min(as.Date(data_store$sample_metadata$Collection_Date)), "to", 
                               max(as.Date(data_store$sample_metadata$Collection_Date)))
    seq_date_range <- paste(min(as.Date(data_store$run_metadata$Sequencing_Date)), "to", 
                            max(as.Date(data_store$run_metadata$Sequencing_Date)))
    
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
    # 
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
  
  ############################################################## 
  
  ## Filter by 1 sample 
  # Update sample selection dropdown 
  observe({
    req(data_store$sample_metadata)
    sample_choices <- data_store$sample_metadata$Sample_ID
    names(sample_choices) <- paste0(data_store$sample_metadata$Sample_ID, 
                                    " (", data_store$sample_metadata$Collection_Date, ")")
    
    updateSelectizeInput(session, "selected_sample_id", 
                         choices = c("All Samples" = "", sample_choices), 
                         selected = "")
  })
  
  # Create reactive for filtered data 
  filtered_data <- reactive({
    req(data_store$sample_metadata, data_store$taxonomy_data)
    
    selected_sample <- input$selected_sample_id
    
    if (selected_sample == "" || is.null(selected_sample)) {
      return(list(
        sample_metadata = data_store$sample_metadata, 
        taxonomy_data = data_store$taxonomy_data, 
        filtered = FALSE
      ))
    } else {
      # Filter the metadata to just selected sample
      filtered_sample_metadata <- data_store$sample_metadata 
      if ("Sample_ID" %in% colnames(filtered_sample_metadata)) {
        filtered_sample_metadata <- filtered_sample_metadata[filtered_sample_metadata$Sample_ID == selected_sample, ]
      } else if ("SampleID" %in% colnames(filtered_sample_metadata)) {
        filtered_sample_metadata <- filtered_sample_metadata[filtered_sample_metadata$SampleID == selected_sample, ]
      }
      
      # Filter taxonomy data to just selected sample
      filtered_taxonomy_data <- data_store$taxonomy_data 
      if ("Sample_ID" %in% colnames(filtered_taxonomy_data)) {
        filtered_taxonomy_data <- filtered_taxonomy_data[filtered_taxonomy_data$Sample_ID == selected_sample, ]
      } else if ("SampleID" %in% colnames(filtered_taxonomy_data)) {
        filtered_taxonomy_data <- filtered_taxonomy_data[filtered_taxonomy_data$SampleID == selected]
      }

      return(list(
        sample_metadata = filtered_sample_metadata,
        taxonomy_data = filtered_taxonomy_data, 
        filtered = TRUE
      ))
    }
  })
  
  # Render information about the selected sample
  output$selected_sample_info <- renderUI({
    req(filtered_data())
    
    if (!filtered_data()$filtered) {
      return(HTML("<p><i>No specific sample selected. Showing data for all samples.</i></p>"))
    }
    sample_data <- filtered_data()$sample_metadata
    
    # Extract relevant information 
    HTML(paste0(
      "<h4>Selected Sample Information</h4>",
      "<table class='table table-condensed table-bordered'>",
      "<tr><th>Sample ID</th><td>", sample_data$Sample_ID, "</td></tr>",
      "<tr><th>Collection Date</th><td>", sample_data$Collection_Date, "</td></tr>",
      "<tr><th>Gender</th><td>", sample_data$Gender, "</td></tr>",
      "<tr><th>Age</th><td>", sample_data$Age, "</td></tr>",
      "<tr><th>BMI</th><td>", sample_data$BMI, "</td></tr>",
      "<tr><th>Ongoing Conditions</th><td>", sample_data$Ongoing_conditions, "</td></tr>",
      "<tr><th>Collection Storage</th><td>", sample_data$Collection_Storage_Temperature, "</td></tr>",
      "</table>"
    ))
  })
  

  # Render metadata table previews with filtering 
  output$sample_metadata_preview <- renderDataTable({
    req(data_store$sample_metadata)
    
    # Get filter condition
    filter_value <- input$sample_metadata_filter 
    selected_sample <- input$selected_sample_id
    
    # Start with base data 
    sample_data <- data_store$sample_metadata
    
    # Apply sample filter if selected 
    if (selected_sample != "" && !is.null(selected_sample)) {
      sample_data <- sample_data %>% filter(Sample_ID == selected_sample)
    }
    
    # Filter data based on condition selection 
    if (!is.null(filter_value)) {
      if (filter_value == "Control") {
        req(data_store$control_sample_metadata)
        return(datatable(data_store$control_sample_metadata, 
                         options = list(pageLength = 10, scrollX = TRUE)))
      } else if (filter_value == "Sample"){
        return(datatable(data_store$sample_metadata, 
                         options = list(pageLength = 10, scrollX = TRUE)))
      } else { # All
        # Combine both datasets if control data exists 
        if (!is.null(data_store$control_sample_metadata)) {
          # Add a column to indicate the source 
          sample_data <- data_store$sample_metadata
          sample_data$Data_Source <- "Sample"
          
          control_data <- data_store$control_sample_metadata
          control_data$Data_Source <- "Control"
          
          # Combine the datasets (assume they have compatible columns)
          combined_data <- rbind(sample_data, control_data)
          return(datatable(combined_data, options = list(pageLength = 10, scrollX = TRUE)))
        } else {
          return(datatable(data_store$sample_metadata, 
                           options = list(pageLength = 10, scrollX = TRUE)))
        }
      }
    } else {
      # Default display (all data)
      return(datatable(data_store$sample_metadata, 
                       options = list(pageLength = 10, scrollX = TRUE)))
    }
  })
  
  output$run_metadata_preview <- renderDataTable({
    req(data_store$run_metadata)
    
    # Get filter condition 
    filter_value <- input$run_metadata_filter
    
    # Filter data based on condition selection
    if (!is.null(filter_value)) {
      if (filter_value == "Control") {
        req(data_store$control_run_metadata) 
        return(datatable(data_store$control_run_metadata, 
                         options = list(pageLength = 10, scrollX = TRUE)))
      } else if (filter_value == "Sample") {
        return(datatable(data_store$run_metadata, 
                          options = list(pageLength = 10, scrollX = TRUE)))
      } else {
        # Combine both datasets 
        if (!is.null(data_store$control_run_metadata)) {
          # Add a column to indicate source 
          sample_data <- data_store$run_metadata
          sample_data$Data_Source <- "Sample"
          
          control_data <- data_store$control_run_metadata
          control_data$Data_Source <- "Control"
          
          # Combine datasets (assume they have compatible columns)
          combined_data <- rbind(sample_data, control_data)
          return(datatable(combined_data, options = list(pageLength = 10, scrollX = TRUE)))
        } else {
          return(datatable(data_store$run_metadata, 
                           options = list(pageLength = 10, scrollX = TRUE)))
        }
      } 
    } else {
      return( datatable(data_store$run_metadata, 
              options = list(pageLength = 10, scrollX = TRUE)))
    }
  })
  
  output$taxonomy_data_preview <- renderDataTable({
    req(data_store$taxonomy_data)
    
    # Get filter condition 
    filter_value <- input$taxonomy_data_filter
    selected_sample <- input$selected_sample_id
    
    # Start with base data 
    taxonomy_data <- data_store$taxonomy_data 
    
    # Apply sample filter if selected
    if (selected_sample != "" && !is.null(selected_sample)) {
      taxonomy_data <- taxonomy_data %>% filter(Sample_ID == selected_sample)
    }
    
    # Filter data based on selection
    if (!is.null(filter_value)) {
      if (filter_value == "Control") {
        req(data_store$control_taxonomy_data)
        return(datatable(data_store$control_taxonomy_data, 
                         options = list(pageLength = 10, scrollX = TRUE)))
      } else if (filter_value == "Sample") {
        return(datatable(data_store$taxonomy_data, 
                         options = list(pageLength = 10, scrollX = TRUE)))
      } else {
        # Combine both datasets
        if (!is.null(data_store$control_taxonomy_data)) {
          # Add a column to indicate the source
          sample_data <- data_store$taxonomy_data
          sample_data$Data_Source <- "Sample"
          
          control_data <- data_store$control_taxonomy_data
          control_data$Data_Source <- "Control"
          
          # Combine the datasets (assuming they have compatible columns)
          combined_data <- rbind(sample_data, control_data)
          return(datatable(combined_data, options = list(pageLength = 10, scrollX = TRUE)))
        } else {
          return(datatable(data_store$taxonomy_data, 
                           options = list(pageLength = 10, scrollX = TRUE)))
        }
      }
    } else {
      return(datatable(data_store$taxonomy_data, 
              options = list(pageLength = 10, scrollX = TRUE)))
    }
  })
  
  
  
 
  
  #### TAXONOMY TAB ######
  
  # Define colors 
  taxonomic_colors <- list(
    "Viridis" = scale_fill_viridis_d(),
    "Set1" = scale_fill_brewer(palette = "Set1"),
    "Set2" = scale_fill_brewer(palette = "Set2"),
    "Paired" = scale_fill_brewer(palette = "Paired"),
    "Default" = scale_fill_discrete()
  )
  
  # Function to extract taxonomic levels from species names
  extract_taxonomy <- reactive({
    req(data_store$taxonomy_data, data_store$control_taxonomy_data)
    
    # Combine patient and control taxonomy data 
    all_taxonomy <- rbind(data_store$taxonomy_data, data_store$control_taxonomy_data)
    
    ### SIMPLIFIED APPROACH. Real data needs to parse proper taxonomic classification
    
    # Get unique species 
    species_list <- unique(all_taxonomy$Species)
    
    # Create a dataframe to store taxonomic hierarchy 
    taxonomy_df <- data.frame(Species = species_list)
    
    # Extract Genus 
    taxonomy_df$Genus <- sapply(strsplit(species_list, " "), function(x) x[1])
    
    # Create Family 
    genera <- unique(taxonomy_df$Genus)
    families <- paste0("Family_", ceiling(seq_along(genera)/2))
    family_map <- setNames(families[1:length(genera)], genera)
    taxonomy_df$Family <- family_map[taxonomy_df$Genus]
    
    # Create Order 
    unique_families <- unique(taxonomy_df$Family)
    orders <- paste0("Order_", ceiling(seq_along(unique_families)/3))
    order_map <- setNames(orders[1:length(unique_families)], unique_families)
    taxonomy_df$Order <- order_map[taxonomy_df$Family]
    
    # Create Class 
    unique_orders <- unique(taxonomy_df$Order)
    classes <- paste0("Class_", ceiling(seq_along(unique_orders)/2))
    class_map <- setNames(classes[1:length(unique_orders)], unique_orders)
    taxonomy_df$Class <- class_map[taxonomy_df$Order]
    
    # Create Phylum 
    unique_classes <- unique(taxonomy_df$Class)
    phyla <- paste0("Phylum_", ceiling(seq_along(unique_classes)/2))
    phylum_map <- setNames(phyla[1:length(unique_classes)], unique_classes)
    taxonomy_df$Phylum <- phylum_map[taxonomy_df$Class]
    
    return(taxonomy_df)
  })
  
  # Update sample selection dropdown when data is available 
  observe({
    req(data_store$sample_metadata, data_store$control_sample_metadata)
    
    # Combine both patient and control sample IDs
    all_samples <- c(data_store$sample_metadata$Sample_ID, 
                     data_store$control_sample_metadata$Sample_ID)
    
    # Update the selectize input with available samples
    updateSelectizeInput(session, "selected_sample_id", 
                         choices = all_samples,
                         selected = NULL)
  })
  
  # Filter taxonomy data based on user selections 
  # Filter taxonomy data based on user selections
  filtered_taxonomy_data <- reactive({
    req(data_store$taxonomy_data, data_store$control_taxonomy_data)
    
    # Apply sample group filter
    if (input$sample_group == "Patients Only") {
      filtered_data <- data_store$taxonomy_data
    } else if (input$sample_group == "Controls Only") {
      filtered_data <- data_store$control_taxonomy_data
    } else { # All Samples
      filtered_data <- rbind(data_store$taxonomy_data, data_store$control_taxonomy_data)
    }
    
    # Apply specific sample filter if selected
    if (!is.null(input$selected_sample_id) && input$selected_sample_id != "") {
      filtered_data <- filtered_data[filtered_data$Sample_ID == input$selected_sample_id, ]
    }
    
    # Merge with taxonomic hierarchy
    taxonomy_hierarchy <- extract_taxonomy()
    filtered_data <- merge(filtered_data, taxonomy_hierarchy, by = "Species")
    
    # Apply abundance threshold filter (will be used in visualization)
    
    return(filtered_data)
  })
  
  # Function to get merged metadata (both patient and control)
  merged_metadata <- reactive({
    req(data_store$sample_metadata, data_store$control_sample_metadata)
    
    # Combine patient and control metadata
    combined_metadata <- rbind(data_store$sample_metadata, data_store$control_sample_metadata)
    
    # Add a group column to distinguish patients from controls
    combined_metadata$Group <- ifelse(grepl("^CTRL_", combined_metadata$Sample_ID), "Control", "Patient")
    
    return(combined_metadata)
  })
  
  # Calculate relative abundance for taxonomy visualization
  taxonomy_relative_abundance <- reactive({
    req(filtered_taxonomy_data())
    
    # Get filtered data
    tax_data <- filtered_taxonomy_data()
    
    # Group by Sample_ID to calculate total abundance per sample
    sample_totals <- aggregate(Abundance ~ Sample_ID, data = tax_data, FUN = sum)
    
    # Merge with original data to calculate relative abundance
    tax_data_with_totals <- merge(tax_data, sample_totals, by = "Sample_ID", 
                                  suffixes = c("", "_total"))
    
    # Calculate relative abundance
    tax_data_with_totals$RelativeAbundance <- (tax_data_with_totals$Abundance / 
                                                 tax_data_with_totals$Abundance_total) * 100
    
    # Filter by abundance threshold
    filtered_by_threshold <- tax_data_with_totals[tax_data_with_totals$RelativeAbundance >= input$abundance_threshold, ]
    
    # Get selected taxonomic level
    tax_level <- input$tax_level
    
    # Group by taxonomic level and aggregate abundances
    if (tax_level != "Species") {
      # Aggregate by the selected taxonomic level
      agg_formula <- as.formula(paste("cbind(Abundance, RelativeAbundance) ~", 
                                      "Sample_ID +", tax_level))
      
      aggregated_data <- aggregate(agg_formula, 
                                   data = filtered_by_threshold, 
                                   FUN = sum)
      
      # Replace column names after aggregation
      colnames(aggregated_data)[colnames(aggregated_data) == tax_level] <- "TaxonomicGroup"
    } else {
      # For species level, just rename the column
      aggregated_data <- filtered_by_threshold
      aggregated_data$TaxonomicGroup <- aggregated_data$Species
      
      # Handle species with low abundance
      if (input$abundance_threshold > 0) {
        sorted_species <- unique(aggregated_data[order(-aggregated_data$RelativeAbundance), "TaxonomicGroup"])
        keep_top_n <- min(15, length(sorted_species))  # Limit to top 15 species
        
        if (length(sorted_species) > keep_top_n) {
          top_species <- sorted_species[1:keep_top_n]
          aggregated_data$TaxonomicGroup <- ifelse(aggregated_data$TaxonomicGroup %in% top_species, 
                                                   as.character(aggregated_data$TaxonomicGroup), 
                                                   "Other")
        }
      }
    }
    # Sort by abundance if requested
    if (input$sort_abundance) {
      # Aggregate by taxonomic group to get total abundance
      group_totals <- aggregate(Abundance ~ TaxonomicGroup, data = aggregated_data, FUN = sum)
      group_totals <- group_totals[order(-group_totals$Abundance), ]
      
      # Create an ordered factor for taxonomic groups based on abundance
      aggregated_data$TaxonomicGroup <- factor(aggregated_data$TaxonomicGroup, 
                                               levels = group_totals$TaxonomicGroup)
    }
    
    return(aggregated_data)
  })
  
  # Add a group column to taxonomy data for comparison
  taxonomy_with_groups <- reactive({
    req(taxonomy_relative_abundance(), merged_metadata())
    
    tax_data <- taxonomy_relative_abundance()
    metadata <- merged_metadata()
    
    # Merge taxonomy data with metadata to get group information
    result <- merge(tax_data, metadata[, c("Sample_ID", "Group")], by = "Sample_ID", all.x = TRUE)
    
    # If metadata_group is selected, add that grouping as well
    if (input$metadata_group != "None") {
      metadata_subset <- metadata[, c("Sample_ID", input$metadata_group)]
      result <- merge(result, metadata_subset, by = "Sample_ID", all.x = TRUE)
    }
    
    return(result)
  })
  
  # Render the quick stat boxes
  output$total_taxa_box <- renderValueBox({
    req(filtered_taxonomy_data())
    
    # Count unique taxa at the selected taxonomic level
    tax_level <- input$tax_level
    unique_taxa <- length(unique(filtered_taxonomy_data()[[tax_level]]))
    
    valueBox(
      value = unique_taxa,
      subtitle = paste("Total", tax_level, "Identified"),
      icon = icon("microbe"),
      color = "purple"
    )
  })
  
  output$dominant_taxa_box <- renderValueBox({
    req(taxonomy_relative_abundance())
    
    # Find dominant taxa (with highest average relative abundance)
    taxa_avg <- aggregate(RelativeAbundance ~ TaxonomicGroup, 
                          data = taxonomy_relative_abundance(), 
                          FUN = mean)
    
    dominant_taxa <- taxa_avg$TaxonomicGroup[which.max(taxa_avg$RelativeAbundance)]
    
    valueBox(
      value = as.character(dominant_taxa),
      subtitle = paste("Dominant", input$tax_level),
      icon = icon("crown"),
      color = "olive"
    )
  })
  
  output$unique_taxa_box <- renderValueBox({
    req(filtered_taxonomy_data())
    
    # Count samples
    sample_count <- length(unique(filtered_taxonomy_data()$Sample_ID))
    
    valueBox(
      value = sample_count,
      subtitle = "Samples Analyzed",
      icon = icon("vials"),
      color = "maroon"
    )
  })
  
  # Get the selected color scheme
  get_color_scheme <- reactive({
    req(input$abundance_color_scheme)
    
    # Return the appropriate color scale based on selection
    if (input$abundance_color_scheme %in% names(taxonomic_colors)) {
      return(taxonomic_colors[[input$abundance_color_scheme]])
    } else {
      return(taxonomic_colors[["Default"]])
    }
  })
  
  # Generate the relative abundance plot
  output$relative_abundance_plot <- renderPlot({
    req(taxonomy_with_groups())
    
    # Get data with groups
    plot_data <- taxonomy_with_groups()
    
    # Base plot with taxonomic groups
    p <- ggplot(plot_data, aes(x = Sample_ID, y = RelativeAbundance, fill = TaxonomicGroup)) +
      geom_bar(stat = "identity") +
      theme_minimal() +
      theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 8),
            strip.text = element_text(size = 10, face = "bold"),
            legend.position = "bottom",
            legend.text = element_text(size = 8)) +
      labs(x = "Sample ID", y = "Relative Abundance (%)", 
           title = paste("Relative Abundance of", input$tax_level),
           fill = input$tax_level) +
      get_color_scheme()
    
    # Determine if we need to split by a metadata group
    if (input$metadata_group != "None") {
      # Get the metadata column name
      metadata_col <- input$metadata_group
      
      # Create plot based on metadata grouping
      if (input$compare_groups && "Group" %in% colnames(plot_data)) {
        # Plot with both metadata group and patient/control group
        p <- p + facet_grid(reformulate("Group", metadata_col), scales = "free_x", space = "free") +
          labs(title = paste("Relative Abundance of", input$tax_level, "by", metadata_col, "and Group"))
      } else {
        # Plot with only metadata grouping
        p <- p + facet_wrap(reformulate(metadata_col), scales = "free_x") +
          labs(title = paste("Relative Abundance of", input$tax_level, "by", metadata_col))
      }
    } else if (input$compare_groups && "Group" %in% colnames(plot_data)) {
      # Plot with only patient/control grouping
      p <- p + facet_wrap(~ Group, scales = "free_x") +
        labs(title = paste("Relative Abundance of", input$tax_level, "by Group"))
    }
    
    return(p)
  })
  
  # Download handler for the abundance plot
  output$download_abundance_plot <- downloadHandler(
    filename = function() {
      paste("taxonomy_", input$tax_level, "_plot_", 
            format(Sys.time(), "%Y%m%d_%H%M%S"), ".png", sep = "")
    },
    content = function(file) {
      # Capture the current plot data
      plot_data <- taxonomy_with_groups()
      
      # Set up the plot with the same parameters as the displayed plot
      p <- ggplot(plot_data, aes(x = Sample_ID, y = RelativeAbundance, fill = TaxonomicGroup)) +
        geom_bar(stat = "identity") +
        theme_minimal() +
        theme(axis.text.x = element_text(angle = 90, hjust = 1),
              legend.position = "bottom") +
        labs(x = "Sample ID", y = "Relative Abundance (%)", 
             title = paste("Relative Abundance of", input$tax_level),
             fill = input$tax_level) +
        get_color_scheme()
      
      # Add faceting if needed
      if (input$metadata_group != "None" && input$compare_groups) {
        p <- p + facet_grid(reformulate("Group", input$metadata_group), scales = "free_x", space = "free")
      } else if (input$metadata_group != "None") {
        p <- p + facet_wrap(reformulate(input$metadata_group), scales = "free_x")
      } else if (input$compare_groups) {
        p <- p + facet_wrap(~ Group, scales = "free_x")
      }
      
      # Save the plot to the file
      ggsave(file, plot = p, width = 10, height = 8, dpi = 300)
    }
  )
  
  ###### DIVERSITY ANALYSIS TAB ######
  output$alpha_plot_type_ui <- renderUI({
    req(input$alpha_metric)
    
    if (input$alpha_metric == "Observed OTUs") {
      selectInput("alpha_plot_type", "Plot Type:",
                  choices = c("Violin", "Jittered Scatterplot"),
                  selected = "Violin")
    } else {
      selectInput("alpha_plot_type", "Plot Type:",
                  choices = c("Boxplot + Dots", "Violin", "Jittered Scatterplot"),
                  selected = "Boxplot + Dots")
    }
  })
  
  
  output$diversity_main_plot <- renderPlotly({
    req(data_store$sample_metadata, data_store$taxonomy_data)
    req(data_store$control_sample_metadata, data_store$control_taxonomy_data)
    req(input$diversity_type == "alpha")
    req(input$alpha_metric, input$alpha_plot_type, input$alpha_color_palette)
    
    set.seed(42)
    jitter_pos <- position_jitter(width = 0.2, height = 0)
    
    sample_meta <- data_store$sample_metadata
    control_meta <- data_store$control_sample_metadata
    
    sample_tax <- data_store$taxonomy_data %>% filter(Sample_ID %in% sample_meta$Sample_ID)
    control_tax <- data_store$control_taxonomy_data %>% filter(Sample_ID %in% control_meta$Sample_ID)
    
    sample_wide <- sample_tax %>%
      select(Sample_ID, Species, Abundance) %>%
      pivot_wider(names_from = Species, values_from = Abundance, values_fill = 0)
    
    sample_div <- sample_wide %>%
      column_to_rownames("Sample_ID") %>%
      as.matrix() %>%
      {
        tibble(
          Sample_ID = rownames(.),
          Observed = rowSums(. > 0),
          Shannon = vegan::diversity(., index = "shannon"),
          Simpson = vegan::diversity(., index = "simpson")
        )
      }
    
    if (input$subdivide_samples) {
      req(input$condition_column)
      selected_col <- input$condition_column
      sample_div <- left_join(sample_div,
                              sample_meta %>% select(Sample_ID, !!sym(selected_col)),
                              by = "Sample_ID")
      sample_div$GroupLabel <- sample_div[[selected_col]]
    } else {
      sample_div$GroupLabel <- "Sample"
    }
    
    control_wide <- control_tax %>%
      select(Sample_ID, Species, Abundance) %>%
      pivot_wider(names_from = Species, values_from = Abundance, values_fill = 0)
    
    control_div <- control_wide %>%
      column_to_rownames("Sample_ID") %>%
      as.matrix() %>%
      {
        tibble(
          Sample_ID = rownames(.),
          Observed = rowSums(. > 0),
          Shannon = vegan::diversity(., index = "shannon"),
          Simpson = vegan::diversity(., index = "simpson"),
          GroupLabel = "Control"
        )
      }
    
    diversity_df <- bind_rows(sample_div, control_div)
    
    group_levels <- c("Control", sort(unique(diversity_df$GroupLabel[diversity_df$GroupLabel != "Control"])))
    diversity_df$GroupLabel <- factor(diversity_df$GroupLabel, levels = group_levels)
    
    metric_col <- switch(input$alpha_metric,
                         "Observed OTUs" = "Observed",
                         "Shannon" = "Shannon",
                         "Simpson" = "Simpson")
    
    p <- ggplot(diversity_df, aes(x = GroupLabel, y = .data[[metric_col]], fill = GroupLabel))
    
    if (input$alpha_plot_type == "Boxplot + Dots") {
      p <- p +
        geom_boxplot(alpha = 0.6, outlier.shape = NA) +
        geom_jitter(aes(text = Sample_ID), shape = 21, alpha = 0.6, color = "black", position = jitter_pos)
    } else if (input$alpha_plot_type == "Violin") {
      p <- p + geom_violin(alpha = 0.7, trim = FALSE)
    } else if (input$alpha_plot_type == "Jittered Scatterplot") {
      p <- p + geom_jitter(aes(text = Sample_ID), shape = 21, alpha = 0.7, color = "black", position = jitter_pos)
    }
    
    if (!is.null(input$selected_sample_id) && input$selected_sample_id %in% diversity_df$Sample_ID) {
      highlighted <- diversity_df %>% filter(Sample_ID == input$selected_sample_id)
      p <- p + geom_point(data = highlighted,
                          aes(x = GroupLabel, y = .data[[metric_col]]),
                          shape = 21, size = 4, stroke = 1.5,
                          color = "red", fill = "yellow")
    }
    
    p <- p +
      labs(title = paste(input$alpha_metric, "Diversity across Groups"),
           x = "Group", y = paste(input$alpha_metric, "Index")) +
      theme_minimal() +
      theme(legend.position = ifelse(input$show_diversity_legend, "right", "none")) +
      scale_fill_brewer(palette = input$alpha_color_palette)
    
    ggplotly(p, tooltip = c("x", "y", "text")) %>%
      layout(dragmode = "zoom")
  })
  
  
  
  output$alpha_stats_table <- renderPrint({
    req(data_store$sample_metadata, data_store$taxonomy_data)
    req(data_store$control_sample_metadata, data_store$control_taxonomy_data)
    
    # Get metadata and taxonomy
    sample_meta <- data_store$sample_metadata
    control_meta <- data_store$control_sample_metadata
    
    sample_tax <- data_store$taxonomy_data %>% filter(Sample_ID %in% sample_meta$Sample_ID)
    control_tax <- data_store$control_taxonomy_data %>% filter(Sample_ID %in% control_meta$Sample_ID)
    
    # Sample diversity
    sample_wide <- sample_tax %>%
      select(Sample_ID, Species, Abundance) %>%
      pivot_wider(names_from = Species, values_from = Abundance, values_fill = 0)
    
    sample_div <- sample_wide %>%
      column_to_rownames("Sample_ID") %>%
      as.matrix() %>%
      {
        tibble(
          Sample_ID = rownames(.),
          Observed = rowSums(. > 0),
          Shannon = vegan::diversity(., index = "shannon"),
          Simpson = vegan::diversity(., index = "simpson")
        )
      }
    
    # Add GroupLabel from condition or default
    if (input$subdivide_samples) {
      req(input$condition_column)
      selected_col <- input$condition_column
      sample_div <- left_join(sample_div,
                              sample_meta %>% select(Sample_ID, !!sym(selected_col)),
                              by = "Sample_ID")
      sample_div$GroupLabel <- sample_div[[selected_col]]
    } else {
      sample_div$GroupLabel <- "Sample"
    }
    
    # Control diversity
    control_wide <- control_tax %>%
      select(Sample_ID, Species, Abundance) %>%
      pivot_wider(names_from = Species, values_from = Abundance, values_fill = 0)
    
    control_div <- control_wide %>%
      column_to_rownames("Sample_ID") %>%
      as.matrix() %>%
      {
        tibble(
          Sample_ID = rownames(.),
          Observed = rowSums(. > 0),
          Shannon = vegan::diversity(., index = "shannon"),
          Simpson = vegan::diversity(., index = "simpson"),
          GroupLabel = "Control"
        )
      }
    
    # Combine
    diversity_df <- bind_rows(sample_div, control_div) %>%
      filter(!is.na(GroupLabel))  # Clean up NAs
    
    # Choose metric
    metric_col <- switch(input$alpha_metric,
                         "Observed OTUs" = "Observed",
                         "Shannon" = "Shannon",
                         "Simpson" = "Simpson")
    
    # Check group count
    n_groups <- dplyr::n_distinct(diversity_df$GroupLabel)
    if (n_groups < 2) {
      cat("Not enough groups to compare.")
      return(NULL)
    }
    
    cat("\n")
    
    if (n_groups == 2) {
      cat("✅ Only 2 groups detected — using Wilcoxon rank-sum test (non-parametric):\n")
      wilcox_result <- wilcox.test(as.formula(paste(metric_col, "~ GroupLabel")), data = diversity_df)
      print(wilcox_result)
      
    } else {
      cat("🔍 Pairwise comparisons (Dunn’s test, BH corrected):\n")
      if (!requireNamespace("FSA", quietly = TRUE)) {
        cat("❌ Please install the 'FSA' package to enable Dunn's test.\n")
      } else {
        suppressPackageStartupMessages(library(FSA))
        dunn_res <- FSA::dunnTest(
          as.formula(paste(metric_col, "~ GroupLabel")),
          data = diversity_df,
          method = "bh"
        )
        print(dunn_res)
      }
    }
  })
  
  output$ordination_plot <- renderPlotly({
    req(input$diversity_type == "beta")
    req(data_store$sample_metadata, data_store$taxonomy_data)
    req(data_store$control_sample_metadata, data_store$control_taxonomy_data)
    req(input$distance_metric, input$ordination_method)
    
    # Combine sample and control taxonomy
    all_tax <- bind_rows(
      data_store$taxonomy_data,
      data_store$control_taxonomy_data
    ) %>%
      select(Sample_ID, Species, Abundance) %>%
      pivot_wider(names_from = Species, values_from = Abundance, values_fill = 0) %>%
      column_to_rownames("Sample_ID")
    
    # Calculate distance matrix
    dist_method <- switch(input$distance_metric,
                          "Bray-Curtis" = "bray",
                          "Euclidean"   = "euclidean",
                          "Jaccard"     = "jaccard")
    dist_matrix <- vegan::vegdist(all_tax, method = dist_method)
    
    # Perform ordination
    ord_result <- if (input$ordination_method == "PCoA") {
      ape::pcoa(dist_matrix)$vectors[, 1:2]
    } else {
      MASS::isoMDS(dist_matrix, k = 2)$points
    }
    
    ord_df <- as.data.frame(ord_result)
    colnames(ord_df)[1:2] <- c("Axis1", "Axis2")  
    ord_df$Sample_ID <- rownames(ord_result)
    
    
    # Add grouping
    all_meta <- bind_rows(data_store$sample_metadata, data_store$control_sample_metadata)
    ord_df <- left_join(ord_df, all_meta, by = "Sample_ID")
    
    # Assign group label
    ord_df$GroupLabel <- if (input$subdivide_samples) {
      req(input$condition_column)
      ord_df[[input$condition_column]]
    } else {
      ifelse(ord_df$Sample_ID %in% data_store$control_sample_metadata$Sample_ID, "Control", "Sample")
    }
    
    # Plot
    p <- ggplot(ord_df, aes(x = Axis1, y = Axis2, color = GroupLabel)) +
      geom_point(size = 3, alpha = 0.8) +
      labs(x = "Axis 1", y = "Axis 2",
           title = paste(input$ordination_method, "on", input$distance_metric, "distance")) +
      theme_minimal()
    
    
    if (input$show_group_ellipses) {
      p <- p + stat_ellipse(type = "norm", linetype = "dashed", alpha = 0.4)
    }
    
    ggplotly(p)
  })
  
  output$beta_stats_table <- renderPrint({
    req(input$diversity_type == "beta")
    req(data_store$sample_metadata, data_store$taxonomy_data)
    req(data_store$control_sample_metadata, data_store$control_taxonomy_data)
    
    # Combine taxonomy
    all_tax <- bind_rows(
      data_store$taxonomy_data,
      data_store$control_taxonomy_data
    ) %>%
      select(Sample_ID, Species, Abundance) %>%
      pivot_wider(names_from = Species, values_from = Abundance, values_fill = 0) %>%
      column_to_rownames("Sample_ID")
    
    # Combine metadata
    all_meta <- bind_rows(data_store$sample_metadata, data_store$control_sample_metadata)
    group_label <- if (input$subdivide_samples) {
      req(input$condition_column)
      all_meta[[input$condition_column]]
    } else {
      ifelse(all_meta$Sample_ID %in% data_store$control_sample_metadata$Sample_ID, "Control", "Sample")
    }
    
    adonis_result <- vegan::adonis2(all_tax ~ group_label, method = tolower(input$distance_metric))
    print(adonis_result)
  })
  
  
  output$clustering_plot <- renderPlot({
    req(input$diversity_type == "beta")
    req(data_store$sample_metadata, data_store$taxonomy_data)
    req(data_store$control_sample_metadata, data_store$control_taxonomy_data)
    req(input$show_clustering_dendrogram)
    
    all_tax <- bind_rows(
      data_store$taxonomy_data,
      data_store$control_taxonomy_data
    ) %>%
      select(Sample_ID, Species, Abundance) %>%
      pivot_wider(names_from = Species, values_from = Abundance, values_fill = 0) %>%
      column_to_rownames("Sample_ID")
    
    dist_matrix <- vegan::vegdist(all_tax, method = tolower(input$distance_metric))
    hc <- hclust(dist_matrix, method = "ward.D2")
    
    plot(hc, main = paste("Hierarchical Clustering (", input$distance_metric, ")"), xlab = "", sub = "")
  })
  
  output$parallel_var_select <- renderUI({
    req(data_store$sample_metadata)
    selectizeInput("parallel_vars", "Select Metadata Variables:",
                   choices = names(data_store$sample_metadata),
                   selected = c("Age", "Gender", "BMI", "Ongoing_conditions"),
                   multiple = TRUE)
  })
  
  output$parallel_plot <- renderPlotly({
    req(data_store$sample_metadata, data_store$taxonomy_data, input$parallel_vars)
    
    # Calcular Total Abundance por muestra
    taxa_data <- data_store$taxonomy_data %>%
      group_by(Sample_ID) %>%
      summarize(Total_Abundance = sum(Abundance), .groups = "drop")
    
    # Unir con metadata
    combined <- merge(data_store$sample_metadata, taxa_data, by = "Sample_ID")
    
    # Variables seleccionadas
    vars <- input$parallel_vars
    
    # Añadir Total_Abundance al final
    all_vars <- c(vars, "Total_Abundance")
    
    # Seleccionar columnas y preparar datos
    plot_data <- combined %>%
      select(all_of(all_vars)) %>%
      na.omit() %>%
      mutate(across(where(is.character), as.factor)) %>%
      mutate(across(where(is.factor), ~ as.numeric(as.factor(.))))
    
    # Crear el gráfico
    plot_ly(
      type = 'parcoords',
      line = list(
        color = plot_data$Total_Abundance,
        colorscale = 'Viridis',
        showscale = TRUE  # 🔹 Oculta leyenda flotante, usamos solo el eje
      ),
      dimensions = lapply(names(plot_data), function(col) {
        if (is.character(combined[[col]]) || is.factor(combined[[col]])) {
          # Es cualitativa, convertir a factor y mostrar etiquetas
          f <- factor(combined[[col]])
          list(
            label = col,
            values = as.numeric(f),
            tickvals = seq_along(levels(f)),
            ticktext = levels(f)
          )
        } else {
          # Variable numérica normal
          list(
            label = col,
            values = plot_data[[col]]
          )
        }
      })
    )
  })
  
}





  
  
  

# Run app
shinyApp(ui, server)





