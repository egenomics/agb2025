
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
                    ),
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
                                         )),
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
      
      # TAB 2: Taxonomy Tab 
      tabItem(tabName = "taxonomy",
              fluidRow(
                box(
                  title = "Taxonomy Settings", 
                  status = "primary", 
                  solidHeader = TRUE,
                  width = 3,
                  collapsible = TRUE, 
                  
                  # Core taxonomy settings 
                  selectInput("tax_level", "Taxonomic Level:",
                              choices = c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"),
                              selected = "Phylum"),
                  
                  # Metadata filtering options 
                  selectInput("group_var", "Primary Group By:",
                              choices = c("None", "Gender", "Geographical_origin", "Age_group", "Institution", 
                                          "Ongoing_conditions", "Antibiotic_intake", "Smoking_status", 
                                          "Alcohol_consumption", "Exercise_frequency"),
                              selected = "None"),
                  
                  # Demographic comparison view 
                  checkboxInput("enable_comparison", "Enable Group Comparison", value = FALSE), 
                  conditionalPanel(
                    condition = "input.enable_comparison == true", 
                    selectInput("compare_var", "Compare By:", 
                                choices = c("Gender", "Geographical_origin", "Age_group", "Institution", 
                                            "Ongoing_conditions", "Antibiotic_intake", "Smoking_status", 
                                            "Alcohol_consumption", "Exercise_frequency"),
                                selected = "Gender"), 
                    selectInput("compare_value1", "Group 1:", choices = NULL), 
                    selectInput("compare_value2", "Group 2:", choices = NULL)
                  ),
                  
                  # Data transformation options
                  selectInput("normalization", "Data Normalization:",
                              choices = c("Relative Abundance (%)" = "percentage", 
                                          "Log10 Transformation" = "log10", 
                                          "Z-Score" = "zscore", 
                                          "Presence/Absence" = "binary"), 
                              selected = "percentage"), 
                  
                  # Common settings 
                  sliderInput("min_abundance", "Minimum Abundance (%):",
                              min = 0, max = 5, value = 1, step = 0.1), 
                  checkboxInput("sort_by_abundance", "Sort by Abundance", value = TRUE),
                  checkboxInput("show_legend", "Show Legend", value = TRUE), 
                  checkboxInput("use_plotly", "Use interactive plot", value = FALSE)
                ), 
                
                # Visualization settings 
                box(
                  title = "Visualization Options", 
                  status = "primary", 
                  solidHeader = TRUE, 
                  width = 3, 
                  collapsible = TRUE, 
                  
                  # Plot types 
                  radioButtons("plot_type", "Plot Type:", 
                               choices = c("Stacked Bar" = "bar", 
                                           "Heatmap" = "heatmap", 
                                           "Bubble Plot" = "bubble", 
                                           "Geographic Map" = "geo_map", 
                                           "Demographic Pattern" = "demo_pattern"), 
                               selected = "bar"), 
                  
                  # Conditional options based on plot type 
                  conditionalPanel(
                    condition = "input.plot_type == 'bar'", 
                    checkboxInput("show_error_bars", "Show Error Bars", value = FALSE), 
                    selectInput("bar_color_scheme", "Color Scheme:",
                                choices = c("Viridis" = "viridis", 
                                            "Set1" = "Set1", 
                                            "Set2" = "Set2",
                                            "Dark2" = "Dark2",
                                            "Paired" = "Paired"),
                                selected = "viridis")
                  ),
                  conditionalPanel(
                    condition = "input.plot_type == 'heatmap'",
                    selectInput("heatmap_palette", "Color Palette:",
                                choices = c("Viridis" = "viridis", 
                                            "Plasma" = "plasma", 
                                            "Magma" = "magma",
                                            "Inferno" = "inferno",
                                            "RdBu" = "RdBu",
                                            "RdYlBu" = "RdYlBu"),
                                selected = "viridis"),
                    selectInput("cluster_method", "Clustering Method:",
                                choices = c("None" = "none",
                                            "Euclidean" = "euclidean",
                                            "Bray-Curtis" = "bray"),
                                selected = "none")
                  ), 
                  conditionalPanel(
                    condition = "input.plot_type == 'geo_map'",
                    selectInput("geo_aggregation", "Geographic Aggregation:",
                                choices = c("Region" = "region", 
                                            "Country" = "country"),
                                selected = "region"),
                    selectInput("map_color_var", "Color By Taxa:",
                                choices = NULL)
                  ),
                  conditionalPanel(
                    condition = "input.plot_type == 'demo_pattern'",
                    selectInput("demo_variable", "Demographic Variable:",
                                choices = c("Age" = "Age", 
                                            "Gender" = "Gender",
                                            "BMI" = "BMI"),
                                selected = "Age"),
                    selectInput("pattern_type", "Pattern View:",
                                choices = c("Trend Analysis" = "trend", 
                                            "Distribution" = "distribution"),
                                selected = "trend"),
                    selectInput("selected_taxa", "Highlight Taxa:", choices = NULL)
                  ),
                  
                  hr(), 
                  downloadButton("download_tax_plot", "Download Plot:", class = "btn-primary")
                ), 
                
                # Main plot area with tabs for multiple views 
                box(
                  title = "Taxonomic Composition", 
                  status = "info", 
                  solidHeader = TRUE,
                  width = 6,
                  tabsetPanel(
                    id = "tax_plot_tabs",
                    tabPanel("Main View", 
                             plotOutput("taxonomy_plot", height = 600)),
                    tabPanel("Comparison View", 
                             conditionalPanel(
                               condition = "input.enable_comparison == true",
                               plotOutput("comparison_plot", height = 600)
                             ),
                             conditionalPanel(
                               condition = "input.enable_comparison == false",
                               div(style = "text-align: center; margin-top: 250px;",
                                   "Enable comparison in settings to view side-by-side comparisons")
                             ))
                  )
                )
              ), 
              fluidRow(
                # Enhanced taxa table with more metrics
                box(
                  title = "Top Taxa Analysis", 
                  status = "info", 
                  solidHeader = TRUE,
                  width = 12,
                  tabsetPanel(
                    id = "taxa_table_tabs",
                    tabPanel("Overall", dataTableOutput("taxa_table")),
                    tabPanel("By Demographics", 
                             fluidRow(
                               column(3,
                                      selectInput("demo_table_var", "Group By:",
                                                  choices = c("Gender", "Geographical_origin", "Age_group", 
                                                              "Ongoing_conditions", "Antibiotic_intake"),
                                                  selected = "Gender")
                               ),
                               column(9, dataTableOutput("taxa_demo_table"))
                             ))
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
      
      
      
      
      # TAB 4: Lifestyle Tab 
      
      tabItem(tabName = "lifestyle",
              fluidRow(
                
                # LEFT PANEL: Lifestyle Settings
                box(
                  title = "Lifestyle Settings",
                  status = "primary",
                  solidHeader = TRUE,
                  width = 3,
                  collapsible = TRUE, 
                  
                  selectInput("tax_level_life", "Taxonomic Level:",
                              choices = c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"),
                              selected = "Phylum"),
                  selectInput("normalization", "Data Normalization:",
                              choices = c("Relative Abundance (%)" = "percentage", 
                                          "Log10 Transformation" = "log10", 
                                          "Z-Score" = "zscore", 
                                          "Presence/Absence" = "binary"), 
                              selected = "percentage"),
                  sliderInput("min_abundance", "Minimum Abundance (%):", min = 0, max = 5, value = 1, step = 0.1),
                  
                  tags$hr(),
                  selectInput("lifestyle_variable", "Group By Lifestyle Variable:",
                              choices = c("Smoking_status", "Alcohol_consumption", 
                                          "Exercise_frequency", "Exercise_intensity", 
                                          "Bowel_movement_frequency", "Bowel_movement_quality", "Diet_type")),
                  checkboxInput("enable_comparison", "Enable Lifestyle Group Comparison", value = FALSE),
                  conditionalPanel(
                    condition = "input.enable_comparison == true", 
                    selectInput("compare_var", "Compare Lifestyle Variable:", 
                                choices = c("Smoking_status", "Alcohol_consumption", 
                                            "Exercise_frequency", "Exercise_intensity", 
                                            "Bowel_movement_frequency", "Bowel_movement_quality", "Diet_type"),
                                selected = "Smoking_status"),
                    selectInput("compare_value1", "Group 1:", choices = NULL), 
                    selectInput("compare_value2", "Group 2:", choices = NULL)
                  ),
                  
                  tags$hr(),
                  checkboxInput("sort_by_abundance", "Sort by Abundance", value = TRUE),
                  checkboxInput("show_legend", "Show Legend", value = TRUE),
                  checkboxInput("use_plotly", "Use interactive plot", value = FALSE),
                  
                  tags$hr(),
                  selectInput("diversity_metric_life", "Diversity Metric:",
                              choices = c("Shannon", "Invsimpson", "Simpson"),
                              selected = "Shannon"),
                  helpText("Used in Diversity View and Combined Score")
                ),
                
                # MIDDLE PANEL: Visualization Options
                box(
                  title = "Visualization Options", 
                  status = "primary", 
                  solidHeader = TRUE, 
                  width = 3, 
                  collapsible = TRUE, 
                  
                  # Plot type selector (rendered dynamically in server)
                  uiOutput("plot_type_life_ui"),
                  
                  # Conditional options
                  conditionalPanel(
                    condition = "input.plot_type_life == 'Stacked Bar'",
                    selectInput("bar_color_life", "Color Scheme:",
                                choices = c("Viridis", "Set1", "Dark2", "Pastel1", "Paired"),
                                selected = "Set1"),
                    checkboxInput("show_error_bars_life", "Show Error Bars", value = FALSE)
                  ),
                  
                  conditionalPanel(
                    condition = "input.plot_type_life == 'Boxplot' || input.plot_type_life == 'Boxplot + Points' || input.plot_type_life == 'Violin' || input.plot_type_life == 'Beeswarm'",
                    selectInput("group_color_life", "Color Palette:",
                                choices = c("Set1", "Dark2", "Pastel1", "Paired", "Viridis"),
                                selected = "Set1")
                  ),
                  
                  conditionalPanel(
                    condition = "input.plot_type_life == 'Trend'",
                    checkboxInput("smooth_trend_life", "Add Smoothed Line", value = TRUE),
                    selectInput("trend_palette_life", "Line Colors:",
                                choices = c("Set1", "Dark2", "Paired", "Pastel1", "Viridis"),
                                selected = "Dark2")
                  ),
                  
                  conditionalPanel(
                    condition = "input.plot_type_life == 'Heatmap'",
                    selectInput("heatmap_palette_life", "Heatmap Color Gradient:",
                                choices = c("YlGnBu", "YlOrRd", "Blues", "Greens", "RdPu", "Oranges", "PuBu", "BuPu"),
                                selected = "YlGnBu")
                  ),
                  
                  
                  hr(),
                  downloadButton("download_lifestyle_plot", "Download Plot", class = "btn-primary")
                ),
                
                # RIGHT PANEL: Output Panel
                box(
                  title = "Lifestyle Impact",
                  status = "info",
                  solidHeader = TRUE,
                  width = 6,
                  
                  tabsetPanel(
                    id = "active_lifestyle_tab",  ## <-- Needed to detect tab
                    tabPanel("Diversity View", 
                             plotOutput("lifestyle_diversity_plot", height = "600px"),
                             verbatimTextOutput("lifestyle_stat_test")),
                    
                    tabPanel("Taxon View", 
                             plotOutput("lifestyle_taxon_plot", height = "600px")),
                    
                    tabPanel("Combined Score", 
                             plotOutput("lifestyle_score_plot", height = "600px"))
                  )
                )
              )
      ),
      
      
      
      tabItem(tabName = "interventions", h3("Medical Interventions - Under Development")),
      tabItem(tabName = "multifactor", h3("Multi-factor Analysis - Under Development")),
      
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
    
    # smoking_counts <- table(data_store$sample_metadata$Smoking)
    # smoking <- paste(names(smoking_counts), paste0("(", smoking_counts, ")"), collapse = ", ")
    # 
    # alcohol_counts <- table(data_store$sample_metadata$Alcohol_Consumption)
    # alcohol <- paste(names(alcohol_counts), paste0("(", alcohol_counts, ")"), collapse = ", ")
    
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
  
  # Create reactive objects to store processed taxonomy data
  taxonomy_data <- reactive({
    req(data_store$taxonomy_data, data_store$sample_metadata)
    
    # Join taxonomy data with sample metadata
    taxonomy_joined <- merge(data_store$taxonomy_data, 
                             data_store$sample_metadata, 
                             by = "Sample_ID")
    
    # Return the processed data
    return(taxonomy_joined)
  })
  
  # Handle the dynamic UI updates for comparison values based on selected variable
  observe({
    req(input$enable_comparison, input$compare_var, data_store$sample_metadata)
    if(input$enable_comparison) {
      # Get unique values for the selected comparison variable
      compare_values <- unique(data_store$sample_metadata[[input$compare_var]])
      
      # Update the select inputs
      updateSelectInput(session, "compare_value1", 
                        choices = compare_values,
                        selected = compare_values[1])
      
      updateSelectInput(session, "compare_value2", 
                        choices = compare_values,
                        selected = if(length(compare_values) > 1) compare_values[2] else compare_values[1])
    }
  })
  
  
  observe({
    req(input$enable_div_comparison, input$grouping_var, data_store$sample_metadata)
    
    # Don't proceed if 'None' is selected
    if (input$grouping_var == "None") return()
    
    group_values <- unique(na.omit(data_store$sample_metadata[[input$grouping_var]]))
    
    updateSelectInput(session, "compare_group1", 
                      choices = group_values,
                      selected = group_values[1])
    
    updateSelectInput(session, "compare_group2", 
                      choices = group_values,
                      selected = if (length(group_values) > 1) group_values[2] else group_values[1])
  })
  
  
  # Update the selected taxa choices for demographic pattern plots
  observe({
    req(data_store$taxonomy_data)
    
    # Get taxonomic data at the selected level
    taxa_names <- unique(data_store$taxonomy_data$Species)
    
    # Update the select inputs
    updateSelectInput(session, "selected_taxa",
                      choices = taxa_names,
                      selected = taxa_names[1])
    
    # Also update map_color_var for geo_map
    updateSelectInput(session, "map_color_var",
                      choices = taxa_names,
                      selected = taxa_names[1])
  })
  
  # Function to process taxonomy data based on selected level and filters
  process_taxonomy_data <- reactive({
    req(taxonomy_data(), input$tax_level, input$normalization)
    
    tax_data <- taxonomy_data()
    
    # Default to Species level
    tax_level_col <- "Species"
    
    # Extract the selected taxonomic level (simulated)
    if(input$tax_level == "Phylum") {
      # For demonstration, extract first word as "phylum"
      tax_data$Phylum <- sapply(strsplit(as.character(tax_data$Species), " "), `[`, 1)
      tax_level_col <- "Phylum"
    } else if(input$tax_level == "Genus") {
      # For demonstration, extract first word as "genus"
      tax_data$Genus <- sapply(strsplit(as.character(tax_data$Species), " "), `[`, 1)
      tax_level_col <- "Genus"
    }
    
    # Apply normalization
    if(input$normalization == "percentage") {
      # Calculate relative abundance per sample
      tax_data <- tax_data %>%
        group_by(Sample_ID) %>%
        mutate(Abundance = (Abundance / sum(Abundance)) * 100) %>%
        ungroup()
    } else if(input$normalization == "log10") {
      # Log10 transformation (add small value to avoid log(0))
      tax_data$Abundance <- log10(tax_data$Abundance + 1)
    } else if(input$normalization == "zscore") {
      # Z-score normalization per sample
      tax_data <- tax_data %>%
        group_by(Sample_ID) %>%
        mutate(Abundance = (Abundance - mean(Abundance)) / max(.Machine$double.eps, sd(Abundance))) %>%
        ungroup()
    } else if(input$normalization == "binary") {
      # Presence/absence
      tax_data$Abundance <- ifelse(tax_data$Abundance > 0, 1, 0)
    }
    
    # Aggregate data by taxonomic level
    tax_data_agg <- tax_data %>%
      group_by(Sample_ID, !!sym(tax_level_col)) %>%
      summarize(Abundance = sum(Abundance),
                .groups = "drop")
    
    # Check if min_abundance is defined before applying filter
    if(!is.null(input$min_abundance) && input$min_abundance > 0 && input$normalization == "percentage") {
      # Calculate mean abundance per taxon
      taxon_means <- tax_data_agg %>%
        group_by(!!sym(tax_level_col)) %>%
        summarize(MeanAbundance = mean(Abundance), .groups = "drop")
      
      # Identify low abundance taxa
      low_abundance_taxa <- taxon_means %>%
        filter(MeanAbundance < input$min_abundance) %>%
        pull(!!sym(tax_level_col))
      
      # Replace low abundance taxa with "Other"
      if(length(low_abundance_taxa) > 0) {
        tax_data_agg[[tax_level_col]] <- ifelse(
          tax_data_agg[[tax_level_col]] %in% low_abundance_taxa,
          "Other",
          tax_data_agg[[tax_level_col]]
        )
        
        # Re-aggregate after replacing with "Other"
        tax_data_agg <- tax_data_agg %>%
          group_by(Sample_ID, !!sym(tax_level_col)) %>%
          summarize(Abundance = sum(Abundance), .groups = "drop")
      }
    }
    
    # Join with metadata for grouping
    tax_data_with_meta <- merge(tax_data_agg, 
                                data_store$sample_metadata, 
                                by = "Sample_ID")
    
    return(list(
      data = tax_data_with_meta,
      tax_level = tax_level_col
    ))
  })
  
  # Generate the main taxonomy plot UI
  output$taxonomy_plot <- renderUI({
    req(data_store$taxonomy_data)
    
    # Check if we should use plotly or ggplot
    if(input$plot_type == "bar" && input$use_plotly) {
      # Use plotly output
      plotlyOutput("plotly_bar", height = 600)
    } else {
      # Use regular plot output
      plotOutput("taxonomy_plot_static", height = 600)
    }
  })
  
  # Add static plot output
  output$taxonomy_plot_static <- renderPlot({
    req(process_taxonomy_data())
    
    # Get processed data
    processed <- process_taxonomy_data()
    tax_data <- processed$data
    tax_level_col <- processed$tax_level
    
    # Check if we have data to plot
    if(nrow(tax_data) == 0) {
      return(ggplot() + 
               geom_text(aes(x = 0.5, y = 0.5, label = "No data available.")) +
               theme_void())
    }
    
    # Initialize plot
    p <- NULL
    
    # Handle different plot types
    if(input$plot_type == "bar") {
      # Stacked bar plot
      if(!is.null(input$group_var) && input$group_var != "None") {
        # Group by selected variable
        if(input$sort_by_abundance) {
          # Sort samples by group and abundance
          sample_order <- tax_data %>%
            group_by(Sample_ID, !!sym(input$group_var)) %>%
            summarize(TotalAbundance = sum(Abundance), .groups = "drop") %>%
            arrange(!!sym(input$group_var), desc(TotalAbundance)) %>%
            pull(Sample_ID)
          
          tax_data$Sample_ID <- factor(tax_data$Sample_ID, levels = sample_order)
        }
        
        # Create grouped stacked bar plot
        p <- ggplot(tax_data, aes(x = Sample_ID, y = Abundance, fill = !!sym(tax_level_col))) +
          geom_bar(stat = "identity") +
          labs(
            title = paste("Taxonomic Composition at", input$tax_level, "Level"),
            x = "Sample",
            y = ifelse(input$normalization == "percentage", "Relative Abundance (%)", "Abundance"),
            fill = input$tax_level
          ) +
          theme_minimal() +
          theme(
            axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
            legend.position = ifelse(input$show_legend, "right", "none")
          ) +
          facet_grid(. ~ !!sym(input$group_var), scales = "free_x", space = "free_x")
        
      } else {
        # Non-grouped bar plot
        if(input$sort_by_abundance) {
          # Sort samples by abundance
          sample_order <- tax_data %>%
            group_by(Sample_ID) %>%
            summarize(TotalAbundance = sum(Abundance), .groups = "drop") %>%
            arrange(desc(TotalAbundance)) %>%
            pull(Sample_ID)
          
          tax_data$Sample_ID <- factor(tax_data$Sample_ID, levels = sample_order)
        }
        
        # Create stacked bar plot
        p <- ggplot(tax_data, aes(x = Sample_ID, y = Abundance, fill = !!sym(tax_level_col))) +
          geom_bar(stat = "identity") +
          labs(
            title = paste("Taxonomic Composition at", input$tax_level, "Level"),
            x = "Sample",
            y = ifelse(input$normalization == "percentage", "Relative Abundance (%)", "Abundance"),
            fill = input$tax_level
          ) +
          theme_minimal() +
          theme(
            axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
            legend.position = ifelse(input$show_legend, "right", "none")
          )
      }
      
      # Add error bars if requested and if grouping is enabled
      if(!is.null(input$show_error_bars) && input$show_error_bars && !is.null(input$group_var) && input$group_var != "None") {
        # Calculate mean and std error for each group
        error_data <- tax_data %>%
          group_by(!!sym(input$group_var), !!sym(tax_level_col)) %>%
          summarize(
            MeanAbundance = mean(Abundance),
            SE = sd(Abundance) / sqrt(n()),
            .groups = "drop"
          )
        
        # Create a new position variable for error bars
        group_counts <- length(unique(tax_data[[input$group_var]]))
        error_data$position <- as.numeric(factor(error_data[[input$group_var]]))
        
        # Add error bars to the plot
        p <- p + 
          geom_errorbar(
            data = error_data,
            aes(
              x = position,  
              y = MeanAbundance,
              ymin = MeanAbundance - SE,
              ymax = MeanAbundance + SE
            ),
            width = 0.2,
            position = position_dodge(0.9)
          )
      }
      
      # Apply color scheme
      if(input$bar_color_scheme == "viridis") {
        p <- p + scale_fill_viridis_d()
      } else {
        p <- p + scale_fill_brewer(palette = input$bar_color_scheme)
      }
    } else if(input$plot_type == "heatmap") {
      # For heatmap plot - simplified implementation here since ComplexHeatmap might not be available
      # Use ggplot2 for a basic heatmap instead
      
      if(!is.null(input$group_var) && input$group_var != "None") {
        # Group samples for heatmap
        heat_data <- tax_data %>%
          group_by(!!sym(tax_level_col), !!sym(input$group_var)) %>%
          summarize(MeanAbundance = mean(Abundance), .groups = "drop")
        
        # Create a basic heatmap using ggplot2
        p <- ggplot(heat_data, aes(x = !!sym(input$group_var), y = !!sym(tax_level_col), fill = MeanAbundance)) +
          geom_tile() +
          labs(
            title = paste("Heatmap of", input$tax_level, "by", input$group_var),
            x = input$group_var,
            y = input$tax_level,
            fill = ifelse(input$normalization == "percentage", "Mean Abundance (%)", "Mean Abundance")
          ) +
          theme_minimal() +
          theme(axis.text.x = element_text(angle = 45, hjust = 1))
      } else {
        # Sample-level heatmap
        # Limit to top taxa for readability
        taxa_abundance <- tax_data %>%
          group_by(!!sym(tax_level_col)) %>%
          summarize(MeanAbundance = mean(Abundance), .groups = "drop") %>%
          arrange(desc(MeanAbundance)) %>%
          head(20)
        
        top_taxa <- taxa_abundance %>% pull(!!sym(tax_level_col))
        
        # Filter for top taxa
        filtered_data <- tax_data %>%
          filter(!!sym(tax_level_col) %in% top_taxa)
        
        p <- ggplot(filtered_data, aes(x = Sample_ID, y = !!sym(tax_level_col), fill = Abundance)) +
          geom_tile() +
          labs(
            title = paste("Heatmap of Top", input$tax_level),
            x = "Sample",
            y = input$tax_level,
            fill = ifelse(input$normalization == "percentage", "Abundance (%)", "Abundance")
          ) +
          theme_minimal() +
          theme(
            axis.text.x = element_text(angle = 90, hjust = 1, size = 7)
          )
      }
      
      # Apply color palette
      if(input$heatmap_palette == "viridis") {
        p <- p + scale_fill_viridis_c()
      } else if(input$heatmap_palette %in% c("plasma", "magma", "inferno")) {
        p <- p + scale_fill_viridis_c(option = input$heatmap_palette)
      } else {
        p <- p + scale_fill_distiller(palette = input$heatmap_palette)
      }
      
    } else if(input$plot_type == "bubble") {
      # Bubble plot implementation
      if(!is.null(input$group_var) && input$group_var != "None") {
        # Group samples
        bubble_data <- tax_data %>%
          group_by(!!sym(tax_level_col), !!sym(input$group_var)) %>%
          summarize(MeanAbundance = mean(Abundance), .groups = "drop")
        
        p <- ggplot(bubble_data, aes(x = !!sym(input$group_var), y = !!sym(tax_level_col), 
                                     size = MeanAbundance, color = !!sym(tax_level_col))) +
          geom_point(alpha = 0.7) +
          scale_size_continuous(range = c(1, 10)) +
          labs(
            title = paste("Bubble Plot of", input$tax_level, "by", input$group_var),
            x = input$group_var,
            y = input$tax_level,
            size = ifelse(input$normalization == "percentage", "Mean Abundance (%)", "Abundance")
          ) +
          theme_minimal() +
          guides(color = guide_legend(title = input$tax_level))
        
      } else {
        # Sample-level bubble plot (limited to top taxa)
        taxa_abundance <- tax_data %>%
          group_by(!!sym(tax_level_col)) %>%
          summarize(MeanAbundance = mean(Abundance), .groups = "drop") %>%
          arrange(desc(MeanAbundance)) %>%
          head(15)
        
        top_taxa <- taxa_abundance %>% pull(!!sym(tax_level_col))
        
        filtered_data <- tax_data %>%
          filter(!!sym(tax_level_col) %in% top_taxa)
        
        p <- ggplot(filtered_data, aes(x = Sample_ID, y = !!sym(tax_level_col), 
                                       size = Abundance, color = !!sym(tax_level_col))) +
          geom_point(alpha = 0.7) +
          scale_size_continuous(range = c(1, 8)) +
          labs(
            title = paste("Bubble Plot of Top", input$tax_level),
            x = "Sample",
            y = input$tax_level,
            size = ifelse(input$normalization == "percentage", "Abundance (%)", "Abundance")
          ) +
          theme_minimal() +
          theme(
            axis.text.x = element_text(angle = 90, hjust = 1, size = 7),
            legend.position = ifelse(input$show_legend, "right", "none")
          )
      }
      
    } else if(input$plot_type == "geo_map") {
      # Geographic map implementation (simplified)
      
      # Verify map_color_var exists
      if(is.null(input$map_color_var) || !input$map_color_var %in% unique(tax_data[[tax_level_col]])) {
        return(ggplot() + 
                 geom_text(aes(0.5, 0.5, label = "Please select a valid taxon for mapping.")) +
                 theme_void())
      }
      
      # Aggregate data by geographic region and selected taxon
      geo_data <- tax_data %>%
        filter(!!sym(tax_level_col) == input$map_color_var) %>%
        group_by(Geographical_origin) %>%
        summarize(MeanAbundance = mean(Abundance), .groups = "drop")
      
      # Check if we have data after filtering
      if(nrow(geo_data) == 0) {
        return(ggplot() + 
                 geom_text(aes(0.5, 0.5, label = "No geographic data available for selected taxon.")) +
                 theme_void())
      }
      
      # This is a simplified placeholder for a geographic map
      p <- ggplot(geo_data, aes(x = Geographical_origin, y = MeanAbundance, fill = Geographical_origin)) +
        geom_bar(stat = "identity") +
        labs(
          title = paste("Geographic Distribution of", input$map_color_var),
          x = "Geographic Region",
          y = ifelse(input$normalization == "percentage", "Mean Abundance (%)", "Mean Abundance"),
          fill = "Region"
        ) +
        theme_minimal() +
        theme(axis.text.x = element_text(angle = 45, hjust = 1))
      
    } else if(input$plot_type == "demo_pattern") {
      # Demographic pattern implementation
      
      # Make sure selected_taxa exists
      if(is.null(input$selected_taxa) || !input$selected_taxa %in% unique(tax_data[[tax_level_col]])) {
        return(ggplot() + 
                 geom_text(aes(0.5, 0.5, label = "Please select a valid taxon.")) +
                 theme_void())
      }
      
      # Filter for selected taxon
      demo_data <- tax_data %>%
        filter(!!sym(tax_level_col) == input$selected_taxa)
      
      # Check if we have any data after filtering
      if(nrow(demo_data) == 0) {
        return(ggplot() + 
                 geom_text(aes(0.5, 0.5, label = "No data available for selected taxon.")) +
                 theme_void())
      }
      
      if(!is.null(input$pattern_type) && input$pattern_type == "trend") {
        # Trend analysis by demographic variable
        if(!is.null(input$demo_variable) && input$demo_variable == "Age") {
          # For age, we can do a scatter plot with trend line
          p <- ggplot(demo_data, aes(x = Age, y = Abundance)) +
            geom_point(alpha = 0.6) +
            geom_smooth(method = "loess", se = TRUE) +
            labs(
              title = paste("Abundance of", input$selected_taxa, "by Age"),
              x = "Age",
              y = ifelse(input$normalization == "percentage", "Abundance (%)", "Abundance")
            ) +
            theme_minimal()
          
        } else if(!is.null(input$demo_variable) && input$demo_variable == "BMI") {
          # Similar approach for BMI
          p <- ggplot(demo_data, aes(x = BMI, y = Abundance)) +
            geom_point(alpha = 0.6) +
            geom_smooth(method = "loess", se = TRUE) +
            labs(
              title = paste("Abundance of", input$selected_taxa, "by BMI"),
              x = "BMI",
              y = ifelse(input$normalization == "percentage", "Abundance (%)", "Abundance")
            ) +
            theme_minimal()
          
        } else {
          # For categorical variables like Gender
          p <- ggplot(demo_data, aes(x = !!sym(input$demo_variable), y = Abundance, fill = !!sym(input$demo_variable))) +
            geom_boxplot() +
            labs(
              title = paste("Abundance of", input$selected_taxa, "by", input$demo_variable),
              x = input$demo_variable,
              y = ifelse(input$normalization == "percentage", "Abundance (%)", "Abundance")
            ) +
            theme_minimal()
        }
        
      } else if(!is.null(input$pattern_type) && input$pattern_type == "distribution") {
        # Distribution analysis
        if(!is.null(input$demo_variable) && input$demo_variable %in% c("Age", "BMI")) {
          # For continuous variables, create density plots
          p <- ggplot(demo_data, aes(x = Abundance, fill = cut_width(!!sym(input$demo_variable), width = 5))) +
            geom_density(alpha = 0.7) +
            labs(
              title = paste("Distribution of", input$selected_taxa, "by", input$demo_variable, "Groups"),
              x = ifelse(input$normalization == "percentage", "Abundance (%)", "Abundance"),
              y = "Density",
              fill = paste(input$demo_variable, "Group")
            ) +
            theme_minimal()
          
        } else {
          # For categorical variables
          p <- ggplot(demo_data, aes(x = Abundance, fill = !!sym(input$demo_variable))) +
            geom_density(alpha = 0.7) +
            labs(
              title = paste("Distribution of", input$selected_taxa, "by", input$demo_variable),
              x = ifelse(input$normalization == "percentage", "Abundance (%)", "Abundance"),
              y = "Density",
              fill = input$demo_variable
            ) +
            theme_minimal()
        }
      }
    }
    
    return(p)
  })
  
  # Plotly output
  output$plotly_bar <- renderPlotly({
    req(process_taxonomy_data())
    
    # Get processed data
    processed <- process_taxonomy_data()
    tax_data <- processed$data
    tax_level_col <- processed$tax_level
    
    # Check if we have data
    if(nrow(tax_data) == 0) {
      # Return an empty plotly plot with a message
      return(plot_ly() %>% 
               add_trace(type = "scatter", mode = "text", text = "No data available") %>%
               layout(title = "No Data"))
    }
    
    # Handle different grouping scenarios
    if(!is.null(input$group_var) && input$group_var != "None") {
      # Group by selected variable
      if(input$sort_by_abundance) {
        # Sort samples by group and abundance
        sample_order <- tax_data %>%
          group_by(Sample_ID, !!sym(input$group_var)) %>%
          summarize(TotalAbundance = sum(Abundance), .groups = "drop") %>%
          arrange(!!sym(input$group_var), desc(TotalAbundance)) %>%
          pull(Sample_ID)
        
        tax_data$Sample_ID <- factor(tax_data$Sample_ID, levels = sample_order)
      }
      
      # Create grouped stacked bar plot with plotly
      p <- plot_ly(tax_data, x = ~Sample_ID, y = ~Abundance, color = ~get(tax_level_col),
                   type = "bar", text = ~paste(get(tax_level_col), ": ", round(Abundance, 2)),
                   hoverinfo = "text") %>%
        layout(
          title = paste("Taxonomic Composition at", input$tax_level, "Level"),
          xaxis = list(title = "Sample", tickangle = 45),
          yaxis = list(title = ifelse(input$normalization == "percentage", "Relative Abundance (%)", "Abundance")),
          barmode = "stack",
          showlegend = input$show_legend,
          margin = list(b = 100)  # Add margin to bottom for rotated labels
        )
      
      # Add faceting by group
      # Since plotly doesn't have direct faceting, we create annotations for groups
      group_levels <- unique(tax_data[[input$group_var]])
      samples_by_group <- split(tax_data$Sample_ID, tax_data[[input$group_var]])
      
      annotations <- list()
      for(i in seq_along(group_levels)) {
        grp <- group_levels[i]
        samples <- unique(samples_by_group[[as.character(grp)]])
        
        # Check if samples exist for this group
        if(length(samples) > 0) {
          # Find the level indices that match our samples
          level_indices <- which(levels(tax_data$Sample_ID) %in% samples)
          
          # Only proceed if we found matching indices
          if(length(level_indices) > 0) {
            midpoint <- mean(level_indices)
            
            annotations[[i]] <- list(
              x = midpoint,
              y = 1.05, 
              text = as.character(grp),
              xref = "x",
              yref = "paper",
              showarrow = FALSE,
              font = list(size = 14)
            )
          }
        }
      }
      
      # Only add annotations if we have some
      if(length(annotations) > 0) {
        p <- p %>% layout(annotations = annotations)
      }
      
    } else {
      # Non-grouped bar plot
      if(input$sort_by_abundance) {
        # Sort samples by abundance
        sample_order <- tax_data %>%
          group_by(Sample_ID) %>%
          summarize(TotalAbundance = sum(Abundance), .groups = "drop") %>%
          arrange(desc(TotalAbundance)) %>%
          pull(Sample_ID)
        
        tax_data$Sample_ID <- factor(tax_data$Sample_ID, levels = sample_order)
      }
      
      # Create stacked bar plot with plotly
      p <- plot_ly(tax_data, x = ~Sample_ID, y = ~Abundance, color = ~get(tax_level_col),
                   type = "bar", text = ~paste(get(tax_level_col), ": ", round(Abundance, 2)),
                   hoverinfo = "text") %>%
        layout(
          title = paste("Taxonomic Composition at", input$tax_level, "Level"),
          xaxis = list(title = "Sample", tickangle = 45),
          yaxis = list(title = ifelse(input$normalization == "percentage", "Relative Abundance (%)", "Abundance")),
          barmode = "stack",
          showlegend = input$show_legend,
          margin = list(b = 100)  # Add margin to bottom for rotated labels
        )
    }
    
    # Apply color scheme through plotly
    # Count unique taxa safely
    n_taxa <- length(unique(tax_data[[tax_level_col]]))
    
    if(input$bar_color_scheme == "viridis") {
      colors <- viridis::viridis(max(n_taxa, 1))
    } else {
      # Use a safe approach to get colors
      if(n_taxa <= 9) {
        colors <- RColorBrewer::brewer.pal(max(n_taxa, 3), input$bar_color_scheme)
        if(n_taxa < 3) colors <- colors[1:n_taxa]
      } else {
        base_colors <- RColorBrewer::brewer.pal(9, input$bar_color_scheme)
        colors <- colorRampPalette(base_colors)(n_taxa)
      }
    }
    
    # Apply colors to plotly
    p <- p %>% layout(colorway = colors)
    
    return(p)
  })

  ##################
  
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
 
  
  
  
  
  
  
  ##################
  
  
  ###### LIFESTYLE TAB ######
  
  process_lifestyle_taxonomy_data <- reactive({
    req(data_store$taxonomy_data, data_store$sample_metadata)
    
    tax_data <- data_store$taxonomy_data
    tax_level_col <- "Species"  # default
    
    # Simulate other levels from Species name
    if (input$tax_level_life == "Phylum") {
      tax_data$Phylum <- sapply(strsplit(as.character(tax_data$Species), " "), `[`, 1)
      tax_level_col <- "Phylum"
    } else if (input$tax_level_life == "Genus") {
      tax_data$Genus <- sapply(strsplit(as.character(tax_data$Species), " "), `[`, 1)
      tax_level_col <- "Genus"
    }
    
    # Normalization
    if (input$normalization == "percentage") {
      tax_data <- tax_data %>%
        group_by(Sample_ID) %>%
        mutate(Abundance = Abundance / sum(Abundance) * 100) %>%
        ungroup()
    } else if (input$normalization == "log10") {
      tax_data$Abundance <- log10(tax_data$Abundance + 1)
    } else if (input$normalization == "zscore") {
      tax_data <- tax_data %>%
        group_by(Sample_ID) %>%
        mutate(Abundance = scale(Abundance)[, 1]) %>%
        ungroup()
    } else if (input$normalization == "binary") {
      tax_data$Abundance <- ifelse(tax_data$Abundance > 0, 1, 0)
    }
    
    # Aggregate by taxonomic level
    tax_data_agg <- tax_data %>%
      group_by(Sample_ID, !!sym(tax_level_col)) %>%
      summarise(Abundance = sum(Abundance), .groups = "drop")
    
    # Join with metadata
    tax_data_merged <- left_join(tax_data_agg, data_store$sample_metadata, by = "Sample_ID")
    
    return(list(data = tax_data_merged, tax_col = tax_level_col))
  })
  
  output$plot_type_life_ui <- renderUI({
    req(input$active_lifestyle_tab)
    
    choices <- if (input$active_lifestyle_tab == "Diversity View") {
      c("Boxplot + Points", "Violin", "Beeswarm", "Trend")
    } else if (input$active_lifestyle_tab == "Taxon View") {
      choices = c("Stacked Bar", "Boxplot + Points", "Violin", "Beeswarm", "Trend", "Heatmap")
    } else {
      character(0)
    }
    
    radioButtons("plot_type_life", "Plot Type:", choices = choices, selected = choices[1])
  })
  
  
  output$lifestyle_diversity_plot <- renderPlot({
    req(data_store$taxonomy_data, data_store$sample_metadata, 
        input$lifestyle_variable, input$diversity_metric_life, 
        input$plot_type_life)
    
    tax_data <- data_store$taxonomy_data
    meta_data <- data_store$sample_metadata
    
    # Reshape to wide format
    wide_data <- tax_data %>%
      select(Sample_ID, Species, Abundance) %>%
      pivot_wider(names_from = Species, values_from = Abundance, values_fill = 0)
    
    # Compute diversity
    diversity_metrics <- wide_data %>%
      column_to_rownames("Sample_ID") %>%
      as.matrix() %>%
      {
        tibble(
          Sample_ID = rownames(.),
          Shannon = vegan::diversity(., index = "shannon"),
          Simpson = vegan::diversity(., index = "simpson"),
          Invsimpson = vegan::diversity(., index = "invsimpson")
        )
      }
    
    # Merge with metadata
    diversity_df <- left_join(diversity_metrics, meta_data, by = "Sample_ID")
    
    metric_col <- input$diversity_metric_life
    group_var <- input$lifestyle_variable
    palette <- input$group_color_life
    
    # Start plot
    p <- ggplot(diversity_df, aes(x = .data[[group_var]], y = .data[[metric_col]], fill = .data[[group_var]]))
    
    # Add selected plot type
    if (input$plot_type_life == "Boxplot + Points") {
      p <- p + 
        geom_boxplot(alpha = 0.6, outlier.shape = NA) +
        geom_jitter(width = 0.2, shape = 21, alpha = 0.6, color = "black")
      
    } else if (input$plot_type_life == "Violin") {
      p <- p + geom_violin(alpha = 0.7, trim = FALSE)
      
    } else if (input$plot_type_life == "Beeswarm") {
      if (!requireNamespace("ggbeeswarm", quietly = TRUE)) {
        showNotification("Install 'ggbeeswarm' package to use beeswarm plots.", type = "error")
        return(NULL)
      }
      p <- p + ggbeeswarm::geom_beeswarm(priority = "density", cex = 1.5, size = 1.5, shape = 21, alpha = 0.6)
      
    } else if (input$plot_type_life == "Trend") {
      add_smooth <- isTRUE(input$smooth_trend_life)
      
      # Handle palette color safely
      trend_color <- if (input$trend_palette_life == "Viridis") {
        viridisLite::viridis(1)
      } else {
        RColorBrewer::brewer.pal(8, input$trend_palette_life)[1]
      }
      
      p <- ggplot(diversity_df, aes(x = as.numeric(factor(.data[[group_var]])), 
                                    y = .data[[metric_col]], group = 1)) +
        geom_point(alpha = 0.5)
      
      if (add_smooth) {
        p <- p + geom_smooth(method = "loess", se = TRUE, color = trend_color)
      }
      
      p <- p +
        scale_x_continuous(breaks = 1:length(unique(diversity_df[[group_var]])),
                           labels = unique(diversity_df[[group_var]])) +
        labs(x = group_var)
    }
    
    
    
    # Labels & Theme
    p <- p +
      labs(
        title = paste(metric_col, "Diversity by", group_var),
        y = paste(metric_col, "Index")
      ) +
      theme_minimal() +
      theme(legend.position = ifelse(input$show_legend, "right", "none"))
    
    # Color palette (except Trend)
    if (!(input$plot_type_life %in% c("Trend")) && !is.null(palette)) {
      p <- p + scale_fill_brewer(palette = palette)
    }
    
    return(p)
  })
  
  
  
  output$lifestyle_stat_test <- renderPrint({
    req(data_store$taxonomy_data, data_store$sample_metadata, input$lifestyle_variable)
    
    tax_data <- data_store$taxonomy_data
    meta_data <- data_store$sample_metadata
    
    # Reshape to wide
    wide_data <- tax_data %>%
      select(Sample_ID, Species, Abundance) %>%
      pivot_wider(names_from = Species, values_from = Abundance, values_fill = 0)
    
    # Compute diversity
    diversity_values <- vegan::diversity(as.matrix(wide_data %>% column_to_rownames("Sample_ID")),
                                         index = tolower(input$diversity_metric_life))
    
    # Merge with metadata
    diversity_df <- tibble(
      Sample_ID = wide_data$Sample_ID,
      Diversity = diversity_values
    ) %>%
      left_join(meta_data, by = "Sample_ID")
    
    # Perform Kruskal-Wallis test
    formula_str <- as.formula(paste("Diversity ~ `", input$lifestyle_variable, "`", sep = ""))
    kruskal.test(formula_str, data = diversity_df)
  })
  
  
  output$lifestyle_taxon_plot <- renderPlot({
    req(input$lifestyle_variable)
    
    processed <- process_lifestyle_taxonomy_data()
    tax_data <- processed$data
    tax_col <- processed$tax_col
    
    # Normalize abundance
    tax_data <- tax_data %>%
      group_by(Sample_ID) %>%
      mutate(NormAbundance = case_when(
        input$normalization == "percentage" ~ Abundance,
        input$normalization == "log10" ~ log10(Abundance + 1),
        input$normalization == "zscore" ~ scale(Abundance)[, 1],
        input$normalization == "binary" ~ as.numeric(Abundance > 0),
        TRUE ~ Abundance
      )) %>%
      ungroup()
    
    # Filter by minimum abundance
    if (input$normalization == "percentage") {
      tax_data <- tax_data %>% filter(NormAbundance >= input$min_abundance)
    }
    
    # Select top N taxa
    top_taxa <- tax_data %>%
      group_by(!!sym(tax_col)) %>%
      summarise(mean_abund = mean(NormAbundance), .groups = "drop") %>%
      arrange(desc(mean_abund)) %>%
      slice_head(n = 5) %>%
      pull(!!sym(tax_col))
    
    filtered_data <- tax_data %>% filter(!!sym(tax_col) %in% top_taxa)
    
    # Plot logic
    if (input$plot_type_life == "Violin") {
      plt <- ggplot(filtered_data, aes(x = .data[[input$lifestyle_variable]], y = NormAbundance, fill = .data[[tax_col]])) +
        geom_violin(alpha = 0.7, trim = FALSE, position = position_dodge()) +
        scale_fill_brewer(palette = input$group_color_life) +
        facet_wrap(as.formula(paste("~", tax_col)), scales = "fixed", ncol = 3)
      
    } else if (input$plot_type_life == "Boxplot + Points") {
      plt <- ggplot(filtered_data, aes(x = .data[[input$lifestyle_variable]], y = NormAbundance, fill = .data[[tax_col]])) +
        geom_boxplot(alpha = 0.6, position = position_dodge(), outlier.shape = NA) +
        geom_jitter(width = 0.2, alpha = 0.5, shape = 21, color = "black") +
        scale_fill_brewer(palette = input$group_color_life) +
        facet_wrap(as.formula(paste("~", tax_col)), scales = "free_y")
      
    } else if (input$plot_type_life == "Beeswarm") {
      if (!requireNamespace("ggbeeswarm", quietly = TRUE)) {
        showNotification("Install 'ggbeeswarm' to use Beeswarm plots.", type = "error")
        return(NULL)
      }
      plt <- ggplot(filtered_data, aes(x = .data[[input$lifestyle_variable]], y = NormAbundance, fill = .data[[tax_col]])) +
        ggbeeswarm::geom_beeswarm(priority = "density", cex = 1.5, size = 1.5, shape = 21, alpha = 0.6) +
        scale_fill_brewer(palette = input$group_color_life) +
        facet_wrap(as.formula(paste("~", tax_col)), scales = "free_y")
      
    } else if (input$plot_type_life == "Trend") {
      plt <- ggplot(filtered_data, aes(x = .data[[input$lifestyle_variable]],
                                       y = NormAbundance,
                                       group = .data[[tax_col]],
                                       color = .data[[tax_col]])) +
        geom_point(alpha = 0.5)
      
      if (input$smooth_trend_life) {
        plt <- plt + geom_smooth(method = "loess", se = FALSE)
      } else {
        plt <- plt + geom_line(linewidth = 1)
      }
      
      plt <- plt +
        scale_color_brewer(palette = input$trend_palette_life) +
        facet_wrap(as.formula(paste("~", tax_col)), scales = "free_y") +
        labs(x = input$lifestyle_variable)
      
    } else if (input$plot_type_life == "Heatmap") {
      # Prepare data for heatmap: average per group x taxon
      heatmap_data <- filtered_data %>%
        group_by(.data[[input$lifestyle_variable]], .data[[tax_col]]) %>%
        summarise(mean_abund = mean(NormAbundance), .groups = "drop") %>%
        rename(Group = 1, Taxon = 2)
      
      plt <- ggplot(heatmap_data, aes(x = Group, y = Taxon, fill = mean_abund)) +
        geom_tile(color = "white") +
        scale_fill_distiller(palette = input$heatmap_palette_life, direction = 1) +
        theme_minimal() +
        labs(x = input$lifestyle_variable, y = tax_col, fill = "Mean Abundance") +
        theme(axis.text.x = element_text(angle = 45, hjust = 1))
      
    } else {
      # Default: Stacked Bar Plot
      plt <- ggplot(filtered_data, aes_string(x = input$lifestyle_variable, y = "NormAbundance", fill = tax_col)) +
        geom_bar(stat = "summary", fun = mean, position = "stack") +
        scale_fill_brewer(palette = input$bar_color_life)
      if (input$show_error_bars_life) {
        plt <- plt + stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.3)
      }
    }
    
    # Final formatting
    plt +
      labs(
        title = paste("Top", input$tax_level_life, "by", input$lifestyle_variable),
        y = paste("Abundance (", input$normalization, ")")
      ) +
      theme(
        panel.background = element_rect(fill = "white"),
        plot.background = element_rect(fill = "white"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
        strip.text = element_text(face = "bold", size = 12),
        strip.background = element_blank(),
        strip.placement = "outside",
        axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
        legend.position = ifelse(input$show_legend, "right", "none"),
        panel.spacing = unit(5, "lines"),
        plot.margin = margin(10, 10, 10, 10)
      )
  })
  
  
  
  
  output$lifestyle_score_plot <- renderPlot({
    req(data_store$taxonomy_data, data_store$sample_metadata, input$lifestyle_variable, input$diversity_metric_life)
    
    merged_data <- merge(data_store$taxonomy_data, data_store$sample_metadata, by = "Sample_ID")
    
    diversity_scores <- merged_data %>%
      group_by(Sample_ID) %>%
      summarize(Diversity = vegan::diversity(Abundance, index = tolower(input$diversity_metric_life)))
    
    faecali_data <- merged_data %>%
      filter(grepl("Faecalibacterium", Species)) %>%
      group_by(Sample_ID) %>%
      summarize(Faecali_Abund = sum(Abundance), .groups = "drop")
    
    score_df <- merge(diversity_scores, faecali_data, by = "Sample_ID", all = TRUE)
    score_df$Faecali_Abund[is.na(score_df$Faecali_Abund)] <- 0
    
    score_df <- merge(score_df, data_store$sample_metadata[, c("Sample_ID", input$lifestyle_variable)], by = "Sample_ID")
    colnames(score_df)[ncol(score_df)] <- "LifestyleGroup"
    
    score_df$Score <- scale(score_df$Diversity) + scale(score_df$Faecali_Abund)
    
    ggplot(score_df, aes(x = LifestyleGroup, y = Score, fill = LifestyleGroup)) +
      geom_boxplot(alpha = 0.7) +
      theme_minimal() +
      labs(title = paste("Combined Score (", input$diversity_metric_life, " + Faecalibacterium)", sep = ""),
           x = input$lifestyle_variable,
           y = "Combined Score") +
      theme(legend.position = ifelse(input$show_legend, "right", "none"))
  })
  
  
  
  

  
  
  
  # # Generate comparison plot
  # output$comparison_plot <- renderPlot({
  #   req(process_taxonomy_data(), input$enable_comparison, input$compare_var, 
  #       input$compare_value1, input$compare_value2)
  #   
  #   # Get processed data
  #   processed <- process_taxonomy_data()
  #   tax_data <- processed$data
  #   tax_level_col <- processed$tax_level
  #   
  #   # Filter data for the two comparison groups
  #   comp_data <- tax_data %>%
  #     filter(!!sym(input$compare_var) %in% c(input$compare_value1, input$compare_value2))
  #   
  #   # Calculate mean abundance by taxon and group
  #   comp_summary <- comp_data %>%
  #     group_by(!!sym(tax_level_col), !!sym(input$compare_var)) %>%
  #     summarize(
  #       MeanAbundance = mean(Abundance),
  #       SE = sd(Abundance) / sqrt(n()),
  #       .groups = "drop"
  #     )
  #   
  #   # Filter to top taxa for readability
  #   if(nrow(comp_summary) > 0) {
  #     # Get top taxa based on overall abundance
  #     top_taxa <- comp_data %>%
  #       group_by(!!sym(tax_level_col)) %>%
  #       summarize(MeanAbundance = mean(Abundance), .groups = "drop") %>%
  #       arrange(desc(MeanAbundance)) %>%
  #       head(15) %>%
  #       pull(!!sym(tax_level_col))
  #     
  #     # Filter for top taxa
  #     comp_summary <- comp_summary %>%
  #       filter(!!sym(tax_level_col) %in% top_taxa)
  #     
  #     # Create comparison plot
  #     p <- ggplot(comp_summary, aes(x = !!sym(tax_level_col), y = MeanAbundance, 
  #                                   fill = !!sym(input$compare_var))) +
  #       geom_bar(stat = "identity", position = position_dodge()) +
  #       geom_errorbar(aes(ymin = MeanAbundance - SE, ymax = MeanAbundance + SE),
  #                     position = position_dodge(0.9), width = 0.25) +
  #       labs(
  #         title = paste("Comparison of", input$tax_level, "between", input$compare_value1, "and", input$compare_value2),
  #         x = input$tax_level,
  #         y = ifelse(input$normalization == "percentage", "Mean Abundance (%)", "Mean Abundance"),
  #         fill = input$compare_var
  #       ) +
  #       theme_minimal() +
  #       theme(axis.text.x = element_text(angle = 45, hjust = 1))
  #     
  #     return(p)
  #   }
  #   
  #   # Return empty plot if no data
  #   ggplot() + 
  #     theme_void() + 
  #     geom_text(aes(0, 0, label = "No data available for comparison"))
  # })
  # 
  # # Generate taxa table by demographics
  # output$taxa_demo_table <- renderDataTable({
  #   req(process_taxonomy_data(), input$demo_table_var)
  #   
  #   # Get processed data
  #   processed <- process_taxonomy_data()
  #   tax_data <- processed$data
  #   tax_level_col <- processed$tax_level
  #   
  #   # Summarize by taxonomic level and demographic variable
  #   taxa_demo_summary <- tax_data %>%
  #     group_by(!!sym(tax_level_col), !!sym(input$demo_table_var)) %>%
  #     summarize(
  #       MeanAbundance = mean(Abundance),
  #       MedianAbundance = median(Abundance),
  #       StdDev = sd(Abundance),
  #       Prevalence = sum(Abundance > 0) / n() * 100,
  #       SampleCount = n(),
  #       .groups = "drop"
  #     ) %>%
  #     arrange(!!sym(input$demo_table_var), desc(MeanAbundance))
  #   
  #   # Render the data table
  #   datatable(taxa_demo_summary,
  #             options = list(pageLength = 10, scrollX = TRUE),
  #             rownames = FALSE) %>%
  #     formatRound(columns = c("MeanAbundance", "MedianAbundance", "StdDev", "Prevalence"), digits = 2)
  # })
  # 
  # # Generate overall taxa table
  # output$taxa_table <- renderDataTable({
  #   req(process_taxonomy_data())
  #   
  #   # Get processed data
  #   processed <- process_taxonomy_data()
  #   tax_data <- processed$data
  #   tax_level_col <- processed$tax_level
  #   
  #   # Summarize by taxonomic level
  #   taxa_summary <- tax_data %>%
  #     group_by(!!sym(tax_level_col)) %>%
  #     summarize(
  #       MeanAbundance = mean(Abundance),
  #       MedianAbundance = median(Abundance),
  #       MaxAbundance = max(Abundance),
  #       StdDev = sd(Abundance),
  #       Prevalence = sum(Abundance > 0) / n() * 100,
  #       SampleCount = n(),
  #       .groups = "drop"
  #     ) %>%
  #     arrange(desc(MeanAbundance))
  #   
  #   # Render the data table
  #   datatable(taxa_summary,
  #             options = list(pageLength = 10, scrollX = TRUE),
  #             rownames = FALSE) %>%
  #     formatRound(columns = c("MeanAbundance", "MedianAbundance", "MaxAbundance", "StdDev", "Prevalence"), digits = 2)
  # })
  # 
  # # Download handler for taxonomy plot
  # output$download_tax_plot <- downloadHandler(
  #   filename = function() {
  #     paste("taxonomy_plot_", input$tax_level, "_", Sys.Date(), ".png", sep = "")
  #   },
  #   content = function(file) {
  #     # Directly render the plot again instead of trying to access output objects
  #     if(input$tax_plot_tabs == "Main View") {
  #       # Get processed data again
  #       processed <- process_taxonomy_data()
  #       tax_data <- processed$data
  #       tax_level_col <- processed$tax_level
  #       
  #       # Recreate the plot logic from taxonomy_plot_static
  #       # Handle different plot types
  #       if(input$plot_type == "bar") {
  #         # Stacked bar plot logic (simplified)
  #         p <- ggplot(tax_data, aes(x = Sample_ID, y = Abundance, fill = !!sym(tax_level_col))) +
  #           geom_bar(stat = "identity") +
  #           labs(title = paste("Taxonomic Composition at", input$tax_level, "Level"))
  #       } else if(input$plot_type == "heatmap") {
  #         # Heatmap logic (simplified)
  #         p <- ggplot(tax_data, aes(x = Sample_ID, y = !!sym(tax_level_col), fill = Abundance)) +
  #           geom_tile() +
  #           labs(title = paste("Heatmap of", input$tax_level))
  #       } else {
  #         # Default plot if other types not handled
  #         p <- ggplot(tax_data, aes(x = Sample_ID, y = Abundance)) +
  #           geom_bar(stat = "identity") +
  #           labs(title = "Taxonomy Plot")
  #       }
  #     } else {
  #       # For comparison plot
  #       processed <- process_taxonomy_data()
  #       tax_data <- processed$data
  #       tax_level_col <- processed$tax_level
  #       
  #       # Filter data for the two comparison groups
  #       comp_data <- tax_data %>%
  #         filter(!!sym(input$compare_var) %in% c(input$compare_value1, input$compare_value2))
  #       
  #       # Check if we have data after filtering
  #       if(nrow(comp_data) > 0) {
  #         # Calculate mean abundance
  #         comp_summary <- comp_data %>%
  #           group_by(!!sym(tax_level_col), !!sym(input$compare_var)) %>%
  #           summarize(
  #             MeanAbundance = mean(Abundance),
  #             SE = sd(Abundance) / sqrt(max(1, n())), # Avoid division by zero
  #             .groups = "drop"
  #           )
  #         
  #         # Check if we have data after summarizing
  #         if(nrow(comp_summary) > 0) {
  #           # Create comparison plot
  #           p <- ggplot(comp_summary, aes(x = !!sym(tax_level_col), y = MeanAbundance, 
  #                                         fill = !!sym(input$compare_var))) +
  #             geom_bar(stat = "identity", position = position_dodge()) +
  #             labs(title = paste("Comparison Plot"))
  #         } else {
  #           # No data after summarizing
  #           p <- ggplot() + 
  #             annotate("text", x = 0, y = 0, label = "No data available for comparison") +
  #             theme_void()
  #         }
  #       } else {
  #         # No data after filtering
  #         p <- ggplot() + 
  #           annotate("text", x = 0, y = 0, label = "No data available for comparison") +
  #           theme_void()
  #       }
  #     }
  #     
  #     # Save the plot safely
  #     tryCatch({
  #       ggsave(file, plot = p, width = 10, height = 8, dpi = 300)
  #     }, error = function(e) {
  #       # Create a simple error message plot if ggsave fails
  #       png(file, width = 800, height = 600)
  #       plot(1, 1, type = "n", axes = FALSE, xlab = "", ylab = "")
  #       text(1, 1, "Error generating plot: No data available")
  #       dev.off()
  #     })
  
  
  # 
  # # Taxonomy plot 
  # output$taxonomy_plot <- renderPlot({
  #   ps <- phyloseq_obj()
  #   if (is.null(ps)) {
  #     return(NULL)
  #   }
  #   
  #   # Get selected taxonomic level 
  #   tax_level <- input$tax_level
  #   
  #   # Agglomerate at the selected taxonomic level 
  #   ps_glom <- tax_glom(ps, taxrank = tax_level)
  #   
  #   # Transform to relavtive abundance
  #   ps_rel <- transform_sample_counts(ps_glom, function(x) x / sum(x) * 100)
  #   
  #   # Melt to long format for ggplot 
  #   ps_melt <- psmelt(ps_rel)
  #   
  #   # Filter low abundance taxa 
  #   ps_melt <- ps_melt %>%
  #     group_by(get(tax_level)) %>%
  #     mutate(MeanAbundance = mean(Abundance)) %>%
  #     ungroup()
  #   
  #   low_abundance_taxa <- ps_melt %>%
  #     filter(MeanAbundance < input$min_abundance) %>%
  #     pull(!!sym(tax_level)) %>%
  #     unique()
  #   
  #   if (length(low_abundance_taxa) > 0) {
  #     ps_melt <- ps_melt %>%
  #       mutate(!!tax_level := ifelse(get(tax_level) %in% low_abundance_taxa, "Other", get(tax_level)))
  #   }
  #   
  #   # Create plot based on plot type
  #   if (input$plot_type == "bar") {
  #     # Group by the selected variable if not "None"
  #     if (input$group_var != "None") {
  #       group_var <- input$group_var
  #       
  #       # Reorder samples by group and abundance if requested 
  #       if (input$sort_by_abundance) {
  #         sample_order <- ps_melt %>%
  #           group_by(Sample, !!sym(group_var)) %>%
  #           summarize(TotalAbundance = sum(Abundance), .groups = "drop") %>%
  #           arrange(!!sym(group_var), desc(TotalAbundance)) %>%
  #           pull(Sample)
  #         ps_melt$Sample <- factor(ps_melt$Sample, levels = sample_order)
  #       }
  #       
  #       # Create the grouped stacked bar plot 
  #       p <- ggplot(ps_melt, aes(x = Sample, y = Abundance, fill = get(tax_level))) +
  #         geom_bar(stat = "identity") +
  #         labs(x = "Sample", y = "Relative Abundance (%)", fill = "tax_level") +
  #         theme_minimal() +
  #         theme(
  #           axis.text.x = element_text(angle = 45, hjust = 1, size = 8), 
  #           legend.position = if (input$show_legend) "right" else "none"
  #         ) +
  #         facet_grid(. ~ get(group_var), scales = "free_x", space = "free_x")
  #     } else {
  #       # Create the non-grouped stacked bar plot. Reorder samples by abundance if requested
  #       if (input$sort_by_abundance) {
  #         sample_order <- ps_melt %>%
  #           group_by(Sample) %>%
  #           summarize(TotalAbundance = sum(Abundance), .groups = "drop") %>%
  #           arrange(desc(TotalAbundance)) %>%
  #           pull(Sample)
  #         ps_melt$Sample <- factor(ps_melt$Sample, levels = sample_order)
  #       }
  #       
  #       p <- ggplot(ps_melt, aes(x = Sample, y = Abundance, fill = get(tax_level))) +
  #         geom_bar(stat = "identity") +
  #         labs(x = "Sample", y = "Relative Abundance (%)", fill = tax_level) +
  #         theme_minimal() +
  #         theme(
  #           axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
  #           legend.position = if(input$show_legend) "right" else "none"
  #         )
  #     }
  #     
  #     # Use a colorful palette 
  #     if (length(unique(ps_melt[[tax_level]])) <= 9) {
  #       p <- p + scale_fill_brewer(palette = "Set1")
  #     } else {
  #       p <- p + scale_fill_viridis_d()
  #     }
  #   } else if (input$plot_type == "heatmap") {
  #     # Create heatmap for taxonomic abundance. Summarize data for heatmap
  #     if (input$group_var != "None") {
  #       # Group samples for heatmap
  #       heat_data <- ps_melt %>%
  #         group_by(get(tax_level), !!sym(input$group_var)) %>%
  #         summarize(MeanAbundance = mean(Abundance), .groups = "drop") %>%
  #         spread(key = !!sym(input$group_var), value = MeanAbundance, fill = 0)
  #       
  #       # Convert back to long format for ggplot 
  #       heat_data_long <- gather(heat_data, key = "Group", value = "Abundance", -1)
  #       colnames(heat_data_long)[1] <- tax_level
  #       
  #       # Order taxa by overall abundance 
  #       taxa_order <- ps_melt %>%
  #         group_by(!!sym(tax_level)) %>%
  #         summarize(MeanAbundance = mean(Abundance), .groups = "drop") %>%
  #         arrange(desc(MeanAbundance)) %>%
  #         pull(!!sym(tax_level))
  #       heat_data_long[[tax_level]] <- factor(heat_data_long[[tax_level]], levels = taxa_order)
  #       
  #       # Create heatmap
  #       p <- ggplot(heat_data_long, aes(x = Group, y =!!sym(tax_level), fill = Abundance)) +
  #         geom_tile() +
  #         scale_fill_viridis_c(option = "plasma") +
  #         labs(x = input$group_var, y = tax_level, fill = "Mean\nAbundance (%)") +
  #         theme_minimal() +
  #         theme(
  #           axis.text.x = element_text(angle = 45, hjust = 1),
  #           legend.position = "right"
  #         )
  #     } else {
  #       # Sample-level heatmap. Order taxa and samples
  #       taxa_order <- ps_melt %>%
  #         group_by(!!sym(tax_level)) %>%
  #         summarize(MeanAbundance = mean(Abundance), .groups = "drop") %>%
  #         arrange(desc(MeanAbundance)) %>%
  #         head(20) %>% # Limit to top taxa for readability 
  #         pull(!!sym(tax_level))
  #       
  #       # Filter for top taxa 
  #       ps_melt_filtered <- ps_melt %>%
  #         filter(!!sym(tax_level) %in% taxa_order)
  #       
  #       ps_melt_filtered[[tax_level]] <- factor(ps_melt_filtered[[tax_level]], levels = taxa_order)
  #       
  #       # Create heatmap
  #       p <- ggplot(ps_melt_filtered, aes(x = Sample, y = !!sym(tax_level), fill = Abundance)) +
  #         geom_tile() +
  #         scale_fill_viridis_c(option = "plasma") +
  #         labs(x = "Sample", y = tax_level, fill = "Abundance (%)") +
  #         theme_minimal() +
  #         theme(
  #           axis.text.x = element_text(angle = 90, hjust = 1, size = 7),
  #           legend.position = "right"
  #         )
  #     }
  #   }
  #   return(p)
  # })
  # 
  # # Top taxa table
  # output$taxa_table <- renderDataTable({
  #   ps <- phyloseq_obj()
  #   if (is.null(ps)) {
  #     return(NULL)
  #   }
  #   # Get selected taxonomic level 
  #   tax_level <- input$tax_level
  #   
  #   # Agglomerate and calculate relative abundance 
  #   ps_glom <- tax_glom(ps, taxrank = tax_level)
  #   ps_rel <- transform_sample_counts(ps_glom, function(x) x / sum(x) * 100)
  #   ps_melt <- psmelt(ps_rel)
  #   
  #   # Summarize by selected taxonomic level 
  #   taxa_summary <- ps_melt %>%
  #     group_by(!!sym(tax_level)) %>%
  #     summarize(
  #       MeanRelativeAbundance = mean(Abundance),
  #       MaxRelativeAbundance = max(Abundance), 
  #       Prevalence = sum(Abundance > 0) / length(unique(Sample)) * 100, 
  #       .groups = "drop"
  #     ) %>%
  #     arrange(desc(MeanRelativeAbundance))
  #   
  #   datatable(taxa_summary, 
  #             options = list(pageLength = 10, scrollX = TRUE),
  #             rownames = FALSE) %>%
  #     formatRound(columns = c("MeanRelativeAbundance", "MaxRelativeAbundance", "Prevalence"), digits = 2)
  # })
  
  ###################################################################
  
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





