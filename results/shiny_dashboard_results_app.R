
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
library(ggtree)
library(treeio)
library(DESeq2)
library(broom)
library(shinyBS)
library(beyonce)
library(tayloRswift)

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
      menuItem(" Phylogenetic tree", tabName = "phylo", icon = icon("tree")), 
      menuItem(" Reports", tabName = "reports", icon = icon("file-export"))
    )
  ),
  
  # Body with dashboard style 
  dashboardBody(
    # Use BS theme for consistent styling 
    tags$head(
      tags$link(rel = "stylesheet", type = "text/css", 
                href = "https://cdnjs.cloudflare.com/ajax/libs/font-awesome/5.15.4/css/all.min.css"), 
      tags$style(HTML("
    /* Existing styles */
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
    
    /* NEW: Summary cards styling for shinydashboard */
    .summary-card {
      background: #fff;
      border: 1px solid #d2d6de;
      border-radius: 3px;
      margin-bottom: 15px;
      box-shadow: 0 1px 3px rgba(0,0,0,0.12);
      transition: box-shadow 0.3s ease;
    }
    
    .summary-card:hover {
      box-shadow: 0 2px 8px rgba(0,0,0,0.15);
    }
    
    .summary-card-header {
      background: #f4f4f4;
      border-bottom: 1px solid #d2d6de;
      padding: 10px 15px;
      font-weight: 600;
      border-radius: 3px 3px 0 0;
    }
    
    .summary-card-header.primary { background: #3c8dbc; color: #fff; }
    .summary-card-header.info { background: #00c0ef; color: #fff; }
    .summary-card-header.success { background: #00a65a; color: #fff; }
    .summary-card-header.warning { background: #f39c12; color: #fff; }
    .summary-card-header.secondary { background: #6c757d; color: #fff; }
    
    .summary-card-body {
      padding: 15px;
    }
    
    .stat-badge {
      display: inline-block;
      padding: 4px 8px;
      margin: 2px;
      background: #f4f4f4;
      border: 1px solid #ddd;
      border-radius: 3px;
      font-size: 12px;
      font-weight: 500;
    }
    
    .stat-badge.primary { background: #3c8dbc; color: #fff; border-color: #3c8dbc; }
    .stat-badge.info { background: #00c0ef; color: #fff; border-color: #00c0ef; }
    .stat-badge.success { background: #00a65a; color: #fff; border-color: #00a65a; }
    .stat-badge.warning { background: #f39c12; color: #fff; border-color: #f39c12; }
    .stat-badge.secondary { background: #6c757d; color: #fff; border-color: #6c757d; }
    
    /* Sample info card */
    .sample-info-card {
      background: #fff;
      border: 1px solid #d2d6de;
      border-radius: 3px;
      margin-top: 10px;
    }
    
    .sample-info-header {
      background: #3c8dbc;
      color: #fff;
      padding: 10px 15px;
      border-radius: 3px 3px 0 0;
      font-weight: 600;
    }
    
    .sample-info-table {
      margin: 0;
    }
    
    .sample-info-table th {
      background: #f9f9f9;
      border-top: none;
      font-weight: 600;
      width: 40%;
      padding: 8px 15px;
    }
    
    .sample-info-table td {
      border-top: 1px solid #f4f4f4;
      padding: 8px 15px;
      width: 60%;
    }
    
    .alert-custom {
      padding: 15px;
      margin-bottom: 20px;
      border: 1px solid transparent;
      border-radius: 4px;
    }
    
    .alert-info {
      color: #31708f;
      background-color: #d9edf7;
      border-color: #bce8f1;
    }
    
    /* Responsive grid */
    .summary-grid {
      display: flex;
      flex-wrap: wrap;
      margin: -7.5px;
    }
    
    .summary-col {
      flex: 0 0 33.333333%;
      max-width: 33.333333%;
      padding: 7.5px;
    }
    
    @media (max-width: 768px) {
      .summary-col {
        flex: 0 0 50%;
        max-width: 50%;
      }
    }
    
    @media (max-width: 576px) {
      .summary-col {
        flex: 0 0 100%;
        max-width: 100%;
      }
    }
  "))
    ), 
    
    tabItems(
      # Home tab 
      tabItem(tabName = "home", 
              # Hero Section
              fluidRow(
                box(
                  title = NULL,
                  status = "primary", 
                  solidHeader = FALSE, 
                  width = 12,
                  style = "background: linear-gradient(135deg, #667eea 0%, #764ba2 100%); color: white; border: none; box-shadow: 0 4px 15px rgba(0,0,0,0.1);",
                  div(
                    style = "text-align: center; padding: 40px 20px;",
                    div(
                      style = "font-size: 3.5em; margin-bottom: 15px;",
                      "🧬"
                    ),
                    h1("AGB Microbiome Analysis Platform", 
                       style = "font-size: 2.8em; font-weight: 300; margin-bottom: 15px; text-shadow: 2px 2px 4px rgba(0,0,0,0.3);"),
                    h3("Comprehensive Analysis Dashboard for Microbial Communities", 
                       style = "font-weight: 300; opacity: 0.9; margin-bottom: 30px;")
                  )
                )
              ),
              
              # Quick Stats Row
              fluidRow(
                valueBoxOutput("home_samples_count", width = 3),
                valueBoxOutput("home_taxa_count", width = 3), 
                valueBoxOutput("home_analyses_count", width = 3),
                valueBoxOutput("home_reports_count", width = 3)
              ),
              
              # Main Features Section
              fluidRow(
                # Data Management
                column(width = 6,
                       box(
                         title = NULL,
                         status = "info",
                         solidHeader = FALSE,
                         width = 12,
                         style = "border-top: 4px solid #3498db; box-shadow: 0 2px 10px rgba(0,0,0,0.08);",
                         div(
                           style = "padding: 20px;",
                           div(
                             style = "display: flex; align-items: center; margin-bottom: 15px;",
                             div(style = "font-size: 2.5em; margin-right: 15px; color: #3498db;", "📊"),
                             h3("Data Management", style = "margin: 0; color: #2c3e50;")
                           ),
                           p("Upload and manage your microbiome datasets with ease", 
                             style = "color: #7f8c8d; font-size: 1.1em; margin-bottom: 20px;"),
                           div(
                             style = "background: #f8f9fa; padding: 15px; border-radius: 8px; margin-bottom: 15px;",
                             h5("✓ Manual Upload", style = "color: #27ae60; margin-bottom: 8px;"),
                             p("Upload sample metadata, run data, and taxonomy files individually", 
                               style = "margin: 0; font-size: 0.95em; color: #6c757d;")
                           ),
                           div(
                             style = "background: #f8f9fa; padding: 15px; border-radius: 8px; margin-bottom: 20px;",
                             h5("✓ Automatic Import", style = "color: #27ae60; margin-bottom: 8px;"),
                             p("Connect to hospital data structures for seamless data loading", 
                               style = "margin: 0; font-size: 0.95em; color: #6c757d;")
                           ),
                         )
                       )
                ),
                
                # Analysis Tools
                column(width = 6,
                       box(
                         title = NULL,
                         status = "success",
                         solidHeader = FALSE,
                         width = 12,
                         style = "border-top: 4px solid #27ae60; box-shadow: 0 2px 10px rgba(0,0,0,0.08);",
                         div(
                           style = "padding: 20px;",
                           div(
                             style = "display: flex; align-items: center; margin-bottom: 15px;",
                             div(style = "font-size: 2.5em; margin-right: 15px; color: #27ae60;", "🔬"),
                             h3("Analysis Suite", style = "margin: 0; color: #2c3e50;")
                           ),
                           p("Comprehensive tools for microbiome data exploration", 
                             style = "color: #7f8c8d; font-size: 1.1em; margin-bottom: 20px;"),
                           div(
                             style = "background: #f8f9fa; padding: 15px; border-radius: 8px; margin-bottom: 15px;",
                             h5("✓ Taxonomic Profiling", style = "color: #27ae60; margin-bottom: 8px;"),
                             p("Visualize microbial composition with interactive plots and heatmaps", 
                               style = "margin: 0; font-size: 0.95em; color: #6c757d;")
                           ),
                           div(
                             style = "background: #f8f9fa; padding: 15px; border-radius: 8px; margin-bottom: 20px;",
                             h5("✓ Diversity Analysis", style = "color: #27ae60; margin-bottom: 8px;"),
                             p("Alpha and beta diversity metrics with statistical testing", 
                               style = "margin: 0; font-size: 0.95em; color: #6c757d;")
                           )
                         )
                       )
                )
              ),
              
              # Analysis Workflow Section
              fluidRow(
                box(
                  title = "Analysis Workflow",
                  status = "warning",
                  solidHeader = TRUE,
                  width = 12,
                  style = "box-shadow: 0 2px 10px rgba(0,0,0,0.08);",
                  div(
                    style = "padding: 20px;",
                    h4("Follow these steps for comprehensive microbiome analysis:", 
                       style = "color: #2c3e50; margin-bottom: 30px; text-align: center;"),
                    
                    # Workflow Steps
                    fluidRow(
                      # Step 1
                      column(width = 2,
                             div(
                               style = "text-align: center; padding: 20px; background: #ecf0f1; border-radius: 10px; margin: 10px;",
                               div(style = "font-size: 3em; color: #3498db; margin-bottom: 10px;", "1️⃣"),
                               h5("Data Upload", style = "color: #2c3e50; margin-bottom: 10px;"),
                               p("Load your microbiome datasets", style = "font-size: 0.9em; color: #7f8c8d; margin: 0;")
                             )
                      ),
                      
                      # Arrow
                      column(width = 1,
                             div(style = "text-align: center; padding-top: 50px; font-size: 1.5em; color: #bdc3c7;", "→")
                      ),
                      
                      # Step 2
                      column(width = 2,
                             div(
                               style = "text-align: center; padding: 20px; background: #ecf0f1; border-radius: 10px; margin: 10px;",
                               div(style = "font-size: 3em; color: #27ae60; margin-bottom: 10px;", "2️⃣"),
                               h5("Taxonomy", style = "color: #2c3e50; margin-bottom: 10px;"),
                               p("Explore microbial composition", style = "font-size: 0.9em; color: #7f8c8d; margin: 0;")
                             )
                      ),
                      
                      # Arrow
                      column(width = 1,
                             div(style = "text-align: center; padding-top: 50px; font-size: 1.5em; color: #bdc3c7;", "→")
                      ),
                      
                      # Step 3
                      column(width = 2,
                             div(
                               style = "text-align: center; padding: 20px; background: #ecf0f1; border-radius: 10px; margin: 10px;",
                               div(style = "font-size: 3em; color: #e74c3c; margin-bottom: 10px;", "3️⃣"),
                               h5("Diversity", style = "color: #2c3e50; margin-bottom: 10px;"),
                               p("Analyze alpha & beta diversity", style = "font-size: 0.9em; color: #7f8c8d; margin: 0;")
                             )
                      ),
                      
                      # Arrow
                      column(width = 1,
                             div(style = "text-align: center; padding-top: 50px; font-size: 1.5em; color: #bdc3c7;", "→")
                      ),
                      
                      # Step 4
                      column(width = 2,
                             div(
                               style = "text-align: center; padding: 20px; background: #ecf0f1; border-radius: 10px; margin: 10px;",
                               div(style = "font-size: 3em; color: #9b59b6; margin-bottom: 10px;", "4️⃣"),
                               h5("Reports", style = "color: #2c3e50; margin-bottom: 10px;"),
                               p("Generate comprehensive reports", style = "font-size: 0.9em; color: #7f8c8d; margin: 0;")
                             )
                      )
                    )
                  )
                )
              ),
              
              # Feature Highlights
              fluidRow(
                # Taxonomy Features
                column(width = 4,
                       box(
                         title = "🦠 Taxonomic Analysis",
                         status = "primary",
                         solidHeader = TRUE,
                         width = 12,
                         style = "box-shadow: 0 2px 10px rgba(0,0,0,0.08);",
                         tags$ul(
                           style = "padding-left: 20px;",
                           tags$li("Interactive relative abundance plots"),
                           tags$li("Taxonomic heatmaps with clustering"),
                           tags$li("Multiple statistical testing options:"),
                           tags$ul(
                             style = "margin-top: 5px; margin-bottom: 10px;",
                             tags$li("PERMANOVA for community structure"),
                             tags$li("DESeq2 for differential abundance"),
                             tags$li("Wilcoxon & t-tests for comparisons")
                           ),
                           tags$li("Customizable color schemes"),
                           tags$li("Patient vs Control comparisons")
                         )
                       )
                ),
                
                # Diversity Features
                column(width = 4,
                       box(
                         title = "📊 Diversity Metrics",
                         status = "success",
                         solidHeader = TRUE,
                         width = 12,
                         style = "box-shadow: 0 2px 10px rgba(0,0,0,0.08);",
                         div(
                           h5("Alpha Diversity:", style = "color: #27ae60; margin-bottom: 8px;"),
                           tags$ul(
                             style = "padding-left: 15px; margin-bottom: 15px;",
                             tags$li("Shannon, Simpson indices"),
                             tags$li("Observed OTUs/Species"),
                             tags$li("Box plots and violin plots")
                           ),
                           h5("Beta Diversity:", style = "color: #27ae60; margin-bottom: 8px;"),
                           tags$ul(
                             style = "padding-left: 15px;",
                             tags$li("Multiple distance metrics"),
                             tags$li("PCoA, NMDS, t-SNE, UMAP"),
                             tags$li("Clustering dendrograms"),
                             tags$li("Group ellipses & statistics")
                           )
                         )
                       )
                ),
                
                # Additional Features
                column(width = 4,
                       box(
                         title = "🔧 Advanced Tools",
                         status = "warning",
                         solidHeader = TRUE,
                         width = 12,
                         style = "box-shadow: 0 2px 10px rgba(0,0,0,0.08);",
                         tags$ul(
                           style = "padding-left: 20px;",
                           tags$li(strong("Phylogenetic Trees:"), "Visualize evolutionary relationships"),
                           tags$li(strong("Multifactor Analysis:"), "Parallel coordinate plots"),
                           tags$li(strong("Metadata Integration:"), "Clinical condition grouping"),
                           tags$li(strong("Control Samples:"), "Automatic comparison analysis"),
                           tags$li(strong("Export Options:"), "High-quality plots and data"),
                           tags$li(strong("Report Generation:"), "PDF, HTML, and Word formats")
                         )
                       )
                )
              ),
              
              # Quick Tips Section
              fluidRow(
                box(
                  title = "💡 Quick Tips for Best Results",
                  status = "info",
                  solidHeader = TRUE,
                  width = 12,
                  style = "box-shadow: 0 2px 10px rgba(0,0,0,0.08);",
                  fluidRow(
                    column(width = 6,
                           div(
                             style = "background: #f8f9fa; padding: 20px; border-radius: 8px; margin: 10px;",
                             h5("📁 Data Preparation", style = "color: #2c3e50; margin-bottom: 15px;"),
                             tags$ul(
                               tags$li("Ensure your CSV/TSV files have proper headers"),
                               tags$li("Include both patient and control samples when possible"),
                               tags$li("Check that sample IDs match across all files"),
                               tags$li("Remove any special characters from column names")
                             )
                           )
                    ),
                    column(width = 6,
                           div(
                             style = "background: #f8f9fa; padding: 20px; border-radius: 8px; margin: 10px;",
                             h5("🎯 Analysis Strategy", style = "color: #2c3e50; margin-bottom: 15px;"),
                             tags$ul(
                               tags$li("Start with taxonomy overview to understand your data"),
                               tags$li("Use appropriate statistical tests for your sample size"),
                               tags$li("Consider clinical metadata for meaningful groupings"),
                               tags$li("Generate reports to document your findings")
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
                # Manual Upload Section 
                box(
                  title = "Upload Microbiome Data Manually", 
                  status = "primary", 
                  solidHeader = TRUE, 
                  width = 12, 
                  collapsible = TRUE, 
                  collapsed = FALSE, 
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
                  )
                )
              ), 
              
              # Manual Statistics and Sample Selection (Only for Manual Upload)
              conditionalPanel(
                condition = "output.manual_data_loaded == true", 
                
                # Manual Data Status and Statistics 
                fluidRow(
                  column(
                    width = 12, 
                    box(
                      title = "Manual Data Loading Status", 
                      status = "info", 
                      solidHeader = TRUE, 
                      width = 12, 
                      verbatimTextOutput("manual_data_loading_status")
                    )
                  )
                ), 
                # Sample Selection Section
                fluidRow(
                  box(
                    title = "Sample Selection", 
                    status = "warning",
                    solidHeader = TRUE, 
                    width = 8, 
                    fluidRow(
                      column(
                        width = 4, 
                        selectizeInput("selected_sample_id", "Select Sample to Analyze:", 
                                       choices = NULL, 
                                       options = list(
                                         placeholder = "Choose a sample...", 
                                         onInitialize = I('function() { this.setValue(""); }')
                                       )), 
                        br(), 
                        actionButton("clear_sample_selection", "Clear Selection", 
                                     icon = icon("times"),
                                     class = "btn-warning")
                      ), 
                      column(
                        width = 8, 
                        htmlOutput("selected_sample_info")
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
                )
              ), 
              
              # Automatic Upload Section
              fluidRow(
                box(
                  title = "Automatic Data Upload", 
                  status = "success",
                  solidHeader = TRUE, 
                  width = 12, 
                  collapsible = TRUE, 
                  fluidRow(
                    column(
                      width = 12, 
                      h4("Hospital Data Structure"), 
                      p("Automatically load data from the hospital's data structure folder."), 
                      br(),
                      textInput("data_folder_path", "Data Folder Path:", 
                                value = "path/to/hospital/data/structure", 
                                placeholder = "Enter the path to your data structure folder"), 
                      br(), 
                      actionButton("scan_folder", "Scan Data Folder", 
                                   icon = icon("folder-open"), 
                                   class = "btn-warning"), 
                      br(), 
                      conditionalPanel(
                        condition = "output.folder_scanned", 
                        h5("Available Run IDs:"), 
                        selectizeInput("available_run_ids", "Select Run ID to Analyze:", 
                                       choices = NULL, 
                                       options = list(
                                         placeholder = "First scan the folder...", 
                                         onInitialize = I('function() { this.setValue(""); }')
                                       )),
                        actionButton("load_automatic_data", "Load Selected Run Data", 
                                     icon = icon("download"), 
                                     class = "btn-success")
                      ), 
                      br(), 
                      h5("Control Data Status:"), 
                      verbatimTextOutput("control_data_status")
                    )
                  )
                )
              ), 
              # Automatic data status and statistics 
              conditionalPanel(
                condition = "output.automatic_data_loaded == true", 
                fluidRow(
                  column(
                    width = 8, 
                    box(
                      title = "Automatic Data Loading Status", 
                      status = "info", 
                      solidHeader = TRUE, 
                      width = 12, 
                      verbatimTextOutput("automatic_data_loading_status")
                    )
                  ), 
                  column(
                    width = 4, 
                    valueBoxOutput("automatic_total_samples_box", width = 12), 
                    valueBoxOutput("automatic_total_species_box", width = 12), 
                    valueBoxOutput("automatic_total_runs_box", width = 12), 
                    valueBoxOutput("automatic_total_control_samples_box", width = 12)
                  )
                )
              ),
              
              # Metadata Previews
              conditionalPanel(
                condition = "output.manual_data_loaded == true || output.automatic_data_loaded == true", 
                fluidRow(
                  box(
                    title = "Metadata & Data Preview", 
                    status = "primary", 
                    solidHeader = TRUE, 
                    width = 12, 
                    tabBox(
                      width = 12, 
                      tabPanel("Sample Metadata", 
                               fluidRow(
                                 column(width = 12, 
                                        radioButtons("sample_metadata_filter", "Filter by Conditiions:", 
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
                )
              ), 
              # Data Overview 
              conditionalPanel(
                condition = "output.manual_data_loaded == true || output.automatic_data_loaded == true", 
                fluidRow(
                  box(
                    title = "Data Summary", 
                    status = "primary", 
                    solidHeader = TRUE, 
                    width = 12, 
                    uiOutput("data_summary_dynamic")
                  )
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
                                     selected = "Genus"), 
                         
                         # Filtering options 
                         radioButtons("sample_group", "Sample Group:", 
                                      choices = c("All Samples", "Patients Only", "Controls Only"), 
                                      selected = "All Samples"),
                         
                         hr(), 
                         
                         # Metadata grouping options 
                         selectInput("metadata_group", "Group by Metadata:", 
                                     choices = c("None", "Age", "Gender", "Body_Mass_Index", "Ongoing_conditions", "Allergies", "Dietary_Information", "Antibiotic_intake", "Exercise_frequency", "Alcohol_consumption"),
                                     selected = "None"), 
                         
                         checkboxInput("compare_groups", "Compare Patients vs Controls", value = TRUE), 
                         
                         hr(), 
                         
                         # Statistical testing options - SIMPLIFIED
                         selectInput("stat_test", "Statistical Test:", 
                                     choices = c("None", "PERMANOVA", "DESeq2", "Wilcoxon", "t-test"), 
                                     selected = "None")
                       ),
                       
                       # Statistical Tests Information Box
                       box(
                         title = "Statistical Tests Guide", 
                         status = "warning", 
                         solidHeader = TRUE, 
                         width = 12,
                         collapsible = TRUE,
                         collapsed = TRUE,
                         
                         div(style = "padding: 10px;",
                             HTML("
                         <div style='margin-bottom: 15px;'>
                           <h5 style='color: #2c3e50; margin-bottom: 8px;'><strong>PERMANOVA</strong></h5>
                           <p style='margin-bottom: 5px; font-size: 13px;'>
                             <strong>Purpose:</strong> Tests for differences in overall microbial community composition between groups.
                           </p>
                           <p style='margin-bottom: 0; font-size: 13px; color: #7f8c8d;'>
                             <strong>Best for:</strong> Comparing entire microbiome profiles, beta-diversity analysis.
                           </p>
                         </div>
                         
                         <div style='margin-bottom: 15px;'>
                           <h5 style='color: #2c3e50; margin-bottom: 8px;'><strong>DESeq2</strong></h5>
                           <p style='margin-bottom: 5px; font-size: 13px;'>
                             <strong>Purpose:</strong> Identifies individual taxa that are significantly different between groups.
                           </p>
                           <p style='margin-bottom: 0; font-size: 13px; color: #7f8c8d;'>
                             <strong>Best for:</strong> Finding specific bacteria/taxa that differ between conditions.
                           </p>
                         </div>
                         
                         <div style='margin-bottom: 15px;'>
                           <h5 style='color: #2c3e50; margin-bottom: 8px;'><strong>Wilcoxon Test</strong></h5>
                           <p style='margin-bottom: 5px; font-size: 13px;'>
                             <strong>Purpose:</strong> Non-parametric test for comparing abundance distributions between two groups.
                           </p>
                           <p style='margin-bottom: 0; font-size: 13px; color: #7f8c8d;'>
                             <strong>Best for:</strong> Robust comparisons when data is not normally distributed.
                           </p>
                         </div>
                         
                         <div style='margin-bottom: 10px;'>
                           <h5 style='color: #2c3e50; margin-bottom: 8px;'><strong>T-test</strong></h5>
                           <p style='margin-bottom: 5px; font-size: 13px;'>
                             <strong>Purpose:</strong> Parametric test comparing means of log-transformed abundance data.
                           </p>
                           <p style='margin-bottom: 0; font-size: 13px; color: #7f8c8d;'>
                             <strong>Best for:</strong> When abundance data is approximately normally distributed after transformation.
                           </p>
                         </div>
                         
                         <hr style='margin: 15px 0; border-color: #bdc3c7;'>
                         <p style='font-style: italic; font-size: 12px; color: #95a5a6; margin-bottom: 0;'>
                           <strong>Note:</strong> All statistical tests include False Discovery Rate (FDR) correction for multiple comparisons.
                         </p>
                       ")
                         )
                       )
                ),
                # Main content area 
                column(width = 9,
                       # Quick stat boxes
                       fluidRow(
                         valueBoxOutput("total_taxa_box", width = 3), 
                         valueBoxOutput("dominant_taxa_box", width = 6),
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
                                                           choices = c("Mix", "Viridis", "Turbo", "Rainbow", "Set3", "Spectral"), 
                                                           selected = "Turbo")
                                             )),
                                      # Plot area 
                                      column(width = 9, 
                                             div(style = "height: 600px; overflow-y: auto;", 
                                                 plotlyOutput("relative_abundance_plot", height = "550px")
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
                           # Tab 2: Heatmap Visualization
                           tabPanel("Heatmap", 
                                    fluidRow(
                                      column(width = 12, 
                                             # Color palette options 
                                             fluidRow(
                                               column(width = 6, 
                                                      offset = 3, 
                                                      selectInput("heatmap_color", "Color Palette:", 
                                                                  choices = c("Viridis", "Magma", "Plasma", "RdBu", "Beyonce", "TayloRSwift"),
                                                                  selected = "RdBu"), 
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
                           )
                         )
                       ), 
                       # Statistical results display box 
                       conditionalPanel(
                         condition = "input.stat_test != 'None' && input.compare_groups == true",
                         box(
                           title = "Statistical Analysis Results", 
                           status = "warning", 
                           solidHeader = TRUE, 
                           width = 12, 
                           collapsible = TRUE, 
                           fluidRow(
                             column(width = 12, 
                                    uiOutput("statistical_results"))
                           )
                         ), 
                         fluidRow(
                           column(width = 12, 
                                  align = "right", 
                                  conditionalPanel(
                                    condition = "output.statistical_results", 
                                    downloadButton("download_statistical_results", 
                                                   "Download Results", 
                                                   class = "btn-warning")
                                  ))
                         )
                       )
                )
              )
      ),
                  
      # TAB 3: Diversity Analysis Tab
      tabItem(tabName = "diversity",
              fluidRow(
                
                #### COLUMN 1: Input Controls
                column(width = 3,
                       
                       # Diversity Settings
                       box(
                         title = "Diversity Settings",
                         status = "primary",
                         solidHeader = TRUE,
                         width = 12,
                         collapsible = TRUE,
                         
                         radioButtons("diversity_type", "Diversity Type:",
                                      choices = c("Alpha Diversity" = "alpha", 
                                                  "Beta Diversity" = "beta"),
                                      selected = "alpha"),
                         checkboxInput("subdivide_samples", "Split samples by condition", value = FALSE),
                         conditionalPanel(
                           condition = "input.subdivide_samples == true",
                           selectInput("condition_column", "Clinical condition to group by:",
                                       choices = c("Ongoing Conditions" = "Ongoing_conditions", 
                                                   "Neurological Disorders" = "Neurological_Disorders", 
                                                   "Allergies" = "Allergies", 
                                                   "Cancer" = "Cancer"),
                                       selected = "Ongoing_conditions")
                           
                         )
                       ),
                       
                       # Visualization Options
                       box(
                         title = "Visualization Options",
                         status = "primary",
                         solidHeader = TRUE,
                         width = 12,
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
                                       choices = c("Bray-Curtis", "Euclidean", "Jaccard",
                                                   "Canberra", "Manhattan", "Kulczynski", "Chord"),
                                       selected = "Bray-Curtis"),
                           selectInput("ordination_method", "Ordination Method:",
                                       choices = c("PCoA", "NMDS", "dbRDA", "t-SNE", "UMAP"),
                                       selected = "PCoA"),
                           checkboxInput("show_group_ellipses", "Show Group Ellipses", value = TRUE)
                         ),
                         checkboxInput("show_diversity_legend", "Show Legend", value = TRUE)
                       ),
                       
                       # Download Button (below the box)
                       fluidRow(
                         column(width = 12, align = "center",
                                
                                conditionalPanel(
                                  condition = "input.diversity_type == 'alpha'",
                                  downloadButton("download_alpha_plot", "Download Alpha Plot", class = "btn-sm btn-info")
                                ),
                                
                                conditionalPanel(
                                  condition = "input.diversity_type == 'beta'",
                                  downloadButton("download_beta_ordination_plot", "Download Beta Plot", class = "btn-sm btn-info"),
                                  downloadButton("download_beta_dendrogram_plot", "Download Dendrogram", class = "btn-sm btn-secondary")
                                )
                         )
                       )
                       
                ),
                
                
                #### COLUMN 2: Output Plots (wider)
                column(width = 6,
                       box(
                         title = "Diversity Analysis",
                         status = "info",
                         solidHeader = TRUE,
                         width = 12,
                         
                         conditionalPanel(
                           condition = "input.diversity_type == 'alpha'",
                           tabsetPanel(
                             tabPanel("Alpha Diversity Plot", plotlyOutput("diversity_main_plot", height = "600px")),
                             tabPanel("Alpha Statistical Summary", verbatimTextOutput("alpha_stats_table"))
                           )
                         ),
                         conditionalPanel(
                           condition = "input.diversity_type == 'beta'",
                           tabsetPanel(
                             tabPanel("Beta Ordination Plot", plotlyOutput("ordination_plot", height = "600px")),
                             tabPanel("Beta Statistical Summary", verbatimTextOutput("beta_stats_table")),
                             tabPanel("Clustering Dendrogram", plotOutput("clustering_plot", height = "600px"))
                           )
                         )
                       )
                ),
                
                #### COLUMN 3: Contextual Explanation (dynamic)
                column(width = 3,
                       conditionalPanel(
                         condition = "input.diversity_type == 'alpha'",
                         box(
                           title = "Explanation of Alpha Metrics",
                           status = "warning",
                           solidHeader = TRUE,
                           width = 12,
                           collapsible = TRUE,
                           HTML("
                       <ul>
                         <li><b>Observed OTUs</b>: The raw count of unique taxa (Operational Taxonomic Units). Reflects richness only, not evenness.</li>
                         <li><b>Shannon Index</b>: Combines richness and evenness. Increases with more species and more uniform distribution.</li>
                         <li><b>Simpson Index</b>: Reflects dominance. Lower values mean a few species dominate; higher values indicate evenness.</li>
                       </ul>
                     ")
                         )
                       ),
                       conditionalPanel(
                         condition = "input.diversity_type == 'beta'",
                         box(
                           title = "Explanation of Beta Metrics",
                           status = "warning",
                           solidHeader = TRUE,
                           width = 12,
                           collapsible = TRUE,
                           HTML("
                       <ul>
                         <li><b>Bray-Curtis</b>: Based on abundance. Sensitive to common species.</li>
                         <li><b>Jaccard</b>: Based on presence/absence. Good for binary comparisons.</li>
                         <li><b>Euclidean</b>: Straight-line distance in multivariate space.</li>
                         <li><b>Manhattan</b>: Absolute differences. Robust to outliers.</li>
                         <li><b>Canberra</b>: Gives weight to rare species. Good for sparse data.</li>
                         <li><b>Kulczynski</b>: Balanced dissimilarity. Accounts for shared vs total taxa.</li>
                         <li><b>Chord</b>: Focuses on proportional differences between samples.</li>
                       </ul>
                       <h4>Ordination Methods</h4>
                       <ul>
                         <li><b>PCoA</b>: Preserves distance dissimilarity. Common and interpretable.</li>
                         <li><b>NMDS</b>: Non-metric, preserves rank of distances. Good for non-linear patterns.</li>
                         <li><b>dbRDA</b>: Allows regression-like constraints in ordination.</li>
                         <li><b>t-SNE</b>: Captures local clusters. Not ideal for global structure.</li>
                         <li><b>UMAP</b>: Preserves both local and global structure. Scalable.</li>
                       </ul>
                     ")
                         )
                       )
                )
              )
      ),
      
      
      
      # TAB 4: Phylogenetic tree
      tabItem(tabName = "phylo",
              fluidRow(
                # Control Panel
                column(width = 3,
                       box(
                         title = "Tree Configuration",
                         status = "primary",
                         solidHeader = TRUE,
                         width = 12,
                         
                         # Button to reload tree if needed
                         actionButton("reload_phylo", "Reload Tree", 
                                      class = "btn-warning btn-sm"),
                         
                         hr(),
                         
                         h4("Visualisation options", style = "font-weight: bold;"),
                         # Tree display options
                         radioButtons("tree_layout", "Tree Layout:",
                                      choices = c("Rectangular" = "rectangular",
                                                  "Circular" = "circular", 
                                                  "Unrooted" = "unrooted"),
                                      selected = "rectangular"),
                         
                         # Tree display options
                         checkboxGroupInput("display_options", "Tree Display Options:",
                                            choices = c("Show Tip Labels" = "show_tip_labels",
                                                        "Show Node Support" = "show_node_labels", 
                                                        "Show Branch Lengths" = "show_branch_length"),
                                            selected = c("show_tip_labels")),
                         
                         # Size controls
                         sliderInput("tree_text_size", "Text Size:",
                                     min = 0.4, max = 5.0, value = 2.8, step = 0.2),
                         
                         hr(),
                         h4("Taxonomic Grouping", style = "font-weight: bold;"),
                         selectInput("tax_level", "Group by Taxonomic Level:",
                                     choices = c("None (show all)" = "none",
                                                 "Kingdom" = "kingdom", 
                                                 "Phylum" = "phylum",
                                                 "Class" = "class",
                                                 "Order" = "order", 
                                                 "Family" = "family",
                                                 "Genus" = "genus"),
                                     selected = "family"),
                         
                         sliderInput("min_group_size", "Minimum group size:",
                                     min = 1, max = 10, value = 2, step = 1),
                         helpText("Groups with fewer members will be shown individually"),
                         
                         checkboxInput("show_group_counts", "Show group counts in labels", value = FALSE),
                         
                         hr(),
                         
                         # Download buttons
                         downloadButton("download_phylo_plot", "Download Tree Plot",
                                        class = "btn-sm btn-success"),
                         br(), br(),
                         downloadButton("download_phylo_data", "Download Tree Data",
                                        class = "btn-sm btn-info")
                       ),
                       
                       # Tree Information Box
                       box(
                         title = "Tree Information",
                         status = "warning",
                         solidHeader = TRUE,
                         width = 12,
                         collapsible = TRUE,
                         
                         htmlOutput("phylo_tree_info")
                       )
                ),
                
                # Main Tree Display
                column(width = 9,
                       box(
                         title = "Phylogenetic Tree Visualization",
                         status = "info",
                         solidHeader = TRUE,
                         width = 12,
                         
                         tabBox(
                           width = 12,
                           height = "980px",
                           
                           # Main tree plot
                           tabPanel("Tree Plot",
                                    div(style = "height: 100%; overflow: auto; text-align: center;",
                                        conditionalPanel(
                                          condition = "output.tree_loaded",
                                          plotOutput("main_phylo_plot", height = "100%", width = "100%")
                                        ),
                                        conditionalPanel(
                                          condition = "!output.tree_loaded",
                                          div(style = "padding: 50px; text-align: center;",
                                              uiOutput("tree_status_message")
                                          )
                                        )
                                    )
                           ),
                           
                           # Tree statistics and data
                           tabPanel("Tree Data",
                                    div(style = "height: 700px; overflow-y: auto;",
                                        conditionalPanel(
                                          condition = "output.tree_loaded",
                                          verbatimTextOutput("tree_summary"),
                                          br(),
                                          h4("Tip Labels:"),
                                          dataTableOutput("tree_tips_table")
                                        ),
                                        conditionalPanel(
                                          condition = "!output.tree_loaded",
                                          div(style = "padding: 50px; text-align: center;",
                                              h4("Tree file not available",
                                                 style = "color: #666;"))
                                        )
                                    )
                           )
                         )
                       )
                )
              )
      ),   
      
      # TAB 5: Lifestyle Tab 
      tabItem(tabName = "multifactor",
              fluidRow(
                box(
                  title = "Parallel Coordinates Plot",
                  status = "primary",
                  solidHeader = TRUE,
                  width = 12,
                  helpText("Visualize the relationship between microbial taxa and selected metadata."),
                  uiOutput("parallel_var_select"),
                  plotlyOutput("parallel_plot", height = "600px"),
                  br(),
                  downloadButton("download_parallel_plot", "Download Plot (PNG)", class = "btn-sm btn-info")
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
    success_message = NULL, 
    available_runs = NULL, 
    folder_scanned = FALSE, 
    automatic_data_loaded = FALSE, 
    manual_data_loaded = FALSE 
  )
  
  ########## AUTOMATIC DATA LOADING FROM HOSPITAL STRUCTURE #######
  # Helper function to extract date from run ID (R01030525 -> 030525)
  extract_date_from_run <- function(run_id) {
    # Extract the date part 
    date_part <- substr(run_id, nchar(run_id) - 5, nchar(run_id))
    # Convert to proper date format
    paste0(substr(date_part, 1, 2), "/", 
           substr(date_part, 3, 4), "/", 
           "20", substr(date_part, 5, 6))
  }
  
  # Helper function to validate directory structure 
  validate_run_directory <- function(run_path) {
    required_dirs <- c("metadata", "qiime_output/relevant_results")
    required_files <- c("metadata/metadata_sample.csv", "metadata/metadata_run.csv")
    
    validation <- list(valid = TRUE, missing = c())
    
    # Check directories 
    for (dir in required_dirs) {
      if (!dir.exists(file.path(run_path, dir))) {
        validation$valid <- FALSE
        validation$missing <- c(validation$missing, paste("Directory:", dir))
      }
    }
    
    # Check files 
    for (file in required_files) {
      if (!file.exists(file.path(run_path, file))) {
        validation$valid <- FALSE
        validation$missing <- c(validation$missing, paste("File:", file))
      }
    }
    return(validation)
  }
  
  # Scan folder for available runs 
  observeEvent(input$scan_folder, {
    data_store$error_message <- NULL
    data_store$folder_scanned <- FALSE 
    
    tryCatch({
      folder_path <- input$data_folder_path
      
      if (is.null(folder_path) || trimws(folder_path) == "" || folder_path == "path/to/hospital/data/structure") {
        data_store$error_message <- "Please enter a valid data folder path"
        showNotification("Please enter a valid data folder path", type = "warning")
        return()
      }
      
      # Clean the path 
      folder_path <- normalizePath(folder_path, mustWork = FALSE)
      
      if (!dir.exists(folder_path)) {
        data_store$error_message <- paste("Directory does not exist:", folder_path)
        showNotification(paste("Directory does not exist:", folder_path), type = "error")
        return()
      }
      
      # Look for runs directory 
      runs_path <- file.path(folder_path, "runs")
      if (!dir.exists(runs_path)) {
        data_store$error_message <- "Runs directory not found. Please check the folder structure."
        showNotification("Runs directory not found in the specified path", type = "error")
        return()
      }
      
      # Find all run directories (pattern: R[01-99][DDMMYY])
      all_dirs <- list.dirs(runs_path, full.names = FALSE, recursive = FALSE)
      run_pattern <- "^R[0-9]+$"
      run_dirs <- all_dirs[grepl(run_pattern, all_dirs)]
      
      if (length(run_dirs) == 0) {
        data_store$error_message <- paste("No valid run directories found in:", runs_path, "\n(expected format: R[01-99][DDMMYY])") 
        showNotification("No valid run directories found", type = "warning")
        return()
      }
      
      # Validate each run directory and create detailed info
      valid_runs <- c()
      run_info <- list()
      
      for (run_id in run_dirs) {
        run_path <- file.path(runs_path, run_id)
        validation <- validate_run_directory(run_path)
        
        if (validation$valid) {
          valid_runs <- c(valid_runs, run_id)
          
          # Try to get additional info from metadata 
          metadata_path <- file.path(run_path, "metadata", "metadata_run.csv")
          run_date <- extract_date_from_run(run_id)
          
          if (file.exists(metadata_path)) {
            tryCatch({
              run_meta <- read.csv(metadata_path, stringsAsFactors = FALSE)
              sample_count <- nrow(run_meta)
              platform <- if("Sequencing_Platform" %in% colnames(run_meta)) {
                paste(unique(run_meta$Sequencing_Platform), collapse = ", ")
              } else "Unknown"
              
              run_info[[run_id]] <- list(
                date = run_date, 
                samples = sample_count,
                platform = platform
              )
            }, error = function(e) {
              run_info[[run_id]] <- list(
                date = run_date, 
                samples = "Unknown", 
                platform = "Unknown"
              )
            })
          } else { 
            run_info[[run_id]] <- list(
              date = run_date, 
              samples = "Unknown", 
              platform = "Unknown"
            )
          }
        } else {
          # Log what's missing for debugging 
          message(paste("Run", run_id, "validation failed. Missing:", paste(validation$missing, collapse = ", ")))
        }
      }
      if (length(valid_runs) == 0) {
        data_store$error_message <- "No valid run directories found with complete data structure"
        showNotification("No complete run directories found. Check that each run has metadata/ and qiime_output/relevant_results/ folders", type = "error")
        return()
      }
      # Create choices with additional information
      run_choices <- valid_runs 
      names(run_choices) <- sapply(valid_runs, function(run) {
        info <- run_info[[run]]
        paste0(run, " | ", info$date, " | ", info$samples, " samples | ", info$platform)
      })
      
      data_store$available_runs <- run_choices
      data_store$folder_scanned <- TRUE
      
      updateSelectizeInput(session, "available_run_ids", 
                           choices = run_choices, 
                           server = TRUE, 
                           options = list(placeholder = "Select a run to load..."))
      
      showNotification(paste("Found", length(valid_runs), "valid run(s)"), type = "message")
    }, error = function(e) {
      data_store$error_message <- paste("Error scanning folder:", e$message)
      showNotification(paste("Scan error:", e$message), type = "error")
      data_store$folder_scanned <- FALSE
    })
  })
  
  # Load automatic data 
  observeEvent(input$load_automatic_data, {
    req(input$available_run_ids)
    
    data_store$error_message <- NULL
    data_store$success_message <- NULL
    
    tryCatch({
      selected_run <- input$available_run_ids 
      folder_path <- input$data_folder_path
      run_path <- file.path(folder_path, "runs", selected_run)
      
      # Load sample metadata 
      sample_meta_path <- file.path(run_path, "metadata", "metadata_sample.csv")
      sample_metadata <- read.csv(sample_meta_path, stringsAsFactors = FALSE)
      
      # Load run metadata 
      run_meta_path <- file.path(run_path, "metadata", "metadata_run.csv")
      run_metadata <- read.csv(run_meta_path, stringsAsFactors = FALSE)
      
      # Try to load taxonomy data from QIIME output 
      taxonomy_paths <- c(
        file.path(run_path, "qiime_output", "relevant_results", "taxonomy_table.tsv"),
        file.path(run_path, "qiime_output", "relevant_results", "feature_table.tsv"),
        file.path(run_path, "qiime_output", "artifacts", "04_taxonomy", "taxonomy.tsv")
      )
      
      taxonomy_data <- NULL
      for(path in taxonomy_paths) {
        if(file.exists(path)) {
          tryCatch({
            taxonomy_data <- read.delim(path, sep = "\t", stringsAsFactors = FALSE)
            break
          }, error = function(e) {
            # Continue to next path
          })
        }
      }
      
      # Try to load control data 
      control_sample_metadata <- NULL
      control_run_metadata <- NULL
      control_taxonomy_data <- NULL
      
      # Look for control data in a separate controls directory
      controls_path <- file.path(run_path, "controls")
      if(dir.exists(controls_path)) {
        tryCatch({
          if(file.exists(file.path(controls_path, "metadata_sample.csv"))) {
            control_sample_metadata <- read.csv(file.path(controls_path, "metadata_sample.csv"),
                                                stringsAsFactors = FALSE)
          }
          if(file.exists(file.path(controls_path, "metadata_run.csv"))) {
            control_run_metadata <- read.csv(file.path(controls_path, "metadata_run.csv"), 
                                             stringsAsFactors = FALSE)
          }
          # Look for control taxonomy data 
          for(path in c("taxonomy_table.tsv", "feature_table.tsv")) {
            control_tax_path <- file.path(controls_path, path)
            if(file.exists(control_tax_path)) {
              control_taxonomy_data <- read.delim(control_tax_path, sep = "\t", stringsAsFactors = FALSE)
              break
            }
          }
        }, error = function(e) {
          # Control data loading failed, continue without controls 
        })
      }
      # Store the loaded data 
      data_store$sample_metadata <- sample_metadata
      data_store$run_metadata <- run_metadata
      data_store$taxonomy_data <- taxonomy_data
      data_store$control_sample_metadata <- control_sample_metadata
      data_store$control_run_metadata <- control_run_metadata
      data_store$control_taxonomy_data <- control_taxonomy_data
      data_store$automatic_data_loaded <- TRUE
      data_store$manual_data_loaded <- FALSE
      
      control_status <- if(!is.null(control_sample_metadata)) {
        paste("Control data loaded:", nrow(control_sample_metadata), "control samples")
      } else {
        "No control data found"
      }
      
      data_store$success_message <- paste0(
        "Automatic data loading successful!\n", 
        "Run: ", selected_run, "\n", 
        "Samples: ", nrow(sample_metadata), "\n",
        "Runs: ", nrow(run_metadata), "\n",
        "Species: ", length(unique(taxonomy_data$Species)), "\n",
        control_status
      )
      showNotification("Automatic data loaded successfully!", type = "message")
    }, error = function(e) {
      data_store$error_message <- paste("Error loading automatic data:", e$message)
      showNotification(paste("Loading error:", e$message), type = "error")
    })
  })
  
  ## Reactive outputs and UI updates (value boxes)
  # Output for folder scan status 
  output$folder_scanned <- reactive({
    data_store$folder_scanned
  })
  outputOptions(output, "folder_scanned", suspendWhenHidden = FALSE)
  
  # Output for automatic data loading status 
  output$automatic_data_loaded <- reactive({
    data_store$automatic_data_loaded 
  })
  outputOptions(output, "automatic_data_loaded", suspendWhenHidden = FALSE)
  
  # Control data status output 
  output$control_data_status <- renderText({
    if(!is.null(data_store$control_sample_metadata)) {
      paste("Control data available:", nrow(data_store$control_sample_metadata), "control samples")
    } else {
      "No control data loaded"
    }
  })
  
  # Automatic data loading status 
  output$automatic_data_loading_status <- renderText({
    if(!is.null(data_store$error_message)) {
      data_store$error_message
    } else if (!is.null(data_store$success_message)) {
      data_store$success_message
    } else {
      "No data loaded yet"
    }
  })
  
  # Value boxes for automatic data 
  output$automatic_total_samples_box <- renderValueBox({
    req(data_store$sample_metadata, data_store$automatic_data_loaded)
    valueBox(
      nrow(data_store$sample_metadata), 
      "Total Samples", 
      icon = icon("users"), 
      color = "blue"
    )
  })
  
  output$automatic_total_species_box <- renderValueBox({
    req(data_store$taxonomy_data, data_store$automatic_data_loaded)
    valueBox(
      length(unique(data_store$taxonomy_data$Species)), 
      "Unique Species", 
      icon = icon("bacteria"), 
      color = "green"
    )
  })
  
  output$automatic_total_runs_box <- renderValueBox({
    req(data_store$run_metadata, data_store$automatic_data_loaded)
    valueBox(
      nrow(data_store$run_metadata), 
      "Total Runs", 
      icon = icon("server"),
      color = "purple"
    )
  })
  
  output$automatic_total_control_samples_box <- renderValueBox({
    req(data_store$automatic_data_loaded)
    control_count <- if(!is.null(data_store$control_sample_metadata)) {
      nrow(data_store$control_sample_metadata)
    } else {
      0
    }
    valueBox(
      control_count, 
      "Control Samples", 
      icon = icon("vial"), 
      color = "yellow"
    )
  })
  
  
  ######### MANUAL DATA UPLOAD FUNCTIONALITY ################
  
  # Generate example data when button is clicked 
  observeEvent(input$generate_example, {
    
    data_store$error_message <- NULL
    
    # Generate sample metadata 
    sample_metadata <- data.frame(
      Sample_ID = paste0("SHM", sprintf("%02d", 1:50), format(Sys.Date(), "%d%m%y")),
      Institution = sample(c("Hospital_del_mar", "Clinical_university", "Medical_center", "Research_institute"), 50, replace = TRUE),
      Department = sample(c("Digestology_unit", "Gastroenterology_department", "Internal_medicine", "Clinical_research"), 50, replace = TRUE),
      Collection_Date = as.character(sample(seq(as.Date('2023-01-01'), as.Date('2025-06-05'), by="day"), 50)),
      `Collection_Storage(Temperature(ºC))` = sample(c("-80", "-20", "4", "Room Temperature"), 50, replace = TRUE),
      Analyst_Processor_Name = sample(c("Sam_1", "Ana_garcia", "John_smith", "Maria_lopez"), 50, replace = TRUE),
      Gender = sample(c("Female", "Male", "Not_specified"), 50, replace = TRUE),
      Age = sample(18:80, 50, replace = TRUE),
      Ongoing_conditions = sample(c("Acid reflux", "Helicobacter pylori–associated gastritis", "Irritable bowel syndrome (IBS)", 
                                    "Diabetes(Type 1)", "Diabetes (Mellitus)", "Metabolic syndromes", 
                                    "Non-alcoholic fatty liver disease (NAFLD)", "Fatty liver disease (FLD)", 
                                    "Inflammatory Bowel Disease (IBD)", "Lupus (SLE)", "Crohn's disease", 
                                    "Rheumatoid arthritis", "Stomach ache"), 50, replace = TRUE),
      Appendix_removed = sample(c("Yes", "No"), 50, replace = TRUE),
      Allergies = sample(c("Peanuts", "Shellfish", "Tree nuts", "Eggs", "Milk", "Unspecified", "Allergy free",
                           "Penicillin", "Amoxicillin", "Asthma", "Seasonal_allergies", "Sulfa_drugs"), 50, replace = TRUE),
      Dietary_Information = sample(c("Omnivore", "Vegetarian", "Vegetarian but eat seafood", "Vegan", 
                                     "Gluten free", "Keto", "Halal", "Kosher", "Paleo"), 50, replace = TRUE),
      Bowel_movement_quality = sample(c("Constipated", "Normal", "Diarrhea"), 50, replace = TRUE),
      Antibiotic_intake = sample(c("Past year", "Year", "Month", "6 months", "Week"), 50, replace = TRUE),
      Medications = sample(c("Antidiabetics", "Probiotics", "Prebiotics", "Laxatives", "Antiemetic", 
                             "PPIs (proton-pump inhibitors)", "Immunosupressors", 
                             "Antidepressors/Antipsicotics/Anxiolytics", "Contraceptives", "Retinoids", 
                             "Antihistamines", "NSAIDs"), 50, replace = TRUE),
      Cancer = sample(c("Yes", "No"), 50, replace = TRUE),
      Body_Mass_Index = sample(c("Underweight<18.5", "Normal Weight 18.5-24.9", "Overweight 25-29.9", "Obese >30"), 50, replace = TRUE),
      Exercise_frequency = sample(c("Never", "Rarely", "1-2 times per week", "3-5 times per week", "Daily"), 50, replace = TRUE),
      Smoking_status = sample(c("Smoker", "Non-smoker"), 50, replace = TRUE),
      Daily_cigarettes = sample(c("1-5", "6-10", "11-15", "16-20", "+20"), 50, replace = TRUE),
      Alcohol_consumption = sample(c("Yes", "No"), 50, replace = TRUE),
      Frequency_of_alcohol_consumption = sample(c("1_per_week", "2_per_week", "3_per_week", "4_per_week", 
                                                  "5_per_week", "6_per_week", "7_per_week", "8_per_week", 
                                                  "9_per_week", "10_per_week"), 50, replace = TRUE),
      Notes_Samples = sample(c("Normal_collection", "Delayed_processing", "Good_quality_sample", 
                               "Standard_procedure", "No_complications"), 50, replace = TRUE),
      stringsAsFactors = FALSE
    )
    
    # Generate run metadata 
    run_metadata <- data.frame(
      RunID = paste0("RHM", sprintf("%02d", 1:30), "_01"),
      Run_accession = paste0("ERR", sample(1000000:9999999, 30)),
      Sequencing_Date = as.character(sample(seq(as.Date('2023-01-15'), as.Date('2025-06-05'), by="day"), 30)),
      Sequencing_Platform = sample(c("MiSeq", "NextSeq", "GrindION", "MinION", "Ion GeneStudio S5"), 30, replace = TRUE),
      Sequencing_Type = sample(c("16S_rRNA", "WGS", "Shotgun metagenomic", "RNA-Seq", "Long-read", "Metatranscriptomics"), 30, replace = TRUE),
      Expected_read_length = sample(c("2x75 bp", "2x100 bp", "2x150 bp", "2x250 bp", "2x300bp"), 30, replace = TRUE),
      Sequencing_depth_target = sample(c("<1 million reads/sample", "1-5 million reads/sample", 
                                         "5-10 million reads/sample", "10-20 million reads/sample", 
                                         "20-50 million reads/sample", "50 million reads/sample"), 30, replace = TRUE),
      Library_preparation_kit = sample(c("Nextera XT", "Nextera DNA Flex", "TruSeq Nano", "NEBNext Ultra II", 
                                         "Swift Biosciences 16S", "QIAseq"), 30, replace = TRUE),
      Technician_name = sample(c("Maria_garcia", "Thomas_johnson", "Sophia_lee", "Michael_brown"), 30, replace = TRUE),
      Notes_Runs = sample(c("Standard_run", "Good_quality_output", "Normal_processing", 
                            "Successful_sequencing", "No_issues_detected"), 30, replace = TRUE),
      Sample_ID = paste0("SHM", sprintf("%02d", sample(1:50, 30, replace = TRUE)), format(Sys.Date(), "%d%m%y")),
      stringsAsFactors = FALSE
    )
    
    # Generate taxonomy data 
    species_names <- c(
      "Bacteroides fragilis", "Escherichia coli", "Lactobacillus acidophilus", "Bifidobacterium longum", 
      "Prevotella copri", "Faecalibacterium prausnitzii", "Akkermansia muciniphila", "Ruminococcus bromii",
      "Clostridium difficile", "Enterococcus faecalis", "Blautia obeum", "Roseburia intestinalis",
      "Streptococcus thermophilus", "Eubacterium rectale", "Methanobrevibacter smithii"
    )
    
    taxonomy_data <- data.frame(
      Sample_ID = rep(sample_metadata$Sample_ID, each = length(species_names)),
      Species = rep(species_names, nrow(sample_metadata)),
      Abundance = round(runif(nrow(sample_metadata) * length(species_names), min = 0, max = 10000))
    )
    
    # Generate CONTROL sample metadata with healthier characteristics 
    control_sample_metadata <- data.frame(
      Sample_ID = paste0("CTRL_SHM", sprintf("%02d", 1:30), format(Sys.Date(), "%d%m%y")),
      Institution = sample(c("Hospital_del_mar", "Clinical_university", "Medical_center"), 30, replace = TRUE),
      Department = sample(c("Digestology_unit", "Gastroenterology_department", "Clinical_research"), 30, replace = TRUE),
      Collection_Date = as.character(sample(seq(as.Date('2023-01-01'), as.Date('2025-06-05'), by="day"), 30)),
      `Collection_Storage(Temperature(ºC))` = sample(c("-80", "-20"), 30, replace = TRUE), # Controls stored at better temps
      Analyst_Processor_Name = sample(c("Sam_1", "Ana_garcia", "John_smith"), 30, replace = TRUE),
      Gender = sample(c("Female", "Male"), 30, replace = TRUE, prob = c(0.5, 0.5)),
      Age = sample(20:65, 30, replace = TRUE), # Healthier age range
      Ongoing_conditions = "Acid reflux", # Controls have minimal conditions
      Appendix_removed = sample(c("Yes", "No"), 30, replace = TRUE, prob = c(0.2, 0.8)),
      Allergies = sample(c("Allergy free", "Seasonal_allergies"), 30, replace = TRUE, prob = c(0.8, 0.2)),
      Dietary_Information = sample(c("Omnivore", "Vegetarian"), 30, replace = TRUE, prob = c(0.7, 0.3)),
      Bowel_movement_quality = sample(c("Normal", "Constipated"), 30, replace = TRUE, prob = c(0.9, 0.1)),
      Antibiotic_intake = sample(c("Past year", "Year"), 30, replace = TRUE, prob = c(0.8, 0.2)),
      Medications = sample(c("Probiotics", "Antihistamines"), 30, replace = TRUE, prob = c(0.3, 0.7)),
      Cancer = "No", # Controls have no cancer
      Body_Mass_Index = sample(c("Normal Weight 18.5-24.9", "Overweight 25-29.9"), 30, replace = TRUE, prob = c(0.8, 0.2)),
      Exercise_frequency = sample(c("3-5 times per week", "Daily"), 30, replace = TRUE, prob = c(0.7, 0.3)),
      Smoking_status = sample(c("Non-smoker"), 30, replace = TRUE), # Controls are non-smokers
      Daily_cigarettes = "1-5", # Placeholder since they don't smoke
      Alcohol_consumption = sample(c("Yes", "No"), 30, replace = TRUE, prob = c(0.3, 0.7)),
      Frequency_of_alcohol_consumption = sample(c("1_per_week", "2_per_week"), 30, replace = TRUE),
      Notes_Samples = sample(c("Control_sample", "Healthy_control", "Standard_control"), 30, replace = TRUE),
      stringsAsFactors = FALSE
    )
    
    # Generate CONTROL run metadata 
    control_run_metadata <- data.frame(
      RunID = paste0("CTRL_RHM", sprintf("%02d", 1:20), "_01"),
      Run_accession = paste0("ERR", sample(1000000:9999999, 20)),
      Sequencing_Date = as.character(sample(seq(as.Date('2023-01-15'), as.Date('2025-06-05'), by="day"), 20)),
      Sequencing_Platform = sample(c("MiSeq", "NextSeq"), 20, replace = TRUE), # Controls use consistent platforms
      Sequencing_Type = sample(c("16S_rRNA", "WGS"), 20, replace = TRUE),
      Expected_read_length = sample(c("2x150 bp", "2x250 bp"), 20, replace = TRUE),
      Sequencing_depth_target = sample(c("10-20 million reads/sample", "20-50 million reads/sample"), 20, replace = TRUE),
      Library_preparation_kit = sample(c("Nextera XT", "TruSeq Nano"), 20, replace = TRUE),
      Technician_name = sample(c("Maria_garcia", "Thomas_johnson"), 20, replace = TRUE),
      Notes_Runs = sample(c("Control_run", "Quality_control", "Standard_control_processing"), 20, replace = TRUE),
      Sample_ID = sample(control_sample_metadata$Sample_ID, 20, replace = TRUE),
      stringsAsFactors = FALSE
    )
    
    # Generate CONTROL taxonomy data with healthier microbiome profile
    control_taxonomy_data <- data.frame(
      Sample_ID = rep(control_sample_metadata$Sample_ID, each = length(species_names)),
      Species = rep(species_names, nrow(control_sample_metadata)),
      Abundance = round(runif(nrow(control_sample_metadata) * length(species_names), min = 0, max = 10000))
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
    data_store$manual_data_loaded <- TRUE
    data_store$automatic_data_loaded <- FALSE
    
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
      data_store$manual_data_loaded <- TRUE
      data_store$automatic_data_loaded <- FALSE 
      
      # Show notification 
      data_store$success_message <- "Data loaded successfully! Ready for analysis."
      showNotification("Data loaded successfully!", type = "message")
    }, error = function(e) {
      # Handle any error that occur during file reading
      data_store$error_message <- paste("Error loading data:", e$message)
      showNotification(paste("Error:", e$message), type = "error")
    })
  })
  
  # Output for manual data loading status
  output$manual_data_loaded <- reactive({
    data_store$manual_data_loaded
  })
  outputOptions(output, "manual_data_loaded", suspendWhenHidden = FALSE)
  
  # Manual data loading status
  output$manual_data_loading_status <- renderText({
    if(!is.null(data_store$error_message)) {
      data_store$error_message
    } else if(!is.null(data_store$success_message)) {
      data_store$success_message
    } else {
      "No data loaded yet"
    }
  })
  
  # Value boxes for manual data 
  output$total_samples_box <- renderValueBox({
    req(data_store$sample_metadata, data_store$manual_data_loaded)
    valueBox(
      nrow(data_store$sample_metadata), 
      "Total Samples", 
      icon = icon("users"), 
      color = "blue"
    )
  })
  
  output$total_species_box <- renderValueBox({
    req(data_store$taxonomy_data, data_store$manual_data_loaded)
    valueBox(
      length(unique(data_store$taxonomy_data$Species)), 
      "Unique Species", 
      icon = icon("bacteria"), 
      color = "green"
    )
  })
  
  output$total_runs_box <- renderValueBox({
    req(data_store$run_metadata, data_store$manual_data_loaded)
    valueBox(
      nrow(data_store$run_metadata), 
      "Total Runs", 
      icon = icon("server"),
      color = "purple"
    )
  })
  
  output$total_control_samples_box <- renderValueBox({
    req(data_store$manual_data_loaded)
    control_count <- if(!is.null(data_store$control_sample_metadata)) {
      nrow(data_store$control_sample_metadata)
    } else {
      0
    }
    valueBox(
      control_count, 
      "Control Samples", 
      icon = icon("vial"), 
      color = "yellow"
    )
  })
  
  # Update sample selection choices when data is loaded
  observe({
    if(data_store$manual_data_loaded && !is.null(data_store$sample_metadata)) {
      sample_choices <- data_store$sample_metadata$Sample_ID
      names(sample_choices) <- paste0(data_store$sample_metadata$Sample_ID, 
                                      " (", data_store$sample_metadata$Collection_Date, ")")
      
      updateSelectizeInput(session, "selected_sample_id", 
                           choices = sample_choices,
                           selected = character(0), 
                           options = list(placeholder = "Choose a sample to analyze..."))
    }
  })
  
  # Clear sample selection button
  observeEvent(input$clear_sample_selection, {
    updateSelectizeInput(session, "selected_sample_id", 
                         selected = character(0))
    showNotification("Sample selection cleared", type = "message")
  })
  
  # Selected sample information display
  output$selected_sample_info <- renderUI({
    # only show info if a sample is actually selected 
    if(is.null(input$selected_sample_id) || input$selected_sample_id == "" ||
       is.null(data_store$sample_metadata)) {
      return(HTML("<p><em>No sample selected</em></p>"))
    }
    
    sample_info <- data_store$sample_metadata[data_store$sample_metadata$Sample_ID == input$selected_sample_id, ]
    
    if(nrow(sample_info) > 0) {
      HTML(paste0(
        "<h4>Sample Information</h4>",
        "<p><strong>Sample ID:</strong> ", sample_info$Sample_ID, "</p>",
        "<p><strong>Collection Date:</strong> ", sample_info$Collection_Date, "</p>",
        "<p><strong>Age:</strong> ", sample_info$Age, "</p>",
        "<p><strong>Gender:</strong> ", sample_info$Gender, "</p>",
        "<p><strong>BMI:</strong> ", sample_info$Body_Mass_Index, "</p>",
        "<p><strong>Conditions:</strong> ", sample_info$Ongoing_conditions, "</p>"
      ))
    }
  })
  
  
  ## DATA PREVIEWS OUTPUTS 
  # Helper function for consistent DataTable styling
  create_enhanced_datatable <- function(data, table_name = "data") {
    # Handle empty or problematic data
    if(is.null(data) || nrow(data) == 0) {
      data <- data.frame(Message = "No data available to display", stringsAsFactors = FALSE)
    }
    
    # Ensure all columns are properly formatted for display
    data <- data.frame(lapply(data, function(x) {
      if(is.numeric(x)) {
        # Format numeric columns for better display
        ifelse(is.na(x), "0", as.character(x))
      } else {
        # Convert factors and other types to character, handle NAs
        ifelse(is.na(x) | x == "", "N/A", as.character(x))
      }
    }), stringsAsFactors = FALSE)
    
    # Create the datatable
    dt <- datatable(
      data, 
      extensions = c('Buttons', 'ColReorder', 'FixedHeader'),
      options = list(
        dom = 'Bfrtip',
        buttons = list(
          list(extend = 'copy', text = 'Copy'),
          list(extend = 'csv', filename = paste0(table_name, '_', Sys.Date())),
          list(extend = 'excel', filename = paste0(table_name, '_', Sys.Date())),
          list(extend = 'pdf', filename = paste0(table_name, '_', Sys.Date())),
          list(extend = 'print', text = 'Print')
        ),
        scrollX = TRUE,
        scrollY = "400px",
        pageLength = 25,
        lengthMenu = list(c(10, 25, 50, 100, -1), c('10', '25', '50', '100', 'All')),
        search = list(regex = TRUE, caseInsensitive = TRUE),
        searchCols = lapply(1:ncol(data), function(x) list(search = "")),
        colReorder = TRUE,
        fixedHeader = TRUE,
        autoWidth = FALSE,
        columnDefs = list(
          list(className = "dt-center", targets = "_all"),
          # Handle wide taxonomy columns better - more restrictive width
          list(width = "150px", targets = which(grepl("taxonomy|Taxon", colnames(data), ignore.case = TRUE)) - 1),
          # Handle numeric abundance columns
          list(className = "dt-right", targets = which(sapply(data, function(x) all(grepl("^[0-9.]+$", x[x != "N/A" & x != "0"])))) - 1),
          # Force word wrapping for long text
          list(className = "dt-body-nowrap", targets = "_all")
        )
      ),
      filter = 'top',  # Add individual column search boxes
      class = 'cell-border stripe hover compact',
      rownames = FALSE
    )
    
    # Apply styling
    dt <- dt %>% formatStyle(columns = colnames(data), fontSize = '12px')
    
    # Apply Type column styling if it exists
    if("Type" %in% colnames(data)) {
      dt <- dt %>% formatStyle('Type', 
                               backgroundColor = styleEqual(c('Sample', 'Control'), 
                                                            c('#e8f4fd', '#fff2cc')),
                               fontWeight = 'bold')
    }
    
    # Apply special styling for taxonomy columns (make them more readable)
    taxonomy_cols <- colnames(data)[grepl("taxonomy|Taxon", colnames(data), ignore.case = TRUE)]
    if(length(taxonomy_cols) > 0) {
      for(col in taxonomy_cols) {
        dt <- dt %>% formatStyle(col, 
                                 fontSize = '10px',
                                 whiteSpace = 'normal',
                                 wordWrap = 'break-word',
                                 maxWidth = '150px',
                                 overflow = 'hidden',
                                 textOverflow = 'ellipsis')
      }
    }
    
    return(dt)
  }
  
  # Sample Metadata Preview with filtering
  output$sample_metadata_preview <- renderDataTable({
    req(data_store$sample_metadata)
    
    data_to_show <- switch(input$sample_metadata_filter,
                           "All" = {
                             combined_data <- data_store$sample_metadata
                             if(!is.null(data_store$control_sample_metadata)) {
                               # Ensure columns match before binding
                               sample_cols <- colnames(data_store$sample_metadata)
                               control_cols <- colnames(data_store$control_sample_metadata)
                               common_cols <- intersect(sample_cols, control_cols)
                               
                               combined_data <- rbind(
                                 cbind(data_store$sample_metadata[, common_cols, drop = FALSE], Type = "Sample"),
                                 cbind(data_store$control_sample_metadata[, common_cols, drop = FALSE], Type = "Control")
                               )
                             } else {
                               combined_data <- cbind(data_store$sample_metadata, Type = "Sample")
                             }
                             combined_data
                           },
                           "Sample" = cbind(data_store$sample_metadata, Type = "Sample"),
                           "Control" = {
                             if(!is.null(data_store$control_sample_metadata)) {
                               cbind(data_store$control_sample_metadata, Type = "Control")
                             } else {
                               data.frame(Message = "No control data available", stringsAsFactors = FALSE)
                             }
                           }
    )
    
    create_enhanced_datatable(data_to_show, "sample_metadata")
  })
  
  # Run Metadata Preview with filtering  
  output$run_metadata_preview <- renderDataTable({
    req(data_store$run_metadata)
    
    data_to_show <- switch(input$run_metadata_filter,
                           "All" = {
                             combined_data <- data_store$run_metadata
                             if(!is.null(data_store$control_run_metadata)) {
                               # Ensure columns match before binding
                               sample_cols <- colnames(data_store$run_metadata)
                               control_cols <- colnames(data_store$control_run_metadata)
                               common_cols <- intersect(sample_cols, control_cols)
                               
                               combined_data <- rbind(
                                 cbind(data_store$run_metadata[, common_cols, drop = FALSE], Type = "Sample"),
                                 cbind(data_store$control_run_metadata[, common_cols, drop = FALSE], Type = "Control")
                               )
                             } else {
                               combined_data <- cbind(data_store$run_metadata, Type = "Sample")
                             }
                             combined_data
                           },
                           "Sample" = cbind(data_store$run_metadata, Type = "Sample"),
                           "Control" = {
                             if(!is.null(data_store$control_run_metadata)) {
                               cbind(data_store$control_run_metadata, Type = "Control")
                             } else {
                               data.frame(Message = "No control data available", stringsAsFactors = FALSE)
                             }
                           }
    )
    
    create_enhanced_datatable(data_to_show, "run_metadata")
  })
  
  # Taxonomy Data Preview with filtering
  output$taxonomy_data_preview <- renderDataTable({
    req(data_store$taxonomy_data)
    
    data_to_show <- switch(input$taxonomy_data_filter,
                           "All" = {
                             combined_data <- data_store$taxonomy_data
                             if(!is.null(data_store$control_taxonomy_data)) {
                               # Ensure columns match before binding
                               sample_cols <- colnames(data_store$taxonomy_data)
                               control_cols <- colnames(data_store$control_taxonomy_data)
                               common_cols <- intersect(sample_cols, control_cols)
                               
                               if(length(common_cols) > 0) {
                                 combined_data <- rbind(
                                   cbind(data_store$taxonomy_data[, common_cols, drop = FALSE], Type = "Sample"),
                                   cbind(data_store$control_taxonomy_data[, common_cols, drop = FALSE], Type = "Control")
                                 )
                               } else {
                                 # If no common columns, show sample data with a warning
                                 combined_data <- cbind(data_store$taxonomy_data, Type = "Sample")
                               }
                             } else {
                               combined_data <- cbind(data_store$taxonomy_data, Type = "Sample")
                             }
                             combined_data[1:min(1000, nrow(combined_data)), ] # Limit for performance
                           },
                           "Sample" = {
                             limited_data <- data_store$taxonomy_data[1:min(1000, nrow(data_store$taxonomy_data)), ]
                             cbind(limited_data, Type = "Sample")
                           },
                           "Control" = {
                             if(!is.null(data_store$control_taxonomy_data)) {
                               limited_data <- data_store$control_taxonomy_data[1:min(1000, nrow(data_store$control_taxonomy_data)), ]
                               cbind(limited_data, Type = "Control")
                             } else {
                               data.frame(Message = "No control data available", stringsAsFactors = FALSE)
                             }
                           }
    )
    
    create_enhanced_datatable(data_to_show, "taxonomy_data")
  })
  
  ##### DATA SUMMARY OUTPUT 
  output$data_summary_dynamic <- renderUI({
    # Check which type of data is loaded
    if(data_store$automatic_data_loaded) {
      create_data_summary("automatic")
    } else if(data_store$manual_data_loaded) {
      create_data_summary("manual")
    } else {
      HTML("<div class='alert alert-warning text-center'>
          <i class='fa fa-exclamation-triangle fa-2x'></i>
          <h4>No Data Loaded</h4>
          <p>Please upload data manually or use automatic data loading to get started.</p>
          </div>")
    }
  })
  
  # Enhanced helper function to create data summary with icons and better styling
  create_data_summary <- function(data_type) {
    req(data_store$sample_metadata, data_store$run_metadata, data_store$taxonomy_data)
    
    # Age distribution with error handling
    age_summary <- tryCatch({
      summary(as.numeric(data_store$sample_metadata$Age))
    }, error = function(e) {
      c(0, 0, 0, 0, 0, 0)
    })
    
    # Gender distribution
    gender_dist <- table(data_store$sample_metadata$Gender)
    
    # Most abundant species (top 10) with error handling - Handle both data formats
    top_species <- tryCatch({
      # Check if this is the long format (generated data) with Abundance and Species columns
      if("Abundance" %in% colnames(data_store$taxonomy_data) && "Species" %in% colnames(data_store$taxonomy_data)) {
        # Long format - aggregate by species
        species_abundance <- aggregate(Abundance ~ Species, data_store$taxonomy_data, sum)
        species_abundance[order(species_abundance$Abundance, decreasing = TRUE)[1:min(10, nrow(species_abundance))], ]
      } 
      # Check if this is OTU table format (wide format with taxonomy column)
      else if("Taxon" %in% colnames(data_store$taxonomy_data) || "taxonomy" %in% colnames(data_store$taxonomy_data)) {
        # Wide format OTU table - columns are samples, rows are OTUs
        taxonomy_col <- if("Taxon" %in% colnames(data_store$taxonomy_data)) "Taxon" else "taxonomy"
        
        # Get sample columns (exclude taxonomy and confidence columns)
        sample_cols <- setdiff(colnames(data_store$taxonomy_data), 
                               c(taxonomy_col, "Confidence", "confidence", "OTU_ID", "otu_id", "#OTU ID"))
        
        if(length(sample_cols) > 0) {
          # Calculate total abundance per OTU across all samples
          abundance_data <- data_store$taxonomy_data[, sample_cols, drop = FALSE]
          # Convert to numeric and handle any non-numeric values
          abundance_data[] <- lapply(abundance_data, function(x) as.numeric(as.character(x)))
          abundance_data[is.na(abundance_data)] <- 0
          
          # Sum across samples for each OTU
          total_abundance <- rowSums(abundance_data, na.rm = TRUE)
          
          # Extract species names from taxonomy string
          taxonomy_strings <- data_store$taxonomy_data[[taxonomy_col]]
          
          # Try to extract species level taxonomy (look for s__ or last level)
          species_names <- sapply(taxonomy_strings, function(tax_str) {
            if(is.na(tax_str) || tax_str == "" || tax_str == "Unassigned") {
              return("Unassigned")
            }
            
            # Split by common delimiters
            tax_levels <- unlist(strsplit(as.character(tax_str), "[;|]"))
            
            # Look for species level (s__)
            species_level <- grep("s__", tax_levels, value = TRUE)
            if(length(species_level) > 0) {
              species_name <- gsub("s__", "", species_level[1])
              species_name <- gsub("_", " ", species_name)
              return(ifelse(species_name == "" || species_name == " ", "Unassigned", species_name))
            }
            
            # If no species level, take the last non-empty level
            tax_levels <- tax_levels[tax_levels != "" & !is.na(tax_levels)]
            if(length(tax_levels) > 0) {
              last_level <- tail(tax_levels, 1)
              # Remove common prefixes
              last_level <- gsub("^[a-z]__", "", last_level)
              last_level <- gsub("_", " ", last_level)
              return(ifelse(last_level == "" || last_level == " ", "Unassigned", last_level))
            }
            
            return("Unassigned")
          })
          
          # Create data frame and aggregate by species (in case multiple OTUs have same species)
          species_abundance <- aggregate(total_abundance, 
                                         by = list(Species = species_names), 
                                         FUN = sum, na.rm = TRUE)
          colnames(species_abundance) <- c("Species", "Abundance")
          
          # Sort and get top 10
          species_abundance <- species_abundance[order(species_abundance$Abundance, decreasing = TRUE), ]
          species_abundance[1:min(10, nrow(species_abundance)), ]
        } else {
          data.frame(Species = "No sample data found", Abundance = 0)
        }
      }
      # Fallback - if neither format is recognized
      else {
        # Try to count occurrences of any species-like column
        possible_species_cols <- c("Species", "species", "Taxon", "taxonomy", "OTU_ID")
        species_col <- intersect(possible_species_cols, colnames(data_store$taxonomy_data))[1]
        
        if(!is.na(species_col)) {
          species_counts <- table(data_store$taxonomy_data[[species_col]])
          top_10 <- head(sort(species_counts, decreasing = TRUE), 10)
          data.frame(Species = names(top_10), Abundance = as.numeric(top_10))
        } else {
          data.frame(Species = "Data format not recognized", Abundance = 0)
        }
      }
    }, error = function(e) {
      data.frame(Species = paste("Error processing data:", e$message), Abundance = 0)
    })
    
    # Sequencing platform distribution
    platform_dist <- if("Sequencing_Platform" %in% colnames(data_store$run_metadata)) {
      table(data_store$run_metadata$Sequencing_Platform)
    } else {
      total_runs <- nrow(data_store$run_metadata)
      c("Unknown" = total_runs)
    }
    
    # date range parsing with multiple format support
    date_range <- tryCatch({
      collection_dates <- data_store$sample_metadata$Collection_Date
      
      # Remove any NA, empty, or null values
      collection_dates <- collection_dates[!is.na(collection_dates) & 
                                             collection_dates != "" & 
                                             !is.null(collection_dates)]
      
      if(length(collection_dates) == 0) {
        return(c("No dates available", "No dates available"))
      }
      
      # Try different date formats commonly used
      date_formats <- c("%Y-%m-%d", "%m/%d/%Y", "%d/%m/%Y", "%Y/%m/%d", 
                        "%m-%d-%Y", "%d-%m-%Y", "%B %d, %Y", "%d %B %Y",
                        "%Y%m%d", "%m/%d/%y", "%d/%m/%y")
      
      parsed_dates <- NULL
      
      # Try each format until one works
      for(fmt in date_formats) {
        parsed_dates <- tryCatch({
          as.Date(collection_dates, format = fmt)
        }, error = function(e) NULL)
        
        # If we got valid dates (not all NA), use this format
        if(!is.null(parsed_dates) && sum(!is.na(parsed_dates)) > 0) {
          break
        }
      }
      
      # If no format worked, try automatic parsing
      if(is.null(parsed_dates) || sum(!is.na(parsed_dates)) == 0) {
        # Try lubridate-style parsing if available, otherwise use as.Date without format
        parsed_dates <- tryCatch({
          # First try base R automatic parsing
          as.Date(collection_dates)
        }, error = function(e) {
          # If that fails, try converting to character first
          tryCatch({
            as.Date(as.character(collection_dates))
          }, error = function(e2) {
            # Last resort: try parsing as numeric (Excel dates)
            tryCatch({
              as.Date(as.numeric(collection_dates), origin = "1899-12-30")
            }, error = function(e3) {
              rep(NA, length(collection_dates))
            })
          })
        })
      }
      
      # Filter out NA dates and get range
      valid_dates <- parsed_dates[!is.na(parsed_dates)]
      
      if(length(valid_dates) == 0) {
        return(c("Invalid date format", "Invalid date format"))
      }
      
      date_range <- range(valid_dates, na.rm = TRUE)
      format(date_range, "%Y-%m-%d")
      
    }, error = function(e) {
      c("Error parsing dates", "Error parsing dates")
    })
    
    # Calculate total samples and runs
    total_samples <- nrow(data_store$sample_metadata)
    total_runs <- nrow(data_store$run_metadata)
    
    # Calculate collection days more safely
    collection_days <- tryCatch({
      if(all(date_range != "Unknown") && 
         all(date_range != "No dates available") && 
         all(date_range != "Invalid date format") &&
         all(date_range != "Error parsing dates")) {
        start_date <- as.Date(date_range[1])
        end_date <- as.Date(date_range[2])
        as.numeric(end_date - start_date + 1)
      } else {
        "N/A"
      }
    }, error = function(e) {
      "N/A"
    })
    
    # Create enhanced HTML summary with icons and cards
    HTML(paste0(
      # Details Section
      "<div class='row'>",
      
      # Left column - Demographics & Collection Info
      "<div class='col-md-6'>",
      "<div class='box box-info'>",
      "<div class='box-header with-border'>",
      "<h3 class='box-title'><i class='fa fa-chart-pie'></i> Demographics & Collection</h3>",
      "</div>",
      "<div class='box-body'>",
      
      # Collection Period
      "<div class='row mb-3'>",
      "<div class='col-xs-12'>",
      "<strong><i class='fa fa-calendar'></i> Collection Period:</strong><br>",
      "<span class='text-muted'>", date_range[1], " to ", date_range[2], "</span>",
      "</div></div>",
      
      # Age Statistics
      "<div class='row mb-3'>",
      "<div class='col-xs-6'>",
      "<strong><i class='fa fa-birthday-cake'></i> Age Range:</strong><br>",
      "<span class='text-muted'>", round(age_summary[1], 1), " - ", round(age_summary[6], 1), " years</span>",
      "</div>",
      "<div class='col-xs-6'>",
      "<strong><i class='fa fa-calculator'></i> Mean Age:</strong><br>",
      "<span class='text-muted'>", round(age_summary[4], 1), " years</span>",
      "</div></div>",
      
      # Gender Distribution
      "<div class='row mb-3'>",
      "<div class='col-xs-12'>",
      "<strong><i class='fa fa-venus-mars'></i> Gender Distribution:</strong><br>",
      paste0("<span class='label label-", 
             c("primary", "success", "warning")[1:length(gender_dist)], 
             " mr-2'>", names(gender_dist), ": ", gender_dist, 
             " (", round(gender_dist/sum(gender_dist)*100, 1), "%)</span>", 
             collapse = " "),
      "</div></div>",
      
      # Sequencing Info
      "<div class='row mb-3'>",
      "<div class='col-xs-12'>",
      "<strong><i class='fa fa-server'></i> Sequencing Platforms:</strong><br>",
      paste0("<span class='label label-info mr-2'>", 
             names(platform_dist), ": ", platform_dist, " runs</span>", 
             collapse = " "),
      "</div></div>",
      
      "</div></div></div>",
      
      # Right column - Top Species
      "<div class='col-md-6'>",
      "<div class='box box-success'>",
      "<div class='box-header with-border'>",
      "<h3 class='box-title'><i class='fa fa-list-ol'></i> Top Abundant Species</h3>",
      "</div>",
      "<div class='box-body'>",
      "<div class='table-responsive'>",
      "<table class='table table-striped table-condensed'>",
      "<thead><tr><th>#</th><th>Species</th><th>Total Abundance</th></tr></thead>",
      "<tbody>",
      
      # Top species list
      paste0(sapply(1:nrow(top_species), function(i) {
        paste0("<tr>",
               "<td><span class='badge bg-blue'>", i, "</span></td>",
               "<td><i class='fa fa-bacteria text-green'></i> ", 
               if(nchar(as.character(top_species$Species[i])) > 35) {
                 paste0(substr(as.character(top_species$Species[i]), 1, 35), "...")
               } else {
                 as.character(top_species$Species[i])
               }, "</td>",
               "<td><strong>", format(top_species$Abundance[i], big.mark = ","), "</strong></td>",
               "</tr>")
      }), collapse = ""),
      
      "</tbody></table>",
      "</div></div></div></div>",
      
      "</div>",
      
      # Additional Info Row
      "<div class='row mt-3'>",
      "<div class='col-md-12'>",
      "<div class='alert alert-info'>",
      "<i class='fa fa-info-circle'></i> ",
      "<strong>Data Summary:</strong> This overview shows key statistics from your microbiome data. ",
      "The taxonomy data has been processed to extract species-level information from ",
      ifelse("Abundance" %in% colnames(data_store$taxonomy_data), 
             "abundance tables", "OTU/ASV tables with taxonomic assignments"), ". ",
      "Use the tabs above to explore detailed metadata and taxonomy information.",
      "</div></div></div>"
    ))
  }
  
  
  ## ADDITIONAL HELPER FUNCTIONS 
  # Function to validate uploaded files
  validate_uploaded_file <- function(file_path, expected_columns = NULL) {
    tryCatch({
      if(grepl("\\.csv$", file_path)) {
        data <- read.csv(file_path, stringsAsFactors = FALSE)
      } else {
        data <- read.delim(file_path, sep = "\t", stringsAsFactors = FALSE)
      }
      
      if(!is.null(expected_columns)) {
        missing_cols <- setdiff(expected_columns, colnames(data))
        if(length(missing_cols) > 0) {
          return(list(valid = FALSE, message = paste("Missing columns:", paste(missing_cols, collapse = ", "))))
        }
      }
      
      return(list(valid = TRUE, data = data))
    }, error = function(e) {
      return(list(valid = FALSE, message = paste("Error reading file:", e$message)))
    })
  }
  
  # Function to check data consistency
  check_data_consistency <- function() {
    req(data_store$sample_metadata, data_store$run_metadata, data_store$taxonomy_data)
    
    issues <- c()
    
    # Check if Sample_IDs in taxonomy data match those in sample metadata
    sample_ids_meta <- data_store$sample_metadata$Sample_ID
    sample_ids_tax <- unique(data_store$taxonomy_data$Sample_ID)
    
    missing_in_taxonomy <- setdiff(sample_ids_meta, sample_ids_tax)
    missing_in_metadata <- setdiff(sample_ids_tax, sample_ids_meta)
    
    if(length(missing_in_taxonomy) > 0) {
      issues <- c(issues, paste("Samples in metadata but not in taxonomy:", paste(missing_in_taxonomy, collapse = ", ")))
    }
    
    if(length(missing_in_metadata) > 0) {
      issues <- c(issues, paste("Samples in taxonomy but not in metadata:", paste(missing_in_metadata, collapse = ", ")))
    }
    
    # Check for duplicate Sample_IDs
    if(any(duplicated(data_store$sample_metadata$Sample_ID))) {
      issues <- c(issues, "Duplicate Sample_IDs found in sample metadata")
    }
    
    if(any(duplicated(data_store$run_metadata$Run_ID))) {
      issues <- c(issues, "Duplicate Run_IDs found in run metadata")
    }
    
    return(issues)
  }
  
  # Reactive expression for data validation
  data_validation <- reactive({
    if(data_store$manual_data_loaded || data_store$automatic_data_loaded) {
      check_data_consistency()
    } else {
      NULL
    }
  })
  
  # Show data validation warnings if any
  observe({
    validation_issues <- data_validation()
    if(!is.null(validation_issues) && length(validation_issues) > 0) {
      showNotification(
        paste("Data validation warnings:", paste(validation_issues, collapse = "; ")),
        type = "warning",
        duration = 10
      )
    }
  })
  
  
  
  #### TAXONOMY TAB ######
  
  # Define colors 
  taxonomic_colors <- list(
    "Viridis" = scale_fill_viridis_d(),
    "Mix" = scale_fill_discrete(), 
    "Turbo" = scale_fill_viridis_d(option = "turbo"), 
    "Rainbow" = scale_fill_manual(values = rainbow(25)), 
    "Set3" = scale_fill_brewer(palette = "Set3"), 
    "Spectral" = scale_fill_brewer(palette = "Spectral")
    
    
  )
  
  # Function to extract taxonomic levels from species names
  extract_taxonomy <- reactive({
    req(data_store$taxonomy_data, data_store$control_taxonomy_data)
    
    # Combine patient and control taxonomy data 
    all_taxonomy <- rbind(data_store$taxonomy_data, data_store$control_taxonomy_data)
    
    # Check if we have real OTU data with taxonomic columns or simplified data
    if ("Taxon" %in% colnames(all_taxonomy) && any(!is.na(all_taxonomy$Taxon))) {
      # REAL OTU DATA APPROACH - Parse actual taxonomic strings
      
      # Get unique taxa strings
      taxa_list <- unique(all_taxonomy$Taxon[!is.na(all_taxonomy$Taxon) & all_taxonomy$Taxon != ""])
      
      # Create a dataframe to store taxonomic hierarchy 
      taxonomy_df <- data.frame(Taxon = taxa_list, stringsAsFactors = FALSE)
      
      # Initialize taxonomic columns
      taxonomy_df$Species <- NA
      taxonomy_df$Genus <- NA
      taxonomy_df$Family <- NA
      taxonomy_df$Order <- NA
      taxonomy_df$Class <- NA
      taxonomy_df$Phylum <- NA
      taxonomy_df$Kingdom <- NA
      
      # Parse taxonomic strings (format: k__Kingdom;p__Phylum;c__Class;o__Order;f__Family;g__Genus;s__Species)
      for (i in 1:nrow(taxonomy_df)) {
        taxon_string <- taxonomy_df$Taxon[i]
        
        if (!is.na(taxon_string) && taxon_string != "") {
          # Split by semicolon
          tax_parts <- strsplit(taxon_string, ";")[[1]]
          
          for (part in tax_parts) {
            if (grepl("^k__", part)) {
              taxonomy_df$Kingdom[i] <- gsub("^k__", "", part)
            } else if (grepl("^p__", part)) {
              taxonomy_df$Phylum[i] <- gsub("^p__", "", part)
            } else if (grepl("^c__", part)) {
              taxonomy_df$Class[i] <- gsub("^c__", "", part)
            } else if (grepl("^o__", part)) {
              taxonomy_df$Order[i] <- gsub("^o__", "", part)
            } else if (grepl("^f__", part)) {
              taxonomy_df$Family[i] <- gsub("^f__", "", part)
            } else if (grepl("^g__", part)) {
              taxonomy_df$Genus[i] <- gsub("^g__", "", part)
            } else if (grepl("^s__", part)) {
              taxonomy_df$Species[i] <- gsub("^s__", "", part)
            }
          }
        }
      }
      
      # Clean up empty or undefined taxonomic assignments
      taxonomy_df[taxonomy_df == "" | is.na(taxonomy_df)] <- "Unassigned"
      
      # For species level, create a more readable name if available
      taxonomy_df$Species <- ifelse(
        taxonomy_df$Species == "Unassigned" & taxonomy_df$Genus != "Unassigned",
        paste(taxonomy_df$Genus, "sp."),
        taxonomy_df$Species
      )
      
      return(taxonomy_df)
      
    } else {
      # SIMPLIFIED APPROACH for generated data (existing code)
      
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
    }
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
    
    # Check if we have real OTU data or simplified data
    if ("Taxon" %in% colnames(filtered_data)) {
      # For real OTU data, merge by Taxon
      filtered_data <- merge(filtered_data, taxonomy_hierarchy, by = "Taxon", all.x = TRUE)
    } else {
      # For simplified data, merge by Species
      filtered_data <- merge(filtered_data, taxonomy_hierarchy, by = "Species", all.x = TRUE)
    }
    
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
    
    # Determine the abundance column name
    abundance_col <- if ("Abundance" %in% colnames(tax_data)) "Abundance" else "Count"
    
    # Group by Sample_ID to calculate total abundance per sample
    sample_totals <- aggregate(tax_data[[abundance_col]], 
                               by = list(Sample_ID = tax_data$Sample_ID), 
                               FUN = sum)
    colnames(sample_totals)[2] <- paste0(abundance_col, "_total")
    
    # Merge with original data to calculate relative abundance
    tax_data_with_totals <- merge(tax_data, sample_totals, by = "Sample_ID")
    
    # Calculate relative abundance
    tax_data_with_totals$RelativeAbundance <- (tax_data_with_totals[[abundance_col]] / 
                                                 tax_data_with_totals[[paste0(abundance_col, "_total")]]) * 100
    
    # Update the Abundance column to use the correct column name for downstream processing
    if (abundance_col != "Abundance") {
      tax_data_with_totals$Abundance <- tax_data_with_totals[[abundance_col]]
    }
    
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
  
  # Categorize categorical variables -> Age 
  categorize_metadata <- reactive({
    req(merged_metadata())
    
    metadata <- merged_metadata()
    
    if ("Age" %in% colnames(metadata) && is.numeric(metadata$Age)) {
      # Create age categories 
      metadata$Age_Category <- cut(metadata$Age, 
                                   breaks = c(-Inf, 30, 45, 60, 75, Inf), 
                                   labels = c("≤30", "31-45", "46-60", "61-75", ">75"), 
                                   right = TRUE)
      # Convert to character to avoid factor issues 
      metadata$Age_Category <- as.character(metadata$Age_Category)
    }
    return(metadata)
  })
  
  # Add a group column to taxonomy data for comparison
  taxonomy_with_groups <- reactive({
    req(taxonomy_relative_abundance(), categorize_metadata())
    
    tax_data <- taxonomy_relative_abundance()
    metadata <- categorize_metadata()
    
    # Merge taxonomy data with metadata to get group information
    result <- merge(tax_data, metadata[, c("Sample_ID", "Group")], by = "Sample_ID", all.x = TRUE)
    
    # If metadata_group is selected, add that grouping as well
    if (input$metadata_group != "None") {
      # Check if we should use the categorical version of Age
      if (input$metadata_group == "Age" && "Age_Category" %in% colnames(metadata)) {
        # Use Age_Category instead of Age
        metadata_subset <- metadata[, c("Sample_ID", "Age_Category")]
        colnames(metadata_subset)[2] <- "Age"  # Rename for consistency
        result <- merge(result, metadata_subset, by = "Sample_ID", all.x = TRUE)
      } else {
        # Use the original metadata column
        metadata_subset <- metadata[, c("Sample_ID", input$metadata_group)]
        result <- merge(result, metadata_subset, by = "Sample_ID", all.x = TRUE)
      }
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
  
  # Generate the interactive relative abundance plot
  output$relative_abundance_plot <- renderPlotly({
    req(taxonomy_with_groups())
    
    # Get data with groups
    plot_data <- taxonomy_with_groups()
    
    # Check if we have a metadata group selected
    if (input$metadata_group != "None") {
      metadata_col <- input$metadata_group
      
      # Create a combined grouping variable for X-axis positioning
      # This will group samples by metadata category and patient/control group if available
      if (input$compare_groups && "Group" %in% colnames(plot_data)) {
        # Create combined grouping: Group_MetadataValue (e.g., "Control_Female", "Patient_Male")
        plot_data$XGroup <- paste(plot_data$Group, plot_data[[metadata_col]], sep = "_")
        plot_data$XGroup <- factor(plot_data$XGroup, 
                                   levels = unique(plot_data$XGroup[order(plot_data$Group, plot_data[[metadata_col]])]))
      } else {
        # Use just the metadata group
        plot_data$XGroup <- plot_data[[metadata_col]]
        plot_data$XGroup <- factor(plot_data$XGroup, levels = unique(sort(plot_data$XGroup)))
      }
      
      # Order samples within each group
      plot_data <- plot_data[order(plot_data$XGroup, plot_data$Sample_ID), ]
      plot_data$Sample_ID <- factor(plot_data$Sample_ID, levels = unique(plot_data$Sample_ID))
      
      # Create the plot with samples on X-axis but grouped by metadata
      p <- ggplot(plot_data, aes(x = Sample_ID, y = RelativeAbundance, fill = TaxonomicGroup,
                                 text = paste("Sample:", Sample_ID,
                                              "<br>Taxonomic Group:", TaxonomicGroup,
                                              "<br>Relative Abundance:", round(RelativeAbundance, 2), "%",
                                              "<br>Metadata Group:", .data[[metadata_col]],
                                              if(input$compare_groups && "Group" %in% colnames(plot_data)) 
                                                paste("<br>Group:", Group) else ""))) +
        geom_bar(stat = "identity") +
        theme_minimal() +
        theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 6),
              strip.text = element_text(size = 10, face = "bold"),
              legend.position = "bottom",
              legend.text = element_text(size = 8),
              panel.spacing = unit(0.5, "lines"),
              panel.border = element_rect(color = "black", fill = NA, size = 0.5)) +
        labs(x = "Sample ID", y = "Relative Abundance (%)", 
             title = paste("Relative Abundance of", input$tax_level, "by", 
                           ifelse(metadata_col == "Age", "Age Category", metadata_col)),
             fill = input$tax_level) +
        get_color_scheme()
      
      # Add faceting to group samples by metadata categories
      if (input$compare_groups && "Group" %in% colnames(plot_data)) {
        # Create nested faceting: Group ~ MetadataCategory
        p <- p + facet_grid(reformulate(metadata_col, "Group"), scales = "free_x", space = "free_x") +
          labs(title = paste("Relative Abundance of", input$tax_level, "by Group and", 
                             ifelse(metadata_col == "Age", "Age Category", metadata_col)))
      } else {
        # Simple faceting by metadata category
        p <- p + facet_wrap(reformulate(metadata_col), scales = "free_x") +
          labs(title = paste("Relative Abundance of", input$tax_level, "by", 
                             ifelse(metadata_col == "Age", "Age Category", metadata_col)))
      }
      
    } else if (input$compare_groups && "Group" %in% colnames(plot_data)) {
      # Only group comparison without metadata grouping
      plot_data <- plot_data[order(plot_data$Group, plot_data$Sample_ID), ]
      plot_data$Sample_ID <- factor(plot_data$Sample_ID, levels = unique(plot_data$Sample_ID))
      
      p <- ggplot(plot_data, aes(x = Sample_ID, y = RelativeAbundance, fill = TaxonomicGroup,
                                 text = paste("Sample:", Sample_ID,
                                              "<br>Taxonomic Group:", TaxonomicGroup,
                                              "<br>Relative Abundance:", round(RelativeAbundance, 2), "%",
                                              "<br>Group:", Group))) +
        geom_bar(stat = "identity") +
        theme_minimal() +
        theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 8),
              legend.position = "bottom",
              legend.text = element_text(size = 8),
              panel.border = element_rect(color = "black", fill = NA, size = 0.5)) +
        labs(x = "Sample ID", y = "Relative Abundance (%)", 
             title = paste("Relative Abundance of", input$tax_level, "by Group"),
             fill = input$tax_level) +
        get_color_scheme() +
        facet_wrap(~ Group, scales = "free_x")
      
    } else {
      # Default case - no grouping
      plot_data$Sample_ID <- factor(plot_data$Sample_ID, levels = unique(sort(plot_data$Sample_ID)))
      
      p <- ggplot(plot_data, aes(x = Sample_ID, y = RelativeAbundance, fill = TaxonomicGroup,
                                 text = paste("Sample:", Sample_ID,
                                              "<br>Taxonomic Group:", TaxonomicGroup,
                                              "<br>Relative Abundance:", round(RelativeAbundance, 2), "%"))) +
        geom_bar(stat = "identity") +
        theme_minimal() +
        theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 8),
              legend.position = "bottom",
              legend.text = element_text(size = 8)) +
        labs(x = "Sample ID", y = "Relative Abundance (%)", 
             title = paste("Relative Abundance of", input$tax_level),
             fill = input$tax_level) +
        get_color_scheme()
    }
    
    # Convert to plotly with custom tooltip and configuration
    plotly_plot <- ggplotly(p, tooltip = "text") %>%
      layout(
        showlegend = TRUE,
        legend = list(
          orientation = "h",
          x = 0,
          y = -0.6,
          font = list(size = 10)
        ),
        margin = list(b = 150, l = 50, r = 50, t = 80),
        hovermode = "closest", 
        xaxis = list(
          title = list( standoff = 30 )
        )
      ) %>%
      config(
        displayModeBar = TRUE,
        modeBarButtonsToRemove = c("select2d", "lasso2d", "autoScale2d", "hoverClosestCartesian", "hoverCompareCartesian"),
        displaylogo = FALSE,
        toImageButtonOptions = list(
          format = "png",
          filename = paste("taxonomy_", input$tax_level, "_plot_", format(Sys.time(), "%Y%m%d_%H%M%S"), sep = ""),
          height = 600,
          width = 1000,
          scale = 2
        )
      )
    
    return(plotly_plot)
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
      
      # Apply the same reordering logic as in the display plot
      if (input$metadata_group != "None") {
        metadata_col <- input$metadata_group
        
        if (input$compare_groups && "Group" %in% colnames(plot_data)) {
          plot_data$XGroup <- paste(plot_data$Group, plot_data[[metadata_col]], sep = "_")
          plot_data$XGroup <- factor(plot_data$XGroup, 
                                     levels = unique(plot_data$XGroup[order(plot_data$Group, plot_data[[metadata_col]])]))
        } else {
          plot_data$XGroup <- plot_data[[metadata_col]]
          plot_data$XGroup <- factor(plot_data$XGroup, levels = unique(sort(plot_data$XGroup)))
        }
        
        plot_data <- plot_data[order(plot_data$XGroup, plot_data$Sample_ID), ]
        plot_data$Sample_ID <- factor(plot_data$Sample_ID, levels = unique(plot_data$Sample_ID))
        
        # Set up the plot with the same parameters as the displayed plot
        p <- ggplot(plot_data, aes(x = Sample_ID, y = RelativeAbundance, fill = TaxonomicGroup)) +
          geom_bar(stat = "identity") +
          theme_minimal() +
          theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 8),
                legend.position = "bottom",
                panel.spacing = unit(0.5, "lines"),
                panel.border = element_rect(color = "black", fill = NA, size = 0.5)) +
          labs(x = "Sample ID", y = "Relative Abundance (%)", 
               title = paste("Relative Abundance of", input$tax_level, "by", 
                             ifelse(metadata_col == "Age", "Age Category", metadata_col)),
               fill = input$tax_level) +
          get_color_scheme()
        
        # Add faceting
        if (input$compare_groups && "Group" %in% colnames(plot_data)) {
          p <- p + facet_grid(reformulate(metadata_col, "Group"), scales = "free_x", space = "free_x")
        } else {
          p <- p + facet_wrap(reformulate(metadata_col), scales = "free_x")
        }
        
      } else if (input$compare_groups && "Group" %in% colnames(plot_data)) {
        plot_data <- plot_data[order(plot_data$Group, plot_data$Sample_ID), ]
        plot_data$Sample_ID <- factor(plot_data$Sample_ID, levels = unique(plot_data$Sample_ID))
        
        p <- ggplot(plot_data, aes(x = Sample_ID, y = RelativeAbundance, fill = TaxonomicGroup)) +
          geom_bar(stat = "identity") +
          theme_minimal() +
          theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 8),
                legend.position = "bottom",
                panel.border = element_rect(color = "black", fill = NA, size = 0.5)) +
          labs(x = "Sample ID", y = "Relative Abundance (%)", 
               title = paste("Relative Abundance of", input$tax_level, "by Group"),
               fill = input$tax_level) +
          get_color_scheme() +
          facet_wrap(~ Group, scales = "free_x")
      } else {
        plot_data$Sample_ID <- factor(plot_data$Sample_ID, levels = unique(sort(plot_data$Sample_ID)))
        
        p <- ggplot(plot_data, aes(x = Sample_ID, y = RelativeAbundance, fill = TaxonomicGroup)) +
          geom_bar(stat = "identity") +
          theme_minimal() +
          theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 8),
                legend.position = "bottom") +
          labs(x = "Sample ID", y = "Relative Abundance (%)", 
               title = paste("Relative Abundance of", input$tax_level),
               fill = input$tax_level) +
          get_color_scheme()
      }
      
      # Save the plot to the file
      ggsave(file, plot = p, width = 14, height = 8, dpi = 300)
    }
  )
  ## STATISTICAL TESTS 
  # Prepare data for statistical testing 
  prepare_statistical_data <- reactive({
    req(taxonomy_with_groups())
    
    # Only proceed if compare_groups is TRUE and we have both groups
    if (!input$compare_groups || input$stat_test == "None") {
      return(NULL)
    }
    
    tax_data <- taxonomy_with_groups()
    
    # Check if we have both Patient and Control groups
    if (!"Group" %in% colnames(tax_data) || 
        length(unique(tax_data$Group)) < 2) {
      return(NULL)
    }
    
    # Create abundance matrix (samples as rows, taxa as columns)
    abundance_matrix <- tax_data %>%
      select(Sample_ID, TaxonomicGroup, Abundance) %>%
      tidyr::pivot_wider(names_from = TaxonomicGroup, 
                         values_from = Abundance, 
                         values_fill = 0) %>%
      column_to_rownames("Sample_ID")
    
    # Get sample metadata
    sample_metadata <- tax_data %>%
      select(Sample_ID, Group) %>%
      distinct() %>%
      column_to_rownames("Sample_ID")
    
    # Ensure the row order matches between abundance matrix and metadata
    sample_metadata <- sample_metadata[rownames(abundance_matrix), , drop = FALSE]
    
    return(list(
      abundance_matrix = abundance_matrix,
      metadata = sample_metadata,
      original_data = tax_data
    ))
  })
  
  # Perform statistical tests
  perform_statistical_test <- reactive({
    req(input$stat_test != "None", prepare_statistical_data())
    
    stat_data <- prepare_statistical_data()
    
    if (is.null(stat_data)) {
      return(NULL)
    }
    
    abundance_matrix <- stat_data$abundance_matrix
    metadata <- stat_data$metadata
    
    tryCatch({
      switch(input$stat_test,
             "Wilcoxon" = {
               # Wilcoxon rank-sum test for each taxonomic group
               results_list <- list()
               
               for (taxon in colnames(abundance_matrix)) {
                 patient_data <- abundance_matrix[metadata$Group == "Patient", taxon]
                 control_data <- abundance_matrix[metadata$Group == "Control", taxon]
                 
                 # Only test if both groups have non-zero values
                 if (sum(patient_data > 0) >= 3 && sum(control_data > 0) >= 3) {
                   test_result <- wilcox.test(patient_data, control_data)
                   
                   results_list[[taxon]] <- data.frame(
                     Taxon = taxon,
                     Test = "Wilcoxon",
                     P_value = test_result$p.value,
                     W_statistic = test_result$statistic,
                     Mean_Patients = mean(patient_data),
                     Mean_Controls = mean(control_data),
                     stringsAsFactors = FALSE
                   )
                 }
               }
               
               if (length(results_list) > 0) {
                 results_df <- do.call(rbind, results_list)
                 results_df$P_adjusted <- p.adjust(results_df$P_value, method = "BH")
                 results_df$Significant <- results_df$P_adjusted < 0.05
                 
                 return(list(
                   test_type = "Wilcoxon Rank-Sum Test",
                   results = results_df,
                   summary = paste("Tested", nrow(results_df), "taxa.",
                                   "Significant differences found in", 
                                   sum(results_df$Significant), "taxa (FDR < 0.05).")
                 ))
               }
             },
             
             "t-test" = {
               # T-test for each taxonomic group (log-transformed data)
               results_list <- list()
               
               # Add small pseudocount to avoid log(0)
               log_abundance <- log(abundance_matrix + 1)
               
               for (taxon in colnames(log_abundance)) {
                 patient_data <- log_abundance[metadata$Group == "Patient", taxon]
                 control_data <- log_abundance[metadata$Group == "Control", taxon]
                 
                 # Only test if both groups have sufficient samples
                 if (length(patient_data) >= 3 && length(control_data) >= 3) {
                   test_result <- t.test(patient_data, control_data)
                   
                   results_list[[taxon]] <- data.frame(
                     Taxon = taxon,
                     Test = "T-test",
                     P_value = test_result$p.value,
                     T_statistic = test_result$statistic,
                     Mean_log_Patients = mean(patient_data),
                     Mean_log_Controls = mean(control_data),
                     stringsAsFactors = FALSE
                   )
                 }
               }
               
               if (length(results_list) > 0) {
                 results_df <- do.call(rbind, results_list)
                 results_df$P_adjusted <- p.adjust(results_df$P_value, method = "BH")
                 results_df$Significant <- results_df$P_adjusted < 0.05
                 
                 return(list(
                   test_type = "T-test (log-transformed)",
                   results = results_df,
                   summary = paste("Tested", nrow(results_df), "taxa.",
                                   "Significant differences found in", 
                                   sum(results_df$Significant), "taxa (FDR < 0.05).")
                 ))
               }
             },
             
             "PERMANOVA" = {
               # PERMANOVA test using vegan package
               if (!requireNamespace("vegan", quietly = TRUE)) {
                 return(list(
                   test_type = "PERMANOVA",
                   error = "vegan package not installed. Please install with: install.packages('vegan')"
                 ))
               }
               
               # Calculate distance matrix (Bray-Curtis dissimilarity)
               dist_matrix <- vegan::vegdist(abundance_matrix, method = "bray")
               
               # Perform PERMANOVA
               permanova_result <- vegan::adonis2(dist_matrix ~ Group, 
                                                  data = metadata, 
                                                  permutations = 999)
               
               return(list(
                 test_type = "PERMANOVA",
                 results = permanova_result,
                 summary = paste("PERMANOVA R² =", round(permanova_result$R2[1], 4),
                                 ", p-value =", round(permanova_result$`Pr(>F)`[1], 4),
                                 ifelse(permanova_result$`Pr(>F)`[1] < 0.05, 
                                        "(Significant)", "(Not significant)"))
               ))
             },
             
             "DESeq2" = {
               # DESeq2 analysis
               if (!requireNamespace("DESeq2", quietly = TRUE)) {
                 return(list(
                   test_type = "DESeq2",
                   error = "DESeq2 package not installed. Please install from Bioconductor."
                 ))
               }
               
               # Prepare count matrix (ensure integers)
               count_matrix <- round(abundance_matrix)
               count_matrix <- count_matrix[, colSums(count_matrix) > 0]  # Remove empty columns
               
               # Create DESeq2 dataset
               sample_data <- data.frame(
                 Group = factor(metadata$Group, levels = c("Control", "Patient"))
               )
               
               dds <- DESeq2::DESeqDataSetFromMatrix(
                 countData = t(count_matrix),  # DESeq2 expects genes as rows
                 colData = sample_data,
                 design = ~ Group
               )
               
               # Filter low count features
               keep <- rowSums(DESeq2::counts(dds)) >= 10
               dds <- dds[keep,]
               
               # Run DESeq2 analysis
               dds <- DESeq2::DESeq(dds, quiet = TRUE)
               results_deseq <- DESeq2::results(dds, contrast = c("Group", "Patient", "Control"))
               
               # Convert to data frame and clean up
               results_df <- as.data.frame(results_deseq)
               results_df$Taxon <- rownames(results_df)
               results_df$Significant <- !is.na(results_df$padj) & results_df$padj < 0.05
               
               # Remove rows with NA p-values
               results_df <- results_df[!is.na(results_df$pvalue), ]
               
               return(list(
                 test_type = "DESeq2",
                 results = results_df,
                 summary = paste("Tested", nrow(results_df), "taxa.",
                                 "Significant differences found in", 
                                 sum(results_df$Significant, na.rm = TRUE), 
                                 "taxa (FDR < 0.05).")
               ))
             }
      )
    }, error = function(e) {
      return(list(
        test_type = input$stat_test,
        error = paste("Error in statistical test:", e$message)
      ))
    })
  })
  
  # Add this output to display statistical results
  output$statistical_results <- renderUI({
    if (input$stat_test == "None" || !input$compare_groups) {
      return(NULL)
    }
    
    stat_results <- perform_statistical_test()
    
    if (is.null(stat_results)) {
      return(div(
        class = "alert alert-warning",
        icon("exclamation-triangle"),
        " No statistical analysis available. Ensure you have both Patient and Control samples."
      ))
    }
    
    if (!is.null(stat_results$error)) {
      return(div(
        class = "alert alert-danger",
        icon("exclamation-circle"),
        " Error: ", stat_results$error
      ))
    }
    
    # Create results display based on test type
    if (stat_results$test_type %in% c("Wilcoxon Rank-Sum Test", "T-test (log-transformed)", "DESeq2")) {
      # For tests that produce per-taxon results
      results_df <- stat_results$results
      
      if (nrow(results_df) == 0) {
        return(div(
          class = "alert alert-info",
          icon("info-circle"),
          " No taxa met the criteria for statistical testing."
        ))
      }
      
      # Sort by significance and p-value
      results_df <- results_df[order(results_df$Significant, results_df$P_value, decreasing = c(TRUE, FALSE)), ]
      
      # Display top significant results
      significant_results <- results_df[results_df$Significant == TRUE, ]
      
      result_content <- div(
        h4(paste("Results:", stat_results$test_type), style = "color: #337ab7;"),
        p(stat_results$summary, style = "font-weight: bold;"),
        
        if (nrow(significant_results) > 0) {
          div(
            h5("Significantly Different Taxa (FDR < 0.05):", style = "color: #d9534f;"),
            renderTable({
              display_cols <- c("Taxon", "P_value", "P_adjusted")
              if ("Mean_Patients" %in% colnames(significant_results)) {
                display_cols <- c(display_cols, "Mean_Patients", "Mean_Controls")
              }
              if ("log2FoldChange" %in% colnames(significant_results)) {
                display_cols <- c("Taxon", "log2FoldChange", "pvalue", "padj")
                significant_results$pvalue <- significant_results$pvalue
                significant_results$padj <- significant_results$padj
              }
              
              significant_results[, intersect(display_cols, colnames(significant_results))]
            }, digits = 4, striped = TRUE, hover = TRUE)
          )
        } else {
          div(
            class = "alert alert-info",
            icon("info-circle"),
            " No taxa showed significant differences between groups."
          )
        }
      )
      
    } else if (stat_results$test_type == "PERMANOVA") {
      # For PERMANOVA results
      permanova_table <- stat_results$results
      
      result_content <- div(
        h4("Results: PERMANOVA", style = "color: #337ab7;"),
        p(stat_results$summary, style = "font-weight: bold;"),
        renderTable({
          # Format PERMANOVA results for display
          display_table <- data.frame(
            Source = rownames(permanova_table),
            Df = permanova_table$Df,
            SumOfSqs = round(permanova_table$SumOfSqs, 4),
            R2 = round(permanova_table$R2, 4),
            F_value = round(permanova_table$`F`, 4),
            P_value = round(permanova_table$`Pr(>F)`, 4)
          )
          display_table[!is.na(display_table$P_value), ]
        }, striped = TRUE, hover = TRUE)
      )
    }
    
    return(result_content)
  })
  
  # Add download handler for statistical results
  output$download_statistical_results <- downloadHandler(
    filename = function() {
      paste("statistical_results_", input$stat_test, "_", 
            format(Sys.time(), "%Y%m%d_%H%M%S"), ".csv", sep = "")
    },
    content = function(file) {
      stat_results <- perform_statistical_test()
      
      if (!is.null(stat_results) && !is.null(stat_results$results) && is.null(stat_results$error)) {
        if (input$stat_test == "PERMANOVA") {
          # For PERMANOVA, save the results table
          permanova_df <- data.frame(
            Source = rownames(stat_results$results),
            stat_results$results,
            stringsAsFactors = FALSE
          )
          write.csv(permanova_df, file, row.names = FALSE)
        } else {
          # For other tests, save the results dataframe
          write.csv(stat_results$results, file, row.names = FALSE)
        }
      }
    }
  )
  
  
  ## HEATMAP 
  # Heatmap color schemes for ComplexHeatmap
  heatmap_colors <- list(
    "Viridis" = viridis::viridis(100),
    "Magma" = viridis::magma(100),
    "Plasma" = viridis::plasma(100),
    "RdBu" = rev(RColorBrewer::brewer.pal(11, "RdBu")), 
    "Beyonce" = beyonce::beyonce_palette(76), 
    "TayloRSwift" = tayloRswift::swift_palettes$speakNow
  )
  
  # Prepare data matrix for heatmap 
  heatmap_matrix_data <- reactive({
    req(taxonomy_with_groups())
    
    # Get taxonomy data with groups
    plot_data <- taxonomy_with_groups()
    
    # Aggregate data by Sample_ID and TaxonomicGroup
    agg_data <- aggregate(RelativeAbundance ~ Sample_ID + TaxonomicGroup, 
                          data = plot_data, FUN = mean)
    
    # Convert to wide format (matrix)
    heatmap_matrix <- dcast(agg_data, TaxonomicGroup ~ Sample_ID, 
                            value.var = "RelativeAbundance", fill = 0)
    
    # Set row names and remove the TaxonomicGroup column
    rownames(heatmap_matrix) <- heatmap_matrix$TaxonomicGroup
    heatmap_matrix <- heatmap_matrix[, -1, drop = FALSE]
    
    # Convert to matrix
    heatmap_matrix <- as.matrix(heatmap_matrix)
    
    # Filter out taxa with very low abundance across all samples
    # Keep only taxa that have at least the minimum abundance threshold in at least one sample
    keep_taxa <- rowSums(heatmap_matrix >= input$abundance_threshold) > 0
    heatmap_matrix <- heatmap_matrix[keep_taxa, , drop = FALSE]
    
    # Limit to top taxa if there are too many (for better visualization)
    if (nrow(heatmap_matrix) > 50) {
      # Calculate row means and keep top 50 most abundant taxa
      row_means <- rowMeans(heatmap_matrix, na.rm = TRUE)
      top_taxa <- names(sort(row_means, decreasing = TRUE))[1:50]
      heatmap_matrix <- heatmap_matrix[top_taxa, , drop = FALSE]
    }
    
    # Log transform the data for better visualization (add pseudocount to avoid log(0))
    heatmap_matrix <- log10(heatmap_matrix + 0.01)
    
    return(heatmap_matrix)
  })
  
  # Prepare column annotations for heatmap
  heatmap_annotations <- reactive({
    req(heatmap_matrix_data(), merged_metadata())
    
    # Get the matrix data
    matrix_data <- heatmap_matrix_data()
    sample_ids <- colnames(matrix_data)
    
    # Get metadata for these samples
    metadata <- merged_metadata()
    sample_metadata <- metadata[metadata$Sample_ID %in% sample_ids, ]
    
    # Ensure the order matches the matrix columns
    sample_metadata <- sample_metadata[match(sample_ids, sample_metadata$Sample_ID), ]
    
    # Create annotation data frame
    annotation_df <- data.frame(
      Sample_ID = sample_metadata$Sample_ID,
      Group = sample_metadata$Group,
      stringsAsFactors = FALSE
    )
    
    # Add metadata grouping if selected
    if (input$metadata_group != "None" && input$metadata_group %in% colnames(sample_metadata)) {
      annotation_df[[input$metadata_group]] <- sample_metadata[[input$metadata_group]]
    }
    
    # Set row names to match matrix columns
    rownames(annotation_df) <- annotation_df$Sample_ID
    annotation_df <- annotation_df[, -1, drop = FALSE]  # Remove Sample_ID column
    
    return(annotation_df)
  })
  
  # Generate the heatmap plot
  output$heatmap_plot <- renderPlot({
    req(heatmap_matrix_data(), heatmap_annotations())
    
    # Get data
    matrix_data <- heatmap_matrix_data()
    annotations <- heatmap_annotations()
    
    # Get selected color scheme
    color_scheme <- input$heatmap_color
    if (color_scheme %in% names(heatmap_colors)) {
      heatmap_col <- heatmap_colors[[color_scheme]]
    } else {
      heatmap_col <- heatmap_colors[["Viridis"]]
    }
    
    # Create color function
    if (color_scheme == "RdBu") {
      # For RdBu, center around 0 (since we're using log-transformed data)
      col_fun <- colorRamp2(c(min(matrix_data, na.rm = TRUE), 
                              mean(range(matrix_data, na.rm = TRUE)), 
                              max(matrix_data, na.rm = TRUE)), 
                            c("blue", "white", "red"))
    } else {
      col_fun <- colorRamp2(seq(min(matrix_data, na.rm = TRUE), 
                                max(matrix_data, na.rm = TRUE), 
                                length.out = length(heatmap_col)), 
                            heatmap_col)
    }
    
    # Create column annotation colors
    annotation_colors <- list()
    
    # Group colors
    if ("Group" %in% colnames(annotations)) {
      group_levels <- unique(annotations$Group)
      annotation_colors$Group <- setNames(c("#E31A1C", "#1F78B4")[1:length(group_levels)], 
                                          group_levels)
    }
    
    # Metadata group colors (if applicable)
    if (input$metadata_group != "None" && input$metadata_group %in% colnames(annotations)) {
      metadata_levels <- unique(annotations[[input$metadata_group]])
      if (is.numeric(metadata_levels)) {
        # For numeric metadata, use a continuous color scale
        annotation_colors[[input$metadata_group]] <- colorRamp2(
          range(metadata_levels, na.rm = TRUE),
          c("lightblue", "darkblue")
        )
      } else {
        # For categorical metadata, use discrete colors
        n_levels <- length(metadata_levels)
        if (n_levels <= 8) {
          colors <- RColorBrewer::brewer.pal(max(3, n_levels), "Set2")[1:n_levels]
        } else {
          colors <- rainbow(n_levels)
        }
        annotation_colors[[input$metadata_group]] <- setNames(colors, metadata_levels)
      }
    }
    
    # Create column annotation
    col_annotation <- HeatmapAnnotation(
      df = annotations,
      col = annotation_colors,
      annotation_name_gp = gpar(fontsize = 10),
      annotation_legend_param = list(
        Group = list(title = "Group", title_gp = gpar(fontsize = 10)),
        labels_gp = gpar(fontsize = 8)
      )
    )
    
    # Create the heatmap
    ht <- Heatmap(
      matrix_data,
      name = "Log10(Abundance + 0.01)",
      col = col_fun,
      
      # Row parameters
      row_title = paste(input$tax_level, "Taxa"),
      row_title_gp = gpar(fontsize = 12, fontface = "bold"),
      row_names_gp = gpar(fontsize = 8),
      row_names_max_width = unit(4, "cm"),
      show_row_dend = TRUE,
      row_dend_width = unit(2, "cm"),
      
      # Column parameters
      column_title = "Samples",
      column_title_gp = gpar(fontsize = 12, fontface = "bold"),
      column_names_gp = gpar(fontsize = 6),
      column_names_rot = 90,
      show_column_dend = TRUE,
      column_dend_height = unit(2, "cm"),
      
      # Clustering
      clustering_distance_rows = "euclidean",
      clustering_method_rows = "complete",
      clustering_distance_columns = "euclidean",
      clustering_method_columns = "complete",
      
      # Annotations
      top_annotation = col_annotation,
      
      # Legend
      heatmap_legend_param = list(
        title = "Relative\nAbundance\n(Log10)",
        title_gp = gpar(fontsize = 10),
        labels_gp = gpar(fontsize = 8),
        legend_direction = "vertical",
        legend_width = unit(4, "cm")
      ),
      
      # Additional parameters
      border = TRUE,
      rect_gp = gpar(col = "white", lwd = 0.5)
    )
    
    # Draw the heatmap
    draw(ht, 
         heatmap_legend_side = "right",
         annotation_legend_side = "right",
         merge_legend = TRUE)
    
  }, height = 600, width = 1000)
  
  # Download handler for the heatmap
  output$download_heatmap <- downloadHandler(
    filename = function() {
      paste("taxonomy_heatmap_", input$tax_level, "_", 
            format(Sys.time(), "%Y%m%d_%H%M%S"), ".png", sep = "")
    },
    content = function(file) {
      # Get data (same as in the render function)
      matrix_data <- heatmap_matrix_data()
      annotations <- heatmap_annotations()
      
      # Get selected color scheme
      color_scheme <- input$heatmap_color
      if (color_scheme %in% names(heatmap_colors)) {
        heatmap_col <- heatmap_colors[[color_scheme]]
      } else {
        heatmap_col <- heatmap_colors[["Viridis"]]
      }
      
      # Create color function
      if (color_scheme == "RdBu") {
        col_fun <- colorRamp2(c(min(matrix_data, na.rm = TRUE), 
                                mean(range(matrix_data, na.rm = TRUE)), 
                                max(matrix_data, na.rm = TRUE)), 
                              c("blue", "white", "red"))
      } else {
        col_fun <- colorRamp2(seq(min(matrix_data, na.rm = TRUE), 
                                  max(matrix_data, na.rm = TRUE), 
                                  length.out = length(heatmap_col)), 
                              heatmap_col)
      }
      
      # Create annotation colors (same logic as above)
      annotation_colors <- list()
      
      if ("Group" %in% colnames(annotations)) {
        group_levels <- unique(annotations$Group)
        annotation_colors$Group <- setNames(c("#E31A1C", "#1F78B4")[1:length(group_levels)], 
                                            group_levels)
      }
      
      if (input$metadata_group != "None" && input$metadata_group %in% colnames(annotations)) {
        metadata_levels <- unique(annotations[[input$metadata_group]])
        if (is.numeric(metadata_levels)) {
          annotation_colors[[input$metadata_group]] <- colorRamp2(
            range(metadata_levels, na.rm = TRUE),
            c("lightblue", "darkblue")
          )
        } else {
          n_levels <- length(metadata_levels)
          if (n_levels <= 8) {
            colors <- RColorBrewer::brewer.pal(max(3, n_levels), "Set2")[1:n_levels]
          } else {
            colors <- rainbow(n_levels)
          }
          annotation_colors[[input$metadata_group]] <- setNames(colors, metadata_levels)
        }
      }
      
      # Create column annotation
      col_annotation <- HeatmapAnnotation(
        df = annotations,
        col = annotation_colors,
        annotation_name_gp = gpar(fontsize = 10)
      )
      
      # Create the heatmap
      ht <- Heatmap(
        matrix_data,
        name = "Log10(Abundance + 0.01)",
        col = col_fun,
        row_title = paste(input$tax_level, "Taxa"),
        row_title_gp = gpar(fontsize = 12, fontface = "bold"),
        row_names_gp = gpar(fontsize = 8),
        column_title = "Samples",
        column_title_gp = gpar(fontsize = 12, fontface = "bold"),
        column_names_gp = gpar(fontsize = 6),
        column_names_rot = 90,
        top_annotation = col_annotation,
        show_row_dend = TRUE,
        show_column_dend = TRUE,
        border = TRUE
      )
      
      # Save the heatmap
      png(file, width = 12, height = 10, units = "in", res = 300)
      draw(ht, heatmap_legend_side = "right", annotation_legend_side = "right")
      dev.off()
    }
  )
  
  ###### DIVERSITY ANALYSIS TAB ######
  rv <- reactiveValues(latest_ggplot = NULL)
  
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
  
  
  # Add this reactive values object at the top of your server function
  rv <- reactiveValues(latest_ggplot = NULL)
  
  # Modified diversity_main_plot with ggplot storage
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
    
    # STORE THE GGPLOT OBJECT FOR DOWNLOAD
    rv$latest_ggplot <- p
    rv$plot_type <- "alpha"
    
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
      filter(!is.na(Sample_ID), !is.na(Species), !is.na(Abundance)) %>%
      pivot_wider(names_from = Species, values_from = Abundance, values_fill = 0)
    
    if (nrow(sample_wide) == 0) {
      cat("❌ No valid sample data after cleaning.\n")
      return(NULL)
    }
    
    sample_div <- sample_wide %>%
      column_to_rownames("Sample_ID") %>%
      as.matrix() %>%
      {
        matrix_data <- .
        matrix_data[is.na(matrix_data)] <- 0
        tibble(
          Sample_ID = rownames(matrix_data),
          Observed = rowSums(matrix_data > 0, na.rm = TRUE),
          Shannon = vegan::diversity(matrix_data, index = "shannon"),
          Simpson = vegan::diversity(matrix_data, index = "simpson")
        )
      }
    
    if (input$subdivide_samples) {
      req(input$condition_column)
      selected_col <- input$condition_column
      
      if (!selected_col %in% names(sample_meta)) {
        cat("❌ Selected condition column not found in sample metadata.\n")
        return(NULL)
      }
      
      sample_div <- left_join(sample_div,
                              sample_meta %>% select(Sample_ID, !!sym(selected_col)),
                              by = "Sample_ID") %>%
        filter(!is.na(.data[[selected_col]]), .data[[selected_col]] != "")
      sample_div$GroupLabel <- as.character(sample_div[[selected_col]])
    } else {
      sample_div$GroupLabel <- "Sample"
    }
    
    # Control diversity
    control_wide <- control_tax %>%
      select(Sample_ID, Species, Abundance) %>%
      filter(!is.na(Sample_ID), !is.na(Species), !is.na(Abundance)) %>%
      pivot_wider(names_from = Species, values_from = Abundance, values_fill = 0)
    
    if (nrow(control_wide) == 0) {
      cat("❌ No valid control data after cleaning.\n")
      return(NULL)
    }
    
    control_div <- control_wide %>%
      column_to_rownames("Sample_ID") %>%
      as.matrix() %>%
      {
        matrix_data <- .
        matrix_data[is.na(matrix_data)] <- 0
        tibble(
          Sample_ID = rownames(matrix_data),
          Observed = rowSums(matrix_data > 0, na.rm = TRUE),
          Shannon = vegan::diversity(matrix_data, index = "shannon"),
          Simpson = vegan::diversity(matrix_data, index = "simpson"),
          GroupLabel = "Control"
        )
      }
    
    diversity_df <- bind_rows(sample_div, control_div)
    
    metric_col <- switch(input$alpha_metric,
                         "Observed OTUs" = "Observed",
                         "Shannon" = "Shannon",
                         "Simpson" = "Simpson")
    
    # Filter invalid entries
    diversity_df <- diversity_df %>%
      filter(
        !is.na(GroupLabel),
        GroupLabel != "",
        !is.na(.data[[metric_col]]),
        is.finite(.data[[metric_col]])
      )
    
    if (nrow(diversity_df) == 0) {
      cat("❌ No valid data remaining after filtering.\n")
      return(NULL)
    }
    
    # Drop groups with <2 samples
    group_counts <- table(diversity_df$GroupLabel)
    valid_groups <- names(group_counts[group_counts >= 2])
    diversity_df <- diversity_df %>% filter(GroupLabel %in% valid_groups)
    diversity_df$GroupLabel <- factor(diversity_df$GroupLabel)
    
    n_groups <- n_distinct(diversity_df$GroupLabel)
    
    if (n_groups < 2) {
      cat("❌ Not enough valid groups to compare after filtering groups with <2 samples.\n")
      cat("Available groups:", paste(unique(diversity_df$GroupLabel), collapse = ", "), "\n")
      return(NULL)
    }
    
    cat("\n")
    
    if (n_groups == 2) {
      cat("✅ Two valid groups detected — using Wilcoxon rank-sum test:\n")
      tryCatch({
        wilcox_result <- wilcox.test(as.formula(paste(metric_col, "~ GroupLabel")), data = diversity_df)
        print(wilcox_result)
      }, error = function(e) {
        cat("❌ Error in Wilcoxon test:", e$message, "\n")
      })
    } else {
      cat("🔍 Performing pairwise Wilcoxon comparisons with BH correction:\n")
      tryCatch({
        pairwise_res <- pairwise.wilcox.test(diversity_df[[metric_col]],
                                             diversity_df$GroupLabel,
                                             p.adjust.method = "BH")
        print(pairwise_res)
      }, error = function(e) {
        cat("❌ Error in pairwise Wilcoxon test:", e$message, "\n")
      })
      
    }
  })
  
  output$ordination_plot <- renderPlotly({
    req(input$diversity_type == "beta")
    req(data_store$sample_metadata, data_store$taxonomy_data)
    req(data_store$control_sample_metadata, data_store$control_taxonomy_data)
    req(input$distance_metric, input$ordination_method)
    
    # Combine taxonomy
    all_tax <- bind_rows(
      data_store$taxonomy_data,
      data_store$control_taxonomy_data
    ) %>%
      select(Sample_ID, Species, Abundance) %>%
      pivot_wider(names_from = Species, values_from = Abundance, values_fill = 0) %>%
      column_to_rownames("Sample_ID")
    
    # Distance
    dist_method <- switch(input$distance_metric,
                          "Bray-Curtis" = "bray",
                          "Euclidean"   = "euclidean",
                          "Jaccard"     = "jaccard",
                          "Canberra"    = "canberra",
                          "Manhattan"   = "manhattan",
                          "Kulczynski"  = "kulczynski",
                          "Chord"       = "chord")
    dist_matrix <- vegan::vegdist(all_tax, method = dist_method)
    
    # Ordination
    coords <- NULL
    title_prefix <- ""
    
    if (input$ordination_method == "PCoA") {
      ord <- ape::pcoa(dist_matrix)
      coords <- ord$vectors[, 1:2]
      title_prefix <- "PCoA"
    } else if (input$ordination_method == "NMDS") {
      ord <- vegan::metaMDS(dist_matrix, k = 2, trymax = 100, autotransform = FALSE, trace = FALSE)
      coords <- ord$points
      title_prefix <- "NMDS"
    } else if (input$ordination_method == "dbRDA") {
      ord <- vegan::capscale(dist_matrix ~ 1)
      coords <- scores(ord, display = "sites")[, 1:2]
      title_prefix <- "dbRDA"
    } else if (input$ordination_method == "t-SNE") {
      if (!requireNamespace("Rtsne", quietly = TRUE)) {
        showNotification("Install 'Rtsne' to use t-SNE", type = "error")
        return(NULL)
      }
      dist_mat <- as.matrix(dist_matrix)
      perplexity_val <- min(30, floor((nrow(dist_mat) - 1) / 3))
      if (perplexity_val < 5) {
        showNotification("Too few samples for t-SNE.", type = "error")
        return(NULL)
      }
      tsne_result <- tryCatch({
        Rtsne::Rtsne(dist_mat, is_distance = TRUE, perplexity = perplexity_val)
      }, error = function(e) {
        showNotification(paste("t-SNE error:", e$message), type = "error")
        return(NULL)
      })
      if (is.null(tsne_result) || is.null(tsne_result$Y)) return(NULL)
      coords <- tsne_result$Y
      rownames(coords) <- rownames(dist_mat)
      title_prefix <- "t-SNE"
    } else if (input$ordination_method == "UMAP") {
      if (!requireNamespace("umap", quietly = TRUE)) {
        showNotification("Install 'umap' to use UMAP", type = "error")
        return(NULL)
      }
      umap_result <- tryCatch({
        umap::umap(as.matrix(all_tax))
      }, error = function(e) {
        showNotification(paste("UMAP error:", e$message), type = "error")
        return(NULL)
      })
      if (is.null(umap_result) || is.null(umap_result$layout)) return(NULL)
      coords <- umap_result$layout[, 1:2]
      rownames(coords) <- rownames(all_tax)
      title_prefix <- "UMAP"
    }
    
    if (is.null(coords)) {
      showNotification("Ordination failed.", type = "error")
      return(NULL)
    }
    
    # Build plotting data
    ord_df <- as.data.frame(coords)
    colnames(ord_df) <- c("Axis1", "Axis2")
    ord_df$Sample_ID <- rownames(coords)
    all_meta <- bind_rows(data_store$sample_metadata, data_store$control_sample_metadata)
    ord_df <- left_join(ord_df, all_meta, by = "Sample_ID")
    
    if (input$subdivide_samples) {
      req(input$condition_column)
      ord_df$GroupLabel <- ifelse(
        ord_df$Sample_ID %in% data_store$control_sample_metadata$Sample_ID,
        "Control",
        ord_df[[input$condition_column]]
      )
    } else {
      ord_df$GroupLabel <- ifelse(ord_df$Sample_ID %in% data_store$control_sample_metadata$Sample_ID,
                                  "Control", "Sample")
    }
    
    # ggplot
    p <- ggplot(ord_df, aes(x = Axis1, y = Axis2, color = GroupLabel, text = Sample_ID)) +
      geom_point(size = 3, alpha = 0.8)
    
    if (input$show_group_ellipses) {
      p <- p + stat_ellipse(aes(group = GroupLabel, color = GroupLabel),
                            type = "norm", linetype = "dashed", alpha = 0.6, size = 1)
    }
    
    if (!is.null(input$selected_sample_id) && input$selected_sample_id %in% ord_df$Sample_ID) {
      p <- p + geom_point(data = ord_df %>% filter(Sample_ID == input$selected_sample_id),
                          aes(x = Axis1, y = Axis2),
                          shape = 21, fill = "yellow", color = "red", size = 5, stroke = 1.5)
    }
    
    p <- p +
      labs(title = paste0(title_prefix, " on ", input$distance_metric, " Distance"),
           x = "Axis 1", y = "Axis 2") +
      theme_minimal()
    
    if (!input$show_diversity_legend) {
      p <- p + theme(legend.position = "none")
    }
    
    rv$latest_ggplot <- p
    rv$plot_type <- "beta"
    
    ggplotly(p, tooltip = c("x", "y", "text"))
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
    
    # Sanity check for Bray-Curtis
    if (input$distance_metric == "Bray-Curtis") {
      if (any(all_tax < 0, na.rm = TRUE)) {
        stop("Abundance matrix contains negative values — Bray-Curtis is not defined.")
      }
    }
    
    # Combine metadata
    all_meta <- bind_rows(data_store$sample_metadata, data_store$control_sample_metadata)
    if (input$subdivide_samples) {
      req(input$condition_column)
      group_label <- ifelse(
        all_meta$Sample_ID %in% data_store$control_sample_metadata$Sample_ID,
        "Control",
        all_meta[[input$condition_column]]
      )
    } else {
      group_label <- ifelse(all_meta$Sample_ID %in% data_store$control_sample_metadata$Sample_ID, "Control", "Sample")
    }
    
    # Distance matrix
    dist_method <- switch(input$distance_metric,
                          "Bray-Curtis" = "bray",
                          "Euclidean"   = "euclidean",
                          "Jaccard"     = "jaccard",
                          "Canberra"    = "canberra",
                          "Manhattan"   = "manhattan",
                          "Kulczynski"  = "kulczynski",
                          "Chord"       = "chord")
    
    
    dist_matrix <- vegan::vegdist(all_tax, method = dist_method)
    
    # Convert to data frame for adonis2
    stat_df <- data.frame(GroupLabel = group_label)
    rownames(stat_df) <- rownames(all_tax)
    
    # Run adonis2
    adonis_result <- vegan::adonis2(dist_matrix ~ GroupLabel, data = stat_df)
    print(adonis_result)
  })
  
  
  rv <- reactiveValues(latest_dendrogram = NULL)
  
  output$clustering_plot <- renderPlot({
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
    
    # Compute distance matrix
    dist_method <- switch(input$distance_metric,
                          "Bray-Curtis" = "bray",
                          "Euclidean"   = "euclidean",
                          "Jaccard"     = "jaccard",
                          "Canberra"    = "canberra",
                          "Manhattan"   = "manhattan",
                          "Kulczynski"  = "kulczynski",
                          "Chord"       = "chord")
    
    dist_matrix <- vegan::vegdist(all_tax, method = dist_method)
    hc <- hclust(dist_matrix, method = "ward.D2")
    dend <- as.dendrogram(hc)
    
    selected_sample <- input$selected_sample_id
    
    # Only change the color of the selected label
    dend <- dendrapply(dend, function(n) {
      if (is.leaf(n)) {
        label <- attr(n, "label")
        attr(n, "nodePar") <- list(
          lab.col = if (!is.null(selected_sample) && label == selected_sample) "red" else "black",
          lab.cex = 0.6
        )
      }
      return(n)
    })
    
    rv$latest_dendrogram <- dend
    
    # Plot dendrogram
    plot(dend,
         main = paste("Clustering Dendrogram (", input$distance_metric, ")"),
         ylab = "Distance")
  })
  

  
  output$download_alpha_plot <- downloadHandler(
    filename = function() {
      metric <- tolower(gsub(" ", "_", input$alpha_metric))
      paste0("diversity_alpha_", metric, "_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".png")
    },
    content = function(file) {
      if (is.null(rv$latest_ggplot)) {
        showNotification("No plot available for download. Please generate a plot first.", type = "error")
        return(NULL)
      }
      
      p <- rv$latest_ggplot +
        theme_minimal() +
        theme(
          legend.position = "right",
          plot.title = element_text(size = 14, hjust = 0.5),
          axis.title = element_text(size = 12),
          axis.text = element_text(size = 10),
          legend.text = element_text(size = 10),
          legend.title = element_text(size = 11)
        )
      
      ggsave(file, plot = p, width = 12, height = 8, dpi = 300, bg = "white")
    },
    contentType = "image/png"
  )
  
  output$download_beta_ordination_plot <- downloadHandler(
  filename = function() {
    method <- tolower(gsub("[^A-Za-z0-9]", "_", input$ordination_method))
    distance <- tolower(gsub("[^A-Za-z0-9]", "_", input$distance_metric))
    paste0("ordination_", method, "_", distance, "_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".png")
  },
  content = function(file) {
    if (is.null(rv$latest_ggplot)) {
      showNotification("No ordination plot available for download.", type = "error")
      return(NULL)
    }

    ggsave(file, plot = rv$latest_ggplot, width = 12, height = 8, dpi = 300, bg = "white")
  },
  contentType = "image/png"
)

  
  output$download_beta_dendrogram_plot <- downloadHandler(
    filename = function() {
      distance <- tolower(gsub("[^A-Za-z0-9]", "_", input$distance_metric))
      paste0("dendrogram_", distance, "_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".png")
    },
    content = function(file) {
      dend <- rv$latest_dendrogram
      if (is.null(dend)) {
        showNotification("No dendrogram available.", type = "error")
        return(NULL)
      }
      
      png(file, width = 1600, height = 1000, res = 150)
      plot(dend,
           main = paste("Clustering Dendrogram (", input$distance_metric, ")"),
           ylab = "Distance")
      dev.off()
    },
    contentType = "image/png"
  )
  

  
  ###### LIFESTYLE TAB ######
  output$parallel_var_select <- renderUI({
    req(data_store$sample_metadata)
    
    var_names <- names(data_store$sample_metadata)
    
    # Create named list: "Ongoing conditions" = "Ongoing_conditions", etc.
    named_vars <- setNames(var_names, gsub("_", " ", var_names))
    
    selectizeInput("parallel_vars", "Select Metadata Variables:",
                   choices = named_vars,
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
  
  output$download_parallel_plot <- downloadHandler(
    filename = function() {
      paste0("parallel_coordinates_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".png")
    },
    content = function(file) {
      req(data_store$sample_metadata, data_store$taxonomy_data, input$parallel_vars)
      
      # Total Abundance per sample
      taxa_data <- data_store$taxonomy_data %>%
        group_by(Sample_ID) %>%
        summarize(Total_Abundance = sum(Abundance), .groups = "drop")
      
      combined <- merge(data_store$sample_metadata, taxa_data, by = "Sample_ID")
      vars <- input$parallel_vars
      all_vars <- c(vars, "Total_Abundance")
      
      # Select, clean and prepare data
      plot_data <- combined %>%
        select(all_of(all_vars)) %>%
        na.omit() %>%
        mutate(across(where(is.character), as.factor)) %>%
        mutate(across(where(is.factor), ~ as.numeric(as.factor(.)))) %>%
        mutate(id = row_number())  # For coloring
      
      # User-friendly labels
      friendly_labels <- gsub("_", " ", colnames(plot_data))
      
      # Create static plot
      p <- GGally::ggparcoord(
        data = plot_data,
        columns = 1:length(all_vars),
        groupColumn = "id",
        scale = "uniminmax",
        showPoints = FALSE,
        alphaLines = 0.3
      ) +
        theme_minimal(base_size = 12) +
        theme(
          plot.background = element_rect(fill = "white", color = NA),
          panel.background = element_rect(fill = "white", color = NA),
          panel.grid.major = element_line(color = "gray90"),
          panel.grid.minor = element_blank()
        ) +
        labs(title = "Parallel Coordinates Plot", x = "Variables", y = "Scaled Value") +
        scale_x_discrete(labels = friendly_labels[1:length(all_vars)])
      
      
      # Save as PNG
      ggsave(filename = file, plot = p, width = 10, height = 6, dpi = 300, bg = "white")
      },
    contentType = "image/png"
  )
  
  
  
  
  
  # ==============================================================================
  # SERVER - PHYLOGENETIC TREE TAB - FIXED VERSION WITH BOTTOM LEGEND
  # ==============================================================================
  
  # Enhanced function to resolve phylo conflicts
  resolve_phylo_conflicts <- function() {
    # Clear any existing phylo class definitions from global environment
    if (exists("phylo", envir = .GlobalEnv)) {
      rm("phylo", envir = .GlobalEnv)
    }
    
    # Force detach and reattach ape to ensure clean phylo class
    if ("package:ape" %in% search()) {
      tryCatch({
        detach("package:ape", unload = TRUE, force = TRUE)
        library(ape, quietly = TRUE)
        cat("Reloaded ape package to resolve phylo conflicts\n")
      }, error = function(e) {
        cat("Could not reload ape package:", e$message, "\n")
      })
    }
    
    # Explicitly clear the S4 class cache
    if (exists(".S4_classes_cache", envir = .GlobalEnv)) {
      .S4_classes_cache <- new.env()
    }
    
    cat("Phylo class conflicts resolved\n")
  }
  
  resolve_phylo_conflicts()
  
  extract_node_support <- function(tree) {
    cat("=== NODE SUPPORT EXTRACTION ===\n")
    
    support_values <- NULL
    support_source <- NULL
    
    # Check multiple possible sources for support values
    if (!is.null(tree$node.label) && length(tree$node.label) > 0) {
      # Clean and filter node labels
      clean_labels <- tree$node.label
      clean_labels <- clean_labels[!is.na(clean_labels)]
      clean_labels <- clean_labels[clean_labels != ""]
      clean_labels <- clean_labels[clean_labels != "NA"]
      
      # Try to convert to numeric (for bootstrap values stored as strings)
      numeric_labels <- suppressWarnings(as.numeric(clean_labels))
      
      if (length(clean_labels) > 0) {
        cat("Found node.label with", length(clean_labels), "values\n")
        cat("Sample values:", head(clean_labels, 3), "\n")
        
        # If they're numeric, use them as support values
        if (any(!is.na(numeric_labels))) {
          support_values <- tree$node.label
          support_source <- "node.label"
          cat("Using node.label as support values\n")
        }
      }
    }
    
    # Check for explicit bootstrap/support attributes
    if (is.null(support_values)) {
      for (attr_name in c("bootstrap", "support", "posterior", "node.support")) {
        if (attr_name %in% names(tree) && !is.null(tree[[attr_name]])) {
          support_values <- tree[[attr_name]]
          support_source <- attr_name
          cat("Found support values in", attr_name, "\n")
          break
        }
      }
    }
    
    # Check edge attributes (sometimes support is stored there)
    if (is.null(support_values) && !is.null(tree$edge.length)) {
      # Sometimes support values are mixed with branch lengths
      if (any(tree$edge.length > 1 & tree$edge.length <= 100)) {
        cat("Potential support values found in edge.length (values > 1)\n")
      }
    }
    
    cat("Final result - Source:", support_source, "Values found:", !is.null(support_values), "\n")
    cat("===============================\n")
    
    return(list(values = support_values, source = support_source))
  }
  
  # Define the tree file path
  PHYLO_TREE_PATH <- "../group2_B/results/qiime_output/relevant_results/phylogenetic_tree.nwk"
  TAXONOMY_PATH   <- "../group2_B/results/qiime_output/relevant_results/taxonomy.tsv"
  
  # Reactive value to store tree data
  phylo_tree_data <- reactiveValues(
    tree = NULL,
    loaded = FALSE,
    file_exists = FALSE,
    error_msg = NULL,
    debug_info = ""
  )
  
  # Enhanced tree cleaning function
  clean_phylo_tree <- function(tree) {
    # Ensure it's a proper ape phylo object
    if (!inherits(tree, "phylo")) {
      tree <- ape::as.phylo(tree)
    }
    
    # Force the class to be exactly what ape expects
    class(tree) <- "phylo"
    
    # Ensure edge lengths are numeric and positive
    if (is.null(tree$edge.length)) {
      tree$edge.length <- rep(1, nrow(tree$edge))
    } else {
      # Convert to numeric and handle any issues
      tree$edge.length <- as.numeric(tree$edge.length)
      tree$edge.length[is.na(tree$edge.length)] <- 1
      tree$edge.length[tree$edge.length <= 0] <- 1
    }
    
    # Ensure edge matrix is integer matrix
    if (!is.null(tree$edge)) {
      tree$edge <- matrix(as.integer(tree$edge), ncol = 2)
    }
    
    # Ensure tip labels are character
    if (!is.null(tree$tip.label)) {
      tree$tip.label <- as.character(tree$tip.label)
    }
    
    # Validate tree structure
    if (!ape::is.rooted(tree)) {
      # For unrooted trees, ensure proper structure
      n_tips <- ape::Ntip(tree)
      n_nodes <- ape::Nnode(tree)
      expected_edges <- 2 * n_tips - 3
      
      if (nrow(tree$edge) != expected_edges) {
        cat("Warning: Tree structure may be invalid for unrooted plotting\n")
      }
    }
    
    return(tree)
  }
  
  # Function to load tree automatically
  load_phylo_tree <- function() {
    tryCatch({
      if (file.exists(PHYLO_TREE_PATH)) {
        tree <- ape::read.tree(PHYLO_TREE_PATH)
        
        # Clean the tree immediately after loading
        tree <- clean_phylo_tree(tree)
        
        if (file.exists(TAXONOMY_PATH)) {
          tax_tbl <- read.delim(
            TAXONOMY_PATH,
            header = TRUE, sep = "\t",
            quote = "", comment.char = "",
            stringsAsFactors = FALSE
          )
          tax_vec  <- setNames(tax_tbl$Taxon, tax_tbl$Feature.ID)
          
          missing <- !(tree$tip.label %in% names(tax_vec))
          if (any(missing)) {
            cat("Warning:", sum(missing),
                "tip labels without taxonomy; keeping their IDs\n")
          }
          tree$tip.label <- ifelse(
            missing,
            tree$tip.label,
            tax_vec[tree$tip.label]
          )
        } else {
          cat("⚠ No taxonomy.tsv found at:", TAXONOMY_PATH, "\n")
        }
        
        phylo_tree_data$tree <- tree
        phylo_tree_data$loaded <- TRUE
        phylo_tree_data$file_exists <- TRUE
        phylo_tree_data$error_msg <- NULL
        phylo_tree_data$debug_info <- paste("Tree loaded with", ape::Ntip(tree), "tips")
      } else {
        phylo_tree_data$tree <- NULL
        phylo_tree_data$loaded <- FALSE
        phylo_tree_data$file_exists <- FALSE
        phylo_tree_data$error_msg <- paste("File not found:", PHYLO_TREE_PATH)
      }
    }, error = function(e) {
      phylo_tree_data$tree <- NULL
      phylo_tree_data$loaded <- FALSE
      phylo_tree_data$file_exists <- file.exists(PHYLO_TREE_PATH)
      phylo_tree_data$error_msg <- paste("Error loading tree:", e$message)
    })
  }
  
  # SIMPLIFIED AND ROBUST function to extract taxonomy
  extract_taxonomy_simple <- function(tip_labels, level) {
    cat("=== EXTRACT TAXONOMY DEBUG ===\n")
    cat("Input level:", level, "\n")
    cat("Number of labels:", length(tip_labels), "\n")
    
    if (level == "none" || is.null(tip_labels) || length(tip_labels) == 0) {
      cat("Returning original labels (level=none or empty input)\n")
      return(as.character(tip_labels))
    }
    
    # Show sample labels
    cat("Sample labels:\n")
    for(i in 1:min(3, length(tip_labels))) {
      cat(i, ":", tip_labels[i], "\n")
    }
    
    # Simple but effective strategy
    result <- switch(level,
                     "kingdom" = extract_by_pattern(tip_labels, c("k__", "D_0__"), 1),
                     "phylum" = extract_by_pattern(tip_labels, c("p__", "D_1__"), 2),
                     "class" = extract_by_pattern(tip_labels, c("c__", "D_2__"), 3),
                     "order" = extract_by_pattern(tip_labels, c("o__", "D_3__"), 4),
                     "family" = extract_by_pattern(tip_labels, c("f__", "D_4__"), 5),
                     "genus" = extract_by_pattern(tip_labels, c("g__", "D_5__"), 6),
                     as.character(tip_labels)  # default
    )
    
    cat("Unique groups found:", length(unique(result)), "\n")
    cat("Sample results:\n")
    for(i in 1:min(3, length(result))) {
      cat(tip_labels[i], " -> ", result[i], "\n")
    }
    cat("==============================\n")
    
    return(as.character(result))
  }
  
  # Helper function to extract by patterns
  extract_by_pattern <- function(labels, patterns, position) {
    result <- character(length(labels))
    
    for(i in seq_along(labels)) {
      label <- as.character(labels[i])
      extracted <- FALSE
      
      # Try with specific patterns
      for(pattern in patterns) {
        if(grepl(pattern, label)) {
          parts <- unlist(strsplit(label, ";"))
          matching_part <- parts[grepl(paste0("^", pattern), parts)]
          
          if(length(matching_part) > 0) {
            clean_name <- gsub(paste0("^", pattern), "", matching_part[1])
            clean_name <- trimws(clean_name)
            
            if(nchar(clean_name) > 0 && !grepl("unidentified|uncultured|unknown", clean_name, ignore.case = TRUE)) {
              result[i] <- clean_name
              extracted <- TRUE
              break
            }
          }
        }
      }
      
      # If nothing found with patterns, use position
      if(!extracted) {
        if(grepl(";", label)) {
          parts <- unlist(strsplit(label, ";"))
          if(length(parts) >= position && nchar(trimws(parts[position])) > 0) {
            # Clean any taxonomic prefix
            clean_part <- gsub("^[a-zA-Z]_*[0-9]*_*", "", parts[position])
            clean_part <- trimws(clean_part)
            if(nchar(clean_part) > 0) {
              result[i] <- clean_part
              extracted <- TRUE
            }
          }
        }
      }
      
      # Last resort: create generic group
      if(!extracted) {
        result[i] <- paste0("Unknown_Group_", i %% 10)
      }
    }
    
    return(result)
  }
  
  # SIMPLIFIED but EFFECTIVE tree collapse function
  collapse_tree_simple <- function(tree, tax_level, min_group_size = 2) {
    cat("=== COLLAPSE TREE DEBUG ===\n")
    cat("Original tree tips:", ape::Ntip(tree), "\n")
    cat("Tax level:", tax_level, "\n")
    cat("Min group size:", min_group_size, "\n")
    
    if(is.null(tree) || tax_level == "none") {
      cat("Returning original tree (null or none)\n")
      return(tree)
    }
    
    # Clean the tree before processing
    tree <- clean_phylo_tree(tree)
    
    # Extract groups
    tip_groups <- extract_taxonomy_simple(tree$tip.label, tax_level)
    
    if(length(tip_groups) != length(tree$tip.label)) {
      cat("ERROR: Group extraction failed\n")
      return(tree)
    }
    
    # Count groups
    group_counts <- table(tip_groups)
    cat("Group distribution:\n")
    print(head(sort(group_counts, decreasing = TRUE), 10))
    
    # Find groups to collapse
    groups_to_collapse <- names(group_counts)[group_counts >= min_group_size]
    cat("Groups meeting criteria:", length(groups_to_collapse), "\n")
    
    if(length(groups_to_collapse) == 0) {
      cat("No groups meet minimum size - returning original tree\n")
      return(tree)
    }
    
    # Create new tree
    new_tree <- tree
    tips_to_remove <- c()
    
    for(group_name in groups_to_collapse) {
      # Find all tips of this group
      group_indices <- which(tip_groups == group_name)
      cat("Group '", group_name, "' has", length(group_indices), "members\n")
      
      if(length(group_indices) >= min_group_size) {
        # Keep the first tip as representative
        representative <- group_indices[1]
        
        # Change its label
        new_tree$tip.label[representative] <- paste0(group_name, " (n=", length(group_indices), ")")
        
        # Mark the rest for removal
        if(length(group_indices) > 1) {
          tips_to_remove <- c(tips_to_remove, group_indices[-1])
        }
      }
    }
    
    cat("Tips to remove:", length(tips_to_remove), "\n")
    
    # Remove duplicate tips
    if(length(tips_to_remove) > 0) {
      # Verify we're not removing all tips
      if(length(tips_to_remove) >= length(new_tree$tip.label)) {
        cat("ERROR: Would remove all tips - returning original\n")
        return(tree)
      }
      
      # Remove duplicate tips
      new_tree <- ape::drop.tip(new_tree, tips_to_remove)
      # Clean the tree after modifications
      new_tree <- clean_phylo_tree(new_tree)
      cat("New tree has", ape::Ntip(new_tree), "tips\n")
    }
    
    cat("==========================\n")
    return(new_tree)
  }
  
  # Function to get fixed color palette (palette 36)
  get_color_palette <- function(n_colors) {
    
    if (requireNamespace("viridis", quietly = TRUE)) {
      # Use viridis discrete colors
      return(viridis::viridis(n_colors, option = "D"))
    } else {
      # Fallback: manual viridis-like colors if package not available
      viridis_colors <- c(
        "#440154", "#482677", "#3F4A8A", "#31678E", "#26838F", "#1F9D8A", 
        "#6CCE5A", "#B6DE2B", "#FEE825", "#FFEA46", "#FCFFA4", "#F0F921"
      )
      
      if (n_colors <= length(viridis_colors)) {
        return(viridis_colors[1:n_colors])
      } else {
        return(colorRampPalette(viridis_colors)(n_colors))
      }
    }
  }
  
  # REACTIVE for processed tree
  processed_tree <- reactive({
    req(phylo_tree_data$loaded, phylo_tree_data$tree)
    
    cat("\n=== PROCESSING TREE ===\n")
    cat("Tax level input:", input$tax_level, "\n")
    
    tree <- phylo_tree_data$tree
    
    # If no taxonomic level selected, return original
    if(is.null(input$tax_level) || input$tax_level == "none") {
      cat("No taxonomic grouping requested\n")
      phylo_tree_data$debug_info <- "Original tree - no grouping"
      # Still clean the tree
      return(clean_phylo_tree(tree))
    }
    
    # Apply collapse
    min_size <- ifelse(is.null(input$min_group_size), 2, as.numeric(input$min_group_size))
    
    cat("Applying collapse with min_size =", min_size, "\n")
    
    collapsed_tree <- collapse_tree_simple(tree, input$tax_level, min_size)
    
    # Update debug info
    phylo_tree_data$debug_info <- paste(
      "Taxonomic level:", input$tax_level,
      "| Original tips:", ape::Ntip(tree),
      "| Final tips:", ape::Ntip(collapsed_tree),
      "| Reduction:", round((1 - ape::Ntip(collapsed_tree)/ape::Ntip(tree)) * 100, 1), "%"
    )
    
    cat("======================\n")
    return(collapsed_tree)
  })
  
  
  # Modified collapse_tree_simple function to respect the group counts option
  collapse_tree_simple <- function(tree, tax_level, min_group_size = 2, show_counts = TRUE) {
    cat("=== COLLAPSE TREE DEBUG ===\n")
    cat("Original tree tips:", ape::Ntip(tree), "\n")
    cat("Tax level:", tax_level, "\n")
    cat("Min group size:", min_group_size, "\n")
    cat("Show group counts:", show_counts, "\n")  # New debug line
    
    if(is.null(tree) || tax_level == "none") {
      cat("Returning original tree (null or none)\n")
      return(tree)
    }
    
    # Clean the tree before processing
    tree <- clean_phylo_tree(tree)
    
    # Extract groups
    tip_groups <- extract_taxonomy_simple(tree$tip.label, tax_level)
    
    if(length(tip_groups) != length(tree$tip.label)) {
      cat("ERROR: Group extraction failed\n")
      return(tree)
    }
    
    # Count groups
    group_counts <- table(tip_groups)
    cat("Group distribution:\n")
    print(head(sort(group_counts, decreasing = TRUE), 10))
    
    # Find groups to collapse
    groups_to_collapse <- names(group_counts)[group_counts >= min_group_size]
    cat("Groups meeting criteria:", length(groups_to_collapse), "\n")
    
    if(length(groups_to_collapse) == 0) {
      cat("No groups meet minimum size - returning original tree\n")
      return(tree)
    }
    
    # Create new tree
    new_tree <- tree
    tips_to_remove <- c()
    
    for(group_name in groups_to_collapse) {
      # Find all tips of this group
      group_indices <- which(tip_groups == group_name)
      cat("Group '", group_name, "' has", length(group_indices), "members\n")
      
      if(length(group_indices) >= min_group_size) {
        # Keep the first tip as representative
        representative <- group_indices[1]
        
        # MODIFIED: Change label based on show_counts option
        if(show_counts) {
          new_tree$tip.label[representative] <- paste0(group_name, " (n=", length(group_indices), ")")
        } else {
          new_tree$tip.label[representative] <- group_name
        }
        
        # Mark the rest for removal
        if(length(group_indices) > 1) {
          tips_to_remove <- c(tips_to_remove, group_indices[-1])
        }
      }
    }
    
    cat("Tips to remove:", length(tips_to_remove), "\n")
    
    # Remove duplicate tips
    if(length(tips_to_remove) > 0) {
      # Verify we're not removing all tips
      if(length(tips_to_remove) >= length(new_tree$tip.label)) {
        cat("ERROR: Would remove all tips - returning original\n")
        return(tree)
      }
      
      # Remove duplicate tips
      new_tree <- ape::drop.tip(new_tree, tips_to_remove)
      # Clean the tree after modifications
      new_tree <- clean_phylo_tree(new_tree)
      cat("New tree has", ape::Ntip(new_tree), "tips\n")
    }
    
    cat("==========================\n")
    return(new_tree)
  }
  
  # Modified processed_tree reactive to pass the show_counts parameter
  processed_tree <- reactive({
    req(phylo_tree_data$loaded, phylo_tree_data$tree)
    
    cat("\n=== PROCESSING TREE ===\n")
    cat("Tax level input:", input$tax_level, "\n")
    cat("Midpoint root option:", input$midpoint_root %||% FALSE, "\n")
    cat("Show group counts:", input$show_group_counts %||% TRUE, "\n")  # New debug line
    
    tree <- phylo_tree_data$tree
    
    # Apply midpoint rooting if requested
    if (!is.null(input$midpoint_root) && input$midpoint_root) {
      cat("Applying midpoint rooting...\n")
      tree <- midpoint_root_tree(tree)
    }
    
    # If no taxonomic level selected, return tree (possibly rooted)
    if(is.null(input$tax_level) || input$tax_level == "none") {
      cat("No taxonomic grouping requested\n")
      if (!is.null(input$midpoint_root) && input$midpoint_root) {
        phylo_tree_data$debug_info <- "Tree with midpoint rooting - no grouping"
      } else {
        phylo_tree_data$debug_info <- "Original tree - no grouping"
      }
      return(clean_phylo_tree(tree))
    }
    
    # Apply collapse with group counts option
    min_size <- ifelse(is.null(input$min_group_size), 2, as.numeric(input$min_group_size))
    show_counts <- ifelse(is.null(input$show_group_counts), TRUE, input$show_group_counts)
    
    cat("Applying collapse with min_size =", min_size, "and show_counts =", show_counts, "\n")
    
    # MODIFIED: Pass the show_counts parameter
    collapsed_tree <- collapse_tree_simple(tree, input$tax_level, min_size, show_counts)
    
    # Update debug info
    rooting_status <- if (!is.null(input$midpoint_root) && input$midpoint_root) "midpoint-rooted" else "original"
    counts_status <- if (show_counts) "with counts" else "without counts"
    
    phylo_tree_data$debug_info <- paste(
      "Taxonomic level:", input$tax_level,
      "| Tree:", rooting_status,
      "| Labels:", counts_status,
      "| Original tips:", ape::Ntip(phylo_tree_data$tree),
      "| Final tips:", ape::Ntip(collapsed_tree),
      "| Reduction:", round((1 - ape::Ntip(collapsed_tree)/ape::Ntip(phylo_tree_data$tree)) * 100, 1), "%"
    )
    
    cat("======================\n")
    return(collapsed_tree)
  })
  
  # Add observer to debug the checkbox state
  observeEvent(input$show_group_counts, {
    cat("\n*** SHOW GROUP COUNTS CHANGED TO:", input$show_group_counts, "***\n")
  }, ignoreInit = TRUE)
  
  
  # Updated plotting section with fixed node support functionality
  output$main_phylo_plot <- renderPlot({
    
    tree <- processed_tree()
    
    cat("Plotting tree with", ape::Ntip(tree), "tips\n")
    
    tryCatch({
      
      # Get layout and parameters
      layout <- input$tree_layout %||% "rectangular"
      text_sz <- as.numeric(input$tree_text_size %||% 1.2)
      line_sz <- as.numeric(input$tree_line_size %||% 0.5)
      
      # Get display options
      show_tips <- "show_tip_labels" %in% input$display_options
      show_nodes <- "show_node_labels" %in% input$display_options
      show_branches <- "show_branch_length" %in% input$display_options
      
      cat("Display options - Tips:", show_tips, "Nodes:", show_nodes, "Branches:", show_branches, "\n")
      cat("Using layout:", layout, "with viridis palette\n")
      
      # Clean tree
      tree <- clean_phylo_tree(tree)
      
      # Check if ggtree is available
      if (requireNamespace("ggtree", quietly = TRUE)) {
        
        # Get taxonomic groups for coloring
        if (!is.null(input$tax_level) && input$tax_level != "none") {
          tip_groups <- extract_taxonomy_simple(tree$tip.label, input$tax_level)
          group_levels <- unique(tip_groups)
          
          # Get color palette (viridis)
          colors <- get_color_palette(length(group_levels))
          names(colors) <- group_levels
          
          # Create metadata data frame
          meta_df <- data.frame(
            label = tree$tip.label,
            grp = tip_groups,
            stringsAsFactors = FALSE
          )
        } else {
          # No grouping - use single color
          meta_df <- data.frame(
            label = tree$tip.label,
            grp = "All",
            stringsAsFactors = FALSE
          )
          colors <- c("All" = "#440154")  # Dark viridis color
        }
        
        # Create base plot with error handling for circular layout
        if (layout == "circular") {
          # Try circular layout with fallback
          tryCatch({
            p <- ggtree::ggtree(tree, layout = "circular", size = line_sz) %<+% meta_df
          }, error = function(e) {
            cat("Circular layout failed, using fan layout instead:", e$message, "\n")
            # Fallback to fan layout which is more stable
            layout <<- "fan"
            p <<- ggtree::ggtree(tree, layout = "fan", size = line_sz) %<+% meta_df
          })
        } else {
          p <- ggtree::ggtree(tree, layout = layout, size = line_sz) %<+% meta_df
        }
        
        # Add branch length labels if requested
        if (show_branches && !is.null(tree$edge.length)) {
          p <- p + ggtree::geom_text2(aes(subset = !isTip, label = round(branch.length, 3)),
                                      hjust = -0.1, vjust = -0.5, size = text_sz * 0.8, color = "gray50")
        }
        
        # FIXED CIRCULAR LAYOUT HANDLING
        if (layout == "circular" || layout == "fan") {
          
          # Calculate optimal text size based on number of tips (INCREASED SIZES)
          n_tips <- ape::Ntip(tree)
          optimal_text_size <- case_when(
            n_tips <= 20 ~ text_sz * 1.2,   # Increased from 0.9
            n_tips <= 50 ~ text_sz * 1.0,   # Increased from 0.7
            n_tips <= 100 ~ text_sz * 0.8,  # Increased from 0.5
            n_tips <= 200 ~ text_sz * 0.6,  # Increased from 0.4
            TRUE ~ text_sz * 0.5            # Increased from 0.3
          )
          
          tryCatch({
            # Create base circular plot
            p <- ggtree::ggtree(tree, layout = "circular", size = line_sz) %<+% meta_df
            
            # Add basic styling first
            p <- p +
              scale_color_manual(values = colors, name = input$tax_level %||% "Group") +
              ggtree::theme_tree() +
              theme(
                plot.margin = margin(60, 60, 80, 60),
                legend.position = "none",
                panel.background = element_rect(fill = "white", color = NA),
                plot.background = element_rect(fill = "white", color = NA)
              )
            
            # Add tip labels/points
            if (show_tips && n_tips <= 100) {
              tryCatch({
                p <- p + ggtree::geom_tiplab(
                  aes(color = grp),
                  size = optimal_text_size,
                  hjust = -0.1,
                  alpha = 0.8
                )
              }, error = function(e) {
                cat("Tip labels failed for circular plot, using points instead\n")
                p <<- p + ggtree::geom_tippoint(aes(color = grp), size = 2, alpha = 0.7)
              })
            } else if (show_tips) {
              p <- p + ggtree::geom_tippoint(aes(color = grp), size = 1.5, alpha = 0.7)
            }
            
            # SIMPLIFIED BRANCH LENGTH SOLUTION FOR CIRCULAR PLOTS
            if (show_branches && !is.null(tree$edge.length)) {
              cat("Adding branch lengths to circular plot...\n")
              
              # Method 1: Try the most reliable approach first
              tryCatch({
                # Only show branch lengths for major branches to avoid clutter
                edge_threshold <- quantile(tree$edge.length[tree$edge.length > 0], 0.8, na.rm = TRUE)
                
                p <- p + ggtree::geom_text2(
                  aes(
                    subset = (!isTip & !is.na(branch.length) & branch.length >= edge_threshold),
                    label = round(branch.length, 3)
                  ),
                  hjust = 0.5,
                  vjust = -0.3,
                  size = optimal_text_size * 0.6,    # Increased from 0.4
                  color = "darkblue",
                  fontface = "bold",
                  alpha = 0.9
                )
                cat("Branch lengths added successfully using geom_text2 with threshold\n")
                
              }, error = function(e1) {
                cat("geom_text2 failed:", e1$message, "\n")
                
                # Method 2: Fallback to showing branch length statistics
                tryCatch({
                  branch_stats <- paste(
                    "Branch lengths: min =", round(min(tree$edge.length[tree$edge.length > 0], na.rm = TRUE), 4),
                    ", max =", round(max(tree$edge.length, na.rm = TRUE), 4),
                    ", mean =", round(mean(tree$edge.length, na.rm = TRUE), 4)
                  )
                  
                  p <<- p + labs(
                    subtitle = branch_stats,
                    caption = "Individual branch lengths filtered for clarity in circular layout"
                  ) + theme(
                    plot.subtitle = element_text(size = 10, color = "darkblue", hjust = 0.5),
                    plot.caption = element_text(size = 9, color = "gray60", hjust = 0.5)
                  )
                  cat("Added branch length summary to plot labels\n")
                  
                }, error = function(e2) {
                  cat("Even summary failed:", e2$message, "\n")
                })
              })
            }
            
            # Add node support (simplified)
            if (show_nodes) {
              cat("Adding node support to circular plot...\n")
              
              support_info <- extract_node_support(tree)
              support_values <- support_info$values
              support_source <- support_info$source
              
              if (!is.null(support_values) && !is.null(support_source)) {
                tryCatch({
                  if (support_source == "node.label") {
                    p <- p + ggtree::geom_nodelab(
                      aes(label = label),
                      subset = (!is.na(label) & label != "" & label != "NA"),
                      hjust = 0.5,
                      vjust = 1.2,
                      size = optimal_text_size * 0.6,
                      color = "red",
                      fontface = "bold",
                      alpha = 0.8
                    )
                  } else {
                    p <- p + ggtree::geom_nodepoint(
                      color = "red", size = 2, alpha = 0.7
                    )
                  }
                }, error = function(e) {
                  p <<- p + ggtree::geom_nodepoint(color = "orange", size = 1.5, alpha = 0.6)
                })
              } else {
                p <- p + ggtree::geom_nodepoint(color = "orange", size = 1.5, alpha = 0.6)
              }
            }
            
            # Print the plot
            print(p)
            
          }, error = function(e) {
            cat("ggtree circular plot failed, using ape fallback:", e$message, "\n")
            
            # ENHANCED APE FALLBACK FOR CIRCULAR PLOTS
            par(mar = c(2, 2, 4, 2), bg = "white")
            
            # Calculate optimal label size for ape (INCREASED SIZES)
            optimal_text_size_ape <- case_when(
              n_tips <= 30 ~ text_sz * 1.0,   # Increased from 0.8
              n_tips <= 60 ~ text_sz * 0.8,   # Increased from 0.6
              n_tips <= 100 ~ text_sz * 0.6,  # Increased from 0.4
              TRUE ~ text_sz * 0.5            # Increased from 0.3
            )
            
            # Plot the tree first
            ape::plot.phylo(
              tree,
              type = "fan",
              cex = optimal_text_size_ape,
              main = paste("Circular Tree (", n_tips, " tips)"),
              show.tip.label = show_tips && n_tips <= 50,
              show.node.label = FALSE,  # We'll add these separately
              edge.color = "darkblue",
              edge.width = line_sz * 2,
              tip.color = "darkred"
            )
            
            # Add branch lengths - this WILL work with ape
            if (show_branches && !is.null(tree$edge.length)) {
              cat("Adding branch lengths using ape::edgelabels\n")
              
              tryCatch({
                # Show only longer branches to avoid clutter
                edge_threshold <- quantile(tree$edge.length[tree$edge.length > 0], 0.7, na.rm = TRUE)
                long_edges <- which(tree$edge.length >= edge_threshold)
                
                if (length(long_edges) > 0 && length(long_edges) <= 30) {
                  # Show individual branch lengths for reasonable number of edges
                  ape::edgelabels(
                    text = round(tree$edge.length[long_edges], 3),
                    edge = long_edges,
                    cex = optimal_text_size_ape * 0.7,    # Increased from 0.5
                    col = "darkblue",
                    bg = "lightblue",
                    frame = "rect"
                  )
                  cat("Individual branch lengths added successfully\n")
                } else {
                  # Show summary statistics
                  mtext(
                    paste("Branch lengths: mean =", round(mean(tree$edge.length, na.rm = TRUE), 4),
                          ", range =", round(min(tree$edge.length[tree$edge.length > 0], na.rm = TRUE), 4),
                          "to", round(max(tree$edge.length, na.rm = TRUE), 4)),
                    side = 3, line = 0, cex = optimal_text_size_ape * 0.8, col = "darkblue"
                  )
                  cat("Branch length summary added\n")
                }
                
              }, error = function(e) {
                cat("ape edgelabels failed:", e$message, "\n")
                # Add basic branch length info as title
                title(sub = paste("Mean branch length:", round(mean(tree$edge.length, na.rm = TRUE), 4)),
                      col.sub = "darkblue", cex.sub = 0.8)
              })
            }
            
            # Add node support if available
            if (show_nodes) {
              support_info <- extract_node_support(tree)
              if (!is.null(support_info$values) && support_info$source == "node.label") {
                tryCatch({
                  ape::nodelabels(
                    tree$node.label,
                    cex = optimal_text_size_ape * 0.6,
                    col = "red",
                    bg = "white",
                    frame = "circle"
                  )
                }, error = function(e) {
                  ape::nodelabels(pch = 19, cex = 0.8, col = "red")
                })
              }
            }
            
            return(NULL)
          })
          
        } else if (layout == "unrooted") {
          # IMPROVED UNROOTED LAYOUT
          p <- p +
            ggtree::geom_tippoint(
              aes(color = grp),
              size = 3, 
              stroke = 0,
              alpha = 0.8
            ) +
            scale_color_manual(values = colors, name = input$tax_level %||% "Group") +
            ggtree::theme_tree() +
            coord_equal(clip = "off") +
            theme(
              legend.position = "bottom",
              legend.title = element_text(size = 12),
              legend.text = element_text(size = 10),
              plot.margin = margin(40, 40, 60, 40),
              panel.background = element_rect(fill = "white", color = NA),
              plot.background = element_rect(fill = "white", color = NA)
            )
          
          # Add tip labels only for small unrooted trees
          if (show_tips && ape::Ntip(tree) <= 30) {
            p <- p + ggtree::geom_tiplab(
              size = text_sz * 0.6, 
              hjust = -0.1,
              alpha = 0.8
            )
          }
          
        } else {
          # RECTANGULAR LAYOUT (existing code with minor improvements)
          p <- p +
            scale_color_manual(values = colors, name = input$tax_level %||% "Group")
          
          # Add tip labels if requested
          if (show_tips) {
            p <- p +
              ggtree::geom_tiplab(
                aes(color = grp),
                size = text_sz,
                align = TRUE,
                linesize = 0.3,
                linetype = "dotted",
                offset = 0.002,
                fontface = "bold"
              ) +
              ggplot2::xlim(NA, max(p$data$x, na.rm = TRUE) * 1.3)
          }
          
          # Add styling
          p <- p +
            ggplot2::coord_cartesian(clip = "off") +
            ggtree::theme_tree() +
            theme(
              plot.margin = margin(20, 20, 80, 20),
              legend.position = "none",
              panel.background = element_rect(fill = "white", color = NA),
              plot.background = element_rect(fill = "white", color = NA)
            )
          
          # Add subtle branch coloring
          if (!is.null(input$tax_level) && input$tax_level != "none") {
            p <- p +
              ggtree::geom_tree(
                aes(color = grp),
                alpha = 0.3,
                size = line_sz * 0.8
              )
          }
        }
        
        # Add node support values (existing code - keeping as is)
        if (show_nodes) {
          cat("Attempting to show node support...\n")
          
          support_info <- extract_node_support(tree)
          support_values <- support_info$values
          support_source <- support_info$source
          
          if (!is.null(support_values) && !is.null(support_source)) {
            cat("Adding node support labels from:", support_source, "\n")
            
            tryCatch({
              if (support_source == "node.label") {
                # Convert to numeric if they look like numbers
                numeric_values <- suppressWarnings(as.numeric(support_values))
                display_values <- ifelse(is.na(numeric_values), support_values, round(numeric_values, 1))
                
                p <- p + ggtree::geom_nodelab(
                  aes(label = label),
                  subset = (!is.na(label) & label != "" & label != "NA"),
                  hjust = 1.2, 
                  vjust = -0.5, 
                  size = text_sz * 0.8, 
                  color = "red",
                  fontface = "bold"
                )
              } else {
                # For other support attributes, we need to manually add them
                # Create a data frame with node support
                n_nodes <- ape::Nnode(tree)
                n_tips <- ape::Ntip(tree)
                
                # Support values usually correspond to internal nodes
                if (length(support_values) == n_nodes) {
                  # Create node support data
                  node_support <- data.frame(
                    node = (n_tips + 1):(n_tips + n_nodes),
                    support = round(as.numeric(support_values), 1)
                  )
                  
                  # Filter out low/invalid support values
                  node_support <- node_support[!is.na(node_support$support) & node_support$support > 0, ]
                  
                  if (nrow(node_support) > 0) {
                    p <- p + ggtree::geom_nodelab(
                      data = node_support,
                      aes(x = x, y = y, label = support),
                      hjust = 1.2,
                      vjust = -0.5,
                      size = text_sz * 0.6,
                      color = "red",
                      fontface = "bold"
                    )
                  }
                }
              }
              
              cat("Node support labels added successfully\n")
              
            }, error = function(e) {
              cat("Failed to add node support labels:", e$message, "\n")
              
              # Fallback: add simple text annotations
              tryCatch({
                if (support_source == "node.label" && length(support_values) > 0) {
                  # Simple approach: just show that support exists
                  p <<- p + ggtree::geom_text2(
                    aes(subset = !isTip, label = "•"),
                    hjust = 0.5,
                    vjust = 0.5,
                    size = text_sz * 0.8,
                    color = "red"
                  )
                  cat("Added simple node markers as fallback\n")
                }
              }, error = function(e2) {
                cat("Even fallback node support failed:", e2$message, "\n")
              })
            })
          } else {
            cat("No node support values found in tree\n")
            
            # Add a visual indicator that no support was found
            if (requireNamespace("ggtree", quietly = TRUE)) {
              tryCatch({
                p <- p + ggtree::geom_nodepoint(
                  color = "orange", 
                  size = 1, 
                  alpha = 0.6
                )
                cat("Added node points to indicate internal nodes\n")
              }, error = function(e) {
                cat("Could not add node points:", e$message, "\n")
              })
            }
          }
        }
        
        print(p)
        
      } else {
        # IMPROVED FALLBACK TO APE PLOTTING
        par(mar = c(1, 1, 2, 1), bg = "white")
        
        tree <- clean_phylo_tree(tree)
        
        if (layout == "circular") {
          # Improved circular plot with ape
          n_tips <- ape::Ntip(tree)
          
          # Calculate optimal label size
          label_cex <- case_when(
            n_tips <= 30 ~ text_sz * 0.8,
            n_tips <= 60 ~ text_sz * 0.6,
            n_tips <= 100 ~ text_sz * 0.4,
            TRUE ~ text_sz * 0.3
          )
          
          ape::plot.phylo(
            tree, 
            type = "fan", 
            cex = label_cex,
            main = paste("Circular Phylogenetic Tree (", n_tips, " tips)"),
            show.tip.label = show_tips && n_tips <= 100,
            show.node.label = show_nodes,
            edge.color = "darkblue",
            edge.width = line_sz * 2,
            tip.color = "darkred"
          )
          
          if (show_branches && !is.null(tree$edge.length)) {
            ape::edgelabels(round(tree$edge.length, 3), cex = 0.5, col = "gray50")
          }
          
        } else if (layout == "unrooted") {
          ape::plot.phylo(
            tree,
            type = "unrooted",
            cex = text_sz * 0.8,
            show.tip.label = show_tips && ape::Ntip(tree) <= 50,
            show.node.label = show_nodes,
            main = paste("Unrooted Phylogenetic Tree (", ape::Ntip(tree), " tips)"),
            edge.color = "darkblue",
            edge.width = line_sz * 2
          )
          
          if (show_branches && !is.null(tree$edge.length)) {
            ape::edgelabels(round(tree$edge.length, 3), cex = 0.6, col = "gray50")
          }
          
        } else {
          ape::plot.phylo(
            tree, 
            cex = text_sz * 0.8,
            main = paste("Phylogenetic Tree (", ape::Ntip(tree), " tips)"),
            show.tip.label = show_tips,
            show.node.label = show_nodes,
            edge.color = "darkblue",
            edge.width = line_sz * 2
          )
          
          if (show_branches && !is.null(tree$edge.length)) {
            ape::edgelabels(round(tree$edge.length, 3), cex = 0.6, col = "gray50")
          }
        }
      }
      
    }, error = function(e) {
      cat("Plot error:", e$message, "\n")
      
      # Enhanced error handling
      par(mar = c(5, 5, 5, 5), bg = "white")
      plot(1, 1, type = "n",
           main = "Tree Visualization Error",
           xlab = "", ylab = "",
           xlim = c(0, 2), ylim = c(0, 2))
      
      text(1, 1.5, paste("Error:", e$message), col = "red", cex = 1)
      
      if (grepl("non-numeric|binary operator", e$message)) {
        text(1, 1.2, "Suggestion: Class conflict detected - try restarting R session",
             col = "blue", cex = 0.9)
        text(1, 1.0, "Alternative: Use rectangular layout for more stable visualization",
             col = "blue", cex = 0.9)
      } else if (grepl("class", e$message)) {
        text(1, 1.2, "Suggestion: Package conflict detected",
             col = "blue", cex = 0.9)
        text(1, 1.0, "Try: Restart R session or use rectangular layout",
             col = "blue", cex = 0.9)
      } else {
        text(1, 1.2, "Suggestion: Try rectangular layout for more stable visualization",
             col = "blue", cex = 0.9)
      }
      
      if (!is.null(tree)) {
        text(1, 0.7, paste("Tree has", ape::Ntip(tree), "tips"),
             col = "gray", cex = 0.8)
      }
    })
    
  }, height = 940)
  
  # Information output
  output$phylo_tree_info <- renderUI({
    if(phylo_tree_data$loaded) {
      
      original_tree <- phylo_tree_data$tree
      current_tree <- processed_tree()
      
      tagList(
        h4("Tree Information"),
        p(strong("Original Tree:")),
        p(paste("• Tips:", ape::Ntip(original_tree))),
        p(paste("• Nodes:", ape::Nnode(original_tree))),
        
        if(!is.null(input$tax_level) && input$tax_level != "none") {
          tagList(
            hr(),
            p(strong("Processed Tree:")),
            p(paste("• Taxonomic level:", input$tax_level)),
            p(paste("• Current tips:", ape::Ntip(current_tree))),
            p(paste("• Reduction:", round((1 - ape::Ntip(current_tree)/ape::Ntip(original_tree)) * 100, 1), "%"))
          )
        } else {
          p("• No taxonomic grouping applied")
        },
        
        hr(),
        p(strong("Debug Info:")),
        p(phylo_tree_data$debug_info,
          style = "font-family: monospace; font-size: 0.9em; color: #666;")
      )
      
    } else {
      p("No tree loaded")
    }
  })
  
  # Status output
  output$tree_status_message <- renderUI({
    if(phylo_tree_data$loaded) {
      div(
        h4("✓ Tree loaded successfully", style = "color: #28a745;"),
        p(paste("File:", basename(PHYLO_TREE_PATH)), style = "color: #666;")
      )
    } else if(phylo_tree_data$file_exists) {
      div(
        h4("⚠ Error loading tree", style = "color: #dc3545;"),
        p(phylo_tree_data$error_msg, style = "color: #666;")
      )
    } else {
      div(
        h4("✗ Tree file not found", style = "color: #dc3545;"),
        p("Check file path and permissions", style = "color: #666;")
      )
    }
  })
  
  # Reactive for status
  output$tree_loaded <- reactive({
    phylo_tree_data$loaded
  })
  outputOptions(output, "tree_loaded", suspendWhenHidden = FALSE)
  
  # Load tree at startup
  observeEvent(TRUE, {
    load_phylo_tree()
  }, once = TRUE)
  
  # Reload button with conflict resolution
  observeEvent(input$reload_phylo, {
    resolve_phylo_conflicts()  # Resolve conflicts before reloading
    load_phylo_tree()
  })
  
  # DEBUG: Show changes in taxonomic level
  observeEvent(input$tax_level, {
    if(phylo_tree_data$loaded) {
      cat("\n*** TAX LEVEL CHANGED TO:", input$tax_level, "***\n")
    }
  })
  
  # DEBUG: Show changes in minimum size
  observeEvent(input$min_group_size, {
    if(phylo_tree_data$loaded) {
      cat("\n*** MIN GROUP SIZE CHANGED TO:", input$min_group_size, "***\n")
    }
  })
  
  ########################################################################
  ########################################################################


  output$generate_report <- downloadHandler(
    filename = function() {
      paste0("Microbiome_Report_", Sys.Date(), ".", input$report_format)
    },
    content = function(file) {
      
      alpha_plot_path <- file.path(tempdir(), "alpha_diversity_plot.png")
      beta_plot_path <- file.path(tempdir(), "beta_diversity_plot.png")
      parallel_plot_path <- file.path(tempdir(), "parallel_plot.png")
      
      try({
        req(data_store$sample_metadata, data_store$taxonomy_data, input$parallel_vars)
        
        taxa_data <- data_store$taxonomy_data %>%
          group_by(Sample_ID) %>%
          summarize(Total_Abundance = sum(Abundance), .groups = "drop")
        
        combined <- merge(data_store$sample_metadata, taxa_data, by = "Sample_ID")
        vars <- input$parallel_vars
        all_vars <- c(vars, "Total_Abundance")
        
        plot_data <- combined %>%
          select(all_of(all_vars)) %>%
          na.omit() %>%
          mutate(across(where(is.character), as.factor)) %>%
          mutate(across(where(is.factor), ~ as.numeric(as.factor(.))))
        
        plot_obj <- plot_ly(
          type = 'parcoords',
          line = list(
            color = plot_data$Total_Abundance,
            colorscale = 'Viridis',
            showscale = TRUE
          ),
          dimensions = lapply(names(plot_data), function(col) {
            if (is.character(combined[[col]]) || is.factor(combined[[col]])) {
              f <- factor(combined[[col]])
              list(
                label = col,
                values = as.numeric(f),
                tickvals = seq_along(levels(f)),
                ticktext = levels(f)
              )
            } else {
              list(
                label = col,
                values = plot_data[[col]]
              )
            }
          })
        )
        # Save as temporary HTML and snapshot as PNG
        html_temp <- tempfile(fileext = ".html")
        saveWidget(as_widget(plot_obj), file = html_temp, selfcontained = TRUE)
        webshot2::webshot(url = html_temp, file = parallel_plot_path, vwidth = 1200, vheight = 800)
        message("✅ Parallel plot saved.")
      })
      
      # === REGENERATE p_alpha ===
      try({
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
        
        diversity_df <- bind_rows(
          sample_div %>%
            left_join(sample_meta, by = "Sample_ID") %>%
            mutate(GroupLabel = if (input$subdivide_samples) .[[input$condition_column]] else "Sample"),
          control_div
        )
        
        metric_col <- switch(input$alpha_metric,
                             "Observed OTUs" = "Observed",
                             "Shannon" = "Shannon",
                             "Simpson" = "Simpson")
        
        p_alpha <- ggplot(diversity_df, aes(x = GroupLabel, y = .data[[metric_col]], fill = GroupLabel)) +
          geom_boxplot(alpha = 0.6, outlier.shape = NA) +
          geom_jitter(shape = 21, alpha = 0.6, color = "black", width = 0.2) +
          labs(title = paste(input$alpha_metric, "Diversity across Groups"),
               x = "Group", y = paste(input$alpha_metric, "Index")) +
          theme_minimal() +
          theme(
            legend.position = "right",
            plot.title = element_text(size = 14, hjust = 0.5),
            axis.title = element_text(size = 12),
            axis.text = element_text(size = 10),
            legend.text = element_text(size = 10),
            legend.title = element_text(size = 11)
          ) +
          scale_fill_brewer(palette = input$alpha_color_palette)
        
        ggsave(alpha_plot_path, plot = p_alpha, width = 10, height = 6, dpi = 300, bg = "white")
        message("✅ Alpha plot saved.")
      })
      
      # === REGENERATE p_beta ===
      try({
        all_tax <- bind_rows(
          data_store$taxonomy_data,
          data_store$control_taxonomy_data
        ) %>%
          select(Sample_ID, Species, Abundance) %>%
          pivot_wider(names_from = Species, values_from = Abundance, values_fill = 0) %>%
          column_to_rownames("Sample_ID")
        
        dist_method <- switch(input$distance_metric,
                              "Bray-Curtis" = "bray",
                              "Euclidean" = "euclidean",
                              "Jaccard" = "jaccard",
                              "Canberra" = "canberra",
                              "Manhattan" = "manhattan",
                              "Kulczynski" = "kulczynski",
                              "Chord" = "chord")
        
        ord <- ape::pcoa(vegan::vegdist(all_tax, method = dist_method))
        coords <- ord$vectors[, 1:2]
        ord_df <- as.data.frame(coords)
        colnames(ord_df) <- c("Axis1", "Axis2")
        ord_df$Sample_ID <- rownames(coords)
        
        all_meta <- bind_rows(data_store$sample_metadata, data_store$control_sample_metadata)
        ord_df <- left_join(ord_df, all_meta, by = "Sample_ID")
        
        ord_df$GroupLabel <- if (input$subdivide_samples) {
          ifelse(ord_df$Sample_ID %in% data_store$control_sample_metadata$Sample_ID,
                 "Control", ord_df[[input$condition_column]])
        } else {
          ifelse(ord_df$Sample_ID %in% data_store$control_sample_metadata$Sample_ID, "Control", "Sample")
        }
        
        p_beta <- ggplot(ord_df, aes(x = Axis1, y = Axis2, color = GroupLabel)) +
          geom_point(size = 3, alpha = 0.8)
        
        if (isTRUE(input$show_group_ellipses)) {
          p_beta <- p_beta +
            stat_ellipse(aes(group = GroupLabel),
                         type = "norm", linetype = "dashed",
                         alpha = 0.6, size = 1, show.legend = FALSE)
        }
        
        p_beta <- p_beta +
          labs(title = paste0("PCoA on ", input$distance_metric, " Distance"),
               x = "Axis 1", y = "Axis 2") +
          theme_minimal() +
          theme(
            legend.position = "right",
            plot.title = element_text(size = 14, hjust = 0.5),
            axis.title = element_text(size = 12),
            axis.text = element_text(size = 10),
            legend.text = element_text(size = 10),
            legend.title = element_text(size = 11)
          )
        
        ggsave(beta_plot_path, plot = p_beta, width = 10, height = 6, dpi = 300, bg = "white")
        message("✅ Beta plot saved.")
      })
      
      # === RENDER REPORT ===
      rmarkdown::render(
        input = "report_template.Rmd",
        output_file = file,
        params = list(
          report_title = input$report_title,
          report_comments = input$report_comments,
          selected_sections = input$report_sections,
          alpha_plot = alpha_plot_path,
          beta_plot = beta_plot_path,
          parallel_plot = parallel_plot_path
        ),
        envir = new.env(parent = globalenv())
      )
    }
  )
  
  
  
  
  
  
  
  
  
}

  

# Run app
shinyApp(ui, server)





