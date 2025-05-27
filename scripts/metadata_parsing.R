# Load necessary libraries
library(readr)
library(dplyr)
library(stringr)

##########################################
# Install libraries if not available
##########################################
if (!require("tidyverse"))    install.packages("tidyverse")
if (!require("lubridate"))    install.packages("lubridate")
if (!require("janitor"))      install.packages("janitor")

library(tidyverse)
library(lubridate)
library(janitor)

# ─────────────────────────────────────────────────────────────────────────────
# 1) Read the sample and run CSV files ---------------------------------------
# ─────────────────────────────────────────────────────────────────────────────
print("Reading metadata files...")

# Read sample metadata
metadata_sample <- read_csv("../metadata/metadata_sample.csv", 
                            locale = locale(encoding = "UTF-8"),
                            show_col_types = FALSE)

# Read run metadata  
metadata_run <- read_csv("../metadata/metadata_run.csv",
                         locale = locale(encoding = "UTF-8"),
                         show_col_types = FALSE)

print(paste("Sample metadata rows:", nrow(metadata_sample)))
print(paste("Run metadata rows:", nrow(metadata_run)))

# ─────────────────────────────────────────────────────────────────────────────
# 2) Merge the two datasets based on RunID -----------------------------------
# ─────────────────────────────────────────────────────────────────────────────
print("Merging datasets by RunID...")

# Check if RunID exists in both datasets
if (!"RunID" %in% colnames(metadata_sample)) {
  stop("RunID column not found in metadata_sample.csv")
}
if (!"RunID" %in% colnames(metadata_run)) {
  stop("RunID column not found in metadata_run.csv")
}

# Perform left join to keep all samples and add run information
metadata <- metadata_sample %>%
  left_join(metadata_run, by = "RunID", suffix = c("", "_run"))

print(paste("Merged metadata rows:", nrow(metadata)))

# ─────────────────────────────────────────────────────────────────────────────
# 3) Clean up column names ---------------------------------------------------
# ─────────────────────────────────────────────────────────────────────────────
clean_names_safe <- function(names_vec) {
  names_vec %>%
    str_trim() %>%
    str_replace_all("[^A-Za-z0-9_()ºC]+", "_") %>%
    str_replace_all("_+", "_") %>%
    str_replace_all("^_|_$", "")
}

colnames(metadata) <- clean_names_safe(colnames(metadata))

# ─────────────────────────────────────────────────────────────────────────────
# 4) Fill empty cells with "N/A" and clean text fields ----------------------
# ─────────────────────────────────────────────────────────────────────────────
print("Cleaning and formatting data...")

metadata <- metadata %>%
  mutate(across(
    .cols = where(is.character),
    .fns  = ~ if_else(is.na(.x) | str_trim(.x) == "", "N/A", .x)
  ))

# ─────────────────────────────────────────────────────────────────────────────
# 4.2) Trim and replace spaces in text fields --------------------------------
# ─────────────────────────────────────────────────────────────────────────────
metadata <- metadata %>%
  mutate(across(
    .cols = where(is.character),
    .fns  = ~ str_trim(.x) %>% str_replace_all("\\s+", "_")
  ))

# ─────────────────────────────────────────────────────────────────────────────
# 4.3) Apply Title Case to selected columns ----------------------------------
# ─────────────────────────────────────────────────────────────────────────────
title_cols <- c("Institution", "Department", "Analyst_Processor_Name", 
                "Notes_Samples", "Technician_name", "Notes_Runs")

metadata <- metadata %>%
  mutate(across(
    .cols = any_of(title_cols),
    .fns  = ~ if_else(.x == "N/A", .x, 
                      .x %>%
                        str_replace_all("_", " ") %>%
                        str_to_title() %>%
                        str_replace_all("\\s+", "_"))
  ))

# ─────────────────────────────────────────────────────────────────────────────
# 5) Field-specific formatting ------------------------------------------------
# ─────────────────────────────────────────────────────────────────────────────


clean_text <- function(x) {
  x <- str_trim(x)
  x <- ifelse(is.na(x) | x == "", "N/A", x)
  x <- str_replace_all(x, " +", "_")
  str_to_title(x)
}


metadata <- metadata %>%
  mutate(across(everything(), ~ ifelse(is.na(.) | . == "", "N/A", .))) %>%
  mutate(
    Sample_ID = str_trim(Sample_ID),
    
    Institution = clean_text(Institution),
    Department = clean_text(Department),
    Collection_Date = dmy(Collection_Date) %>% format("%d/%m/%Y"),
    
    Collection_Storage_Temperature = case_when(
      str_detect(Collection_Storage_Temperature, "-?80") ~ "-80",
      str_detect(Collection_Storage_Temperature, "-?20") ~ "-20",
      Collection_Storage_Temperature == "4" ~ "4",
      str_to_lower(Collection_Storage_Temperature) %in% c("room temperature", "room_temperature") ~ "Room_Temperature",
      TRUE ~ "N/A"
    ),
    
    Analyst_Processor_Name = clean_text(Analyst_Processor_Name),
    
    Gender = case_when(
      str_to_lower(Gender) == "female" ~ "Female",
      str_to_lower(Gender) == "male" ~ "Male",
      TRUE ~ "Not_specified"
    ),
    
    Age = suppressWarnings(as.integer(Age)),
    Age = if_else(is.na(Age), "N/A", as.character(Age)),
    
    Ongoing_conditions = case_when(
      str_detect(str_to_lower(Ongoing_conditions), "acid") ~ "Acid_reflux",
      str_detect(str_to_lower(Ongoing_conditions), "h[.]? pylori|gastritis") ~ "Helicobacter_pylori_associated_gastritis",
      str_detect(str_to_lower(Ongoing_conditions), "ibs") ~ "Irritable_bowel_syndrome(IBS)",
      str_detect(str_to_lower(Ongoing_conditions), "type[ _]?1") ~ "Diabetes(Type 1)",
      str_detect(str_to_lower(Ongoing_conditions), "diabetes") ~ "Diabetes(Mellitus)",
      str_detect(str_to_lower(Ongoing_conditions), "metabolic") ~ "Metabolic_syndromes",
      str_detect(str_to_lower(Ongoing_conditions), "nafld") ~ "Non_alcoholic_fatty_liver_disease(NAFLD)",
      str_detect(str_to_lower(Ongoing_conditions), "fld|fatty") ~ "Fatty_liver_disease(FLD)",
      str_detect(str_to_lower(Ongoing_conditions), "ibd") ~ "Inflammatory_Bowel_Disease(IBD)",
      str_detect(str_to_lower(Ongoing_conditions), "lupus|sle") ~ "Lupus(SLE)",
      str_detect(str_to_lower(Ongoing_conditions), "crohn") ~ "Crohn's_disease",
      str_detect(str_to_lower(Ongoing_conditions), "rheumatoid") ~ "Rheumatoid_arthritis",
      str_detect(str_to_lower(Ongoing_conditions), "ache") ~ "Stomach_ache",
      TRUE ~ Ongoing_conditions
    ),
    
    Appendix_removed = case_when(
      str_to_lower(Appendix_removed) == "yes" ~ "Yes",
      str_to_lower(Appendix_removed) == "no" ~ "No",
      TRUE ~ "N/A"
    ),
    
    Allergies = case_when(
      str_detect(Allergies, "peanut") ~ "Peanuts",
      str_detect(Allergies, "shellfish") ~ "Shellfish",
      str_detect(Allergies, "tree") ~ "Tree_nuts",
      str_detect(Allergies, "egg") ~ "Eggs",
      str_detect(Allergies, "milk") ~ "Milk",
      str_detect(Allergies, "free") ~ "Allergy free",
      str_detect(Allergies, "unspecified") ~ "Unspecified",
      TRUE ~ "N/A"
    ),
    
    Dietary_Information = case_when(
      str_detect(Dietary_Information, "vegetarian.*seafood") ~ "Vegetarian_but_eat_seafood",
      str_detect(Dietary_Information, "vegetarian") ~ "Vegetarian",
      str_detect(Dietary_Information, "vegan") ~ "Vegan",
      str_detect(Dietary_Information, "gluten") ~ "Gluten_free",
      str_detect(Dietary_Information, "keto") ~ "Keto",
      str_detect(Dietary_Information, "halal") ~ "Halal",
      str_detect(Dietary_Information, "omnivore") ~ "Omnivore",
      TRUE ~ "N/A"
    ),
    
    Bowel_movement_quality = case_when(
      str_detect(Bowel_movement_quality, "constipat") ~ "Constipated",
      str_detect(Bowel_movement_quality, "normal") ~ "Normal",
      str_detect(Bowel_movement_quality, "diarrh") ~ "Diarrhea",
      TRUE ~ "N/A"
    ),
    
    Antibiotic_intake = case_when(
      str_detect(Antibiotic_intake, "week") ~ "Week",
      str_detect(Antibiotic_intake, "month") ~ "Month",
      str_detect(Antibiotic_intake, "6") ~ "6 months",
      str_detect(Antibiotic_intake, "year") & str_detect(Antibiotic_intake, "past") ~ "Past_year",
      str_detect(Antibiotic_intake, "year") ~ "Year",
      TRUE ~ "N/A"
    ),
    
    Medications = case_when(
      str_detect(Medications, "diabet") ~ "Antidiabetics",
      str_detect(Medications, "probiot") ~ "Probiotics",
      str_detect(Medications, "prebiot") ~ "Prebiotics",
      str_detect(Medications, "laxative") ~ "Laxatives",
      str_detect(Medications, "mimetic") ~ "Antimimetic",
      str_detect(Medications, "ppi|pump") ~ "PPIs(proton_pump_inhibitors)",
      str_detect(Medications, "immuno") ~ "Immunosupressors",
      str_detect(Medications, "antidepress|antipsych|anxio") ~ "Antidepressors/Antipsicotics/Anxiolytics",
      str_detect(Medications, "contracep") ~ "Contraceptives",
      TRUE ~ "N/A"
    ),
    
    Cancer = case_when(
      str_to_lower(Cancer) == "yes" ~ "Yes",
      str_to_lower(Cancer) == "no" ~ "No",
      TRUE ~ "N/A"
    ),
    
    Body_Mass_Index = case_when(
      str_detect(Body_Mass_Index, "under") ~ "Underweight<18.5",
      str_detect(Body_Mass_Index, "normal") ~ "Normal_Weight_18.5-24.9",
      str_detect(Body_Mass_Index, "over") ~ "Overweight_25-29.9",
      str_detect(Body_Mass_Index, "obese") ~ "Obese_30<",
      TRUE ~ "N/A"
    ),
    
    Exercise_frequency = case_when(
      str_detect(Exercise_frequency, "daily") ~ "Daily",
      str_detect(Exercise_frequency, "1-2") ~ "1-2_times_per_week",
      str_detect(Exercise_frequency, "3-5") ~ "3-5 times_per_week",
      str_detect(Exercise_frequency, "rare") ~ "Rarely",
      str_detect(Exercise_frequency, "never") ~ "Never",
      TRUE ~ "N/A"
    ),
    
    Smoking_status = case_when(
      str_detect(Smoking_status, "non") ~ "Non-smoker",
      str_detect(Smoking_status, "smok") ~ "Smoker",
      TRUE ~ "N/A"
    ),
    
    Daily_cigarettes = case_when(
      str_detect(Daily_cigarettes, "1[-–]5") ~ "1-5",
      str_detect(Daily_cigarettes, "6[-–]10") ~ "6-10",
      str_detect(Daily_cigarettes, "11[-–]15") ~ "11-15",
      str_detect(Daily_cigarettes, "16[-–]20") ~ "16-20",
      str_detect(Daily_cigarettes, "20\\+|\\+20") ~ "+20",
      TRUE ~ "N/A"
    ),
    
    Alcohol_consumption = case_when(
      str_to_lower(Alcohol_consumption) == "yes" ~ "Yes",
      str_to_lower(Alcohol_consumption) == "no" ~ "No",
      TRUE ~ "N/A"
    ),
    
    Frequency_of_alcohol_consumption = case_when(
      str_detect(Frequency_of_alcohol_consumption, "(?i)daily") ~ "Daily",
      str_detect(Frequency_of_alcohol_consumption, "(?i)3[-–]5") ~ "3-5_per_week",
      str_detect(Frequency_of_alcohol_consumption, "(?i)1[-–]2") ~ "1-2_per_week",
      str_detect(Frequency_of_alcohol_consumption, "(?i)1[-–]3") ~ "1-3_per_month",
      TRUE ~ "N/A"
    ),
    
    Sequencing_Platform = case_when(
      str_detect(Sequencing_Platform, "(?i)miseq") ~ "MiSeq",
      str_detect(Sequencing_Platform, "(?i)nextseq") ~ "NextSeq",
      str_detect(Sequencing_Platform, "(?i)grindion") ~ "GrindION",
      str_detect(Sequencing_Platform, "(?i)minion") ~ "MinION",
      str_detect(Sequencing_Platform, "(?i)ion.*genestudio") ~ "Ion_GeneStudio_S5",
      TRUE ~ "N/A"
    ),
    
    Sequencing_Type = case_when(
      str_detect(Sequencing_Type, "(?i)16s") ~ "16S_rRNA",
      str_detect(Sequencing_Type, "(?i)wgs") ~ "WGS",
      str_detect(Sequencing_Type, "(?i)shotgun") ~ "Shotgun_metagenomic",
      str_detect(Sequencing_Type, "(?i)rna") ~ "RNA-Seq",
      str_detect(Sequencing_Type, "(?i)long") ~ "Long-read",
      str_detect(Sequencing_Type, "(?i)metatrans") ~ "Metatranscriptomics",
      TRUE ~ "N/A"
    ),
    
    Sequencing_depth_target = case_when(
      str_detect(Sequencing_depth_target, "(?i)<1") ~ "<1_million_reads/sample",
      str_detect(Sequencing_depth_target, "(?i)1.*5") ~ "1-5_million_reads/sample",
      str_detect(Sequencing_depth_target, "(?i)5.*10") ~ "5-10_million_reads/sample",
      str_detect(Sequencing_depth_target, "(?i)10.*20") ~ "10-20_million_reads/sample",
      str_detect(Sequencing_depth_target, "(?i)20.*50") ~ "20-50_million_reads/sample",
      str_detect(Sequencing_depth_target, "(?i)50") ~ "50_million_reads/sample",
      TRUE ~ "N/A"
    ),
    
    Library_preparation_kit = case_when(
      str_detect(Library_preparation_kit, "(?i)nextera.*xt") ~ "Nextera_XT",
      str_detect(Library_preparation_kit, "(?i)nextera.*flex") ~ "Nextera_DNA_Flex",
      str_detect(Library_preparation_kit, "(?i)truseq") ~ "TruSeq_Nano",
      str_detect(Library_preparation_kit, "(?i)nebnext") ~ "NEBNext_Ultra_II",
      str_detect(Library_preparation_kit, "(?i)swift") ~ "Swift_Biosciences_16S",
      str_detect(Library_preparation_kit, "(?i)qia") ~ "QIAseq",
      TRUE ~ "N/A"
    ),
    
    Technician_name = str_to_title(str_replace_all(Technician_name, "_", " ")) %>%
      str_replace_all("\\s+", "_"),
    
    Notes_Runs = str_to_title(str_replace_all(Notes_Runs, "_", " ")) %>%
      str_replace_all("\\s+", "_")
  )
                 
# ─────────────────────────────────────────────────────────────────────────────
# 6) Write out the cleaned CSV ------------------------------------------------
# ─────────────────────────────────────────────────────────────────────────────
print("Writing cleaned metadata...")
write_csv(metadata, "../metadata/template_metadata_merged/metadata_cleaned.csv", na = "")

print("✓ Successfully created metadata_cleaned.csv")
print(paste("Final dataset contains", nrow(metadata), "rows and", ncol(metadata), "columns"))

# Display summary
print("\n=== SUMMARY ===")
print(paste("Sample metadata entries:", nrow(metadata_sample)))
print(paste("Run metadata entries:", nrow(metadata_run))) 
print(paste("Merged entries:", nrow(metadata)))
print("\nColumn names in final dataset:")
print(colnames(metadata))