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
metadata_sample_raw <- read_csv("metadata/healthy_controls/non_curated/metadata_sample.csv", 
                            locale = locale(encoding = "UTF-8"),
                            show_col_types = FALSE)

# Read run metadata  
metadata_run_raw <- read_csv("metadata/healthy_controls/non_curated/metadata_run.csv",
                         locale = locale(encoding = "UTF-8"),
                         show_col_types = FALSE)

print(paste("Sample metadata rows:", nrow(metadata_sample_raw)))
print(paste("Run metadata rows:", nrow(metadata_run_raw)))

# ─────────────────────────────────────────────────────────────────────────────
# 1.5) Remove Duplicated entries ----------------------------------------------
# ─────────────────────────────────────────────────────────────────────────────

print("Deleting unconsistent records...")


metadata_sample <- metadata_sample_raw %>%
  group_by(Sample_Accession) %>%
  filter(!(n() > 1 & (n_distinct(Age) > 1 | n_distinct(Gender) > 1))) %>%
  ungroup()

metadata_run <- metadata_run_raw %>%
  filter(Sample_Accession %in% metadata_sample$Sample_Accession)

# ─────────────────────────────────────────────────────────────────────────────
# 2) Merge the two datasets based on SampleID -----------------------------------
# ─────────────────────────────────────────────────────────────────────────────
print("Merging datasets by Sample_Accession...")

# Check if SampleID exists in both datasets
if (!"Sample_Accession" %in% colnames(metadata_sample)) {
  stop("Sample_Accession column not found in metadata_sample.csv")
}
if (!"Sample_Accession" %in% colnames(metadata_run)) {
  stop("Sample_Accession column not found in metadata_run.csv")
}

# Perform left join to keep all samples and add run information
metadata <- metadata_sample %>%
  left_join(metadata_run, by = "Sample_Accession", suffix = c("", "_run"))

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
names(metadata)[names(metadata) == "Ongoing_conditions_"] <- "Ongoing_conditions"
names(metadata)[names(metadata) == "Sample_ID"] <- "Sample ID"
# ─────────────────────────────────────────────────────────────────────────────
# 4) Helper functions for consistent formatting -------------------------------
# ─────────────────────────────────────────────────────────────────────────────

# Function to clean text with proper title case and underscores
clean_text <- function(x) {
  x <- str_trim(as.character(x))
  x <- ifelse(is.na(x) | x == "" | x == "NA", "N/A", x)
  # Don't process if already N/A
  x <- ifelse(x == "N/A", x, 
              x %>%
                str_replace_all("\\s+", "_") %>%
                str_to_title())
  return(x)
}

# Function to handle numeric BMI conversion to categories
convert_bmi <- function(bmi_value) {
  # Try to convert to numeric first
  bmi_num <- suppressWarnings(as.numeric(as.character(bmi_value)))
  
  if (is.na(bmi_num)) {
    return("N/A")
  } else if (bmi_num < 18.5) {
    return("Underweight<18.5")
  } else if (bmi_num >= 18.5 & bmi_num <= 24.9) {
    return("Normal_Weight_18.5-24.9")
  } else if (bmi_num >= 25 & bmi_num <= 29.9) {
    return("Overweight_25-29.9")
  } else if (bmi_num >= 30) {
    return("Obese_>30")
  } else {
    return("N/A")
  }
}

convert_alcohol_frequency <- function(freq_vector) {
  sapply(freq_vector, function(freq) {
    if (is.na(freq) || freq == "N/A" || freq == "") {
      return("N/A")
    }
    
    # Remove quotes and whitespace
    freq <- trimws(gsub('"', '', as.character(freq)))
    
    # Handle range patterns (e.g., "3-5", "7-10")
    if (grepl('-', freq)) {
      # Extract numbers from range
      numbers <- as.numeric(unlist(regmatches(freq, gregexpr('\\d+', freq))))
      
      if (length(numbers) == 2) {
        start <- numbers[1]
        end <- numbers[2] 
        mean_val <- mean(c(start, end))
        
        # If mean > 7, classify as daily
        if (mean_val >= 7) {
          return("Daily")
        } else {
          # Round to nearest integer to avoid decimals
          return(paste0(round(mean_val), "_per_week"))
        }
      }
    }
    # Handle "< X" patterns
    else if (grepl('<', freq)) {
      # Extract the number after <
      number <- as.numeric(gsub('[^0-9.]', '', freq))
      if (!is.na(number)) {
        # Use 1 instead of 0.5 for "< X" patterns, or round up
        val <- max(1, ceiling(number / 2))
        return(paste0(val, "_per_week"))
      }
    }
    # Handle single numbers
    else {
      number <- as.numeric(freq)
      if (!is.na(number)) {
        if (number >= 7) {
          return("Daily")
        } else {
          return(paste0(number, "_per_week"))
        }
      }
    }
    
    return("N/A")
  }, USE.NAMES = FALSE)
}
# ─────────────────────────────────────────────────────────────────────────────
# 5) Fill empty cells with "N/A" and apply field-specific formatting ---------
# ─────────────────────────────────────────────────────────────────────────────
print("Cleaning and formatting data...")

metadata <- metadata %>%
  # First, replace all empty/NA values with "N/A"
  mutate(across(-Collection_Date, ~ ifelse(is.na(.) | . == "" | . == "NA", "N/A", as.character(.)))) %>%
  mutate(
    # Sample_Accession - keep as is (should already be in correct format)
    Sample_Accession = str_trim(Sample_Accession),
    
    # Institution - First letter uppercase with underscores
    Institution = clean_text(Institution),
    
    # Department - First letter uppercase with underscores  
    Department = clean_text(Department),
    
    # Collection_Date can be processed but is giving errors so we omit
    
    # Collection_Storage_Temperature - Specific parameters
    Collection_Storage_Temperature = case_when(
      Collection_Storage_Temperature == "N/A" ~ "N/A",
      str_detect(Collection_Storage_Temperature, "-?80") ~ "-80",
      str_detect(Collection_Storage_Temperature, "-?20") ~ "-20",
      str_detect(Collection_Storage_Temperature, "^4$") ~ "4",
      str_detect(str_to_lower(Collection_Storage_Temperature), "room") ~ "Room_Temperature",
      TRUE ~ "N/A"
    ),
    
    # Analyst_Processor_Name - First letter uppercase with underscores
    Analyst_Processor_Name = clean_text(Analyst_Processor_Name),
    
    # Gender - Specific parameters
    Gender = case_when(
      Gender == "N/A" ~ "Not_specified",
      str_to_lower(Gender) %in% c("f", "female") ~ "Female",
      str_to_lower(Gender) %in% c("m", "male") ~ "Male",
      TRUE ~ "Not_specified"
    ),
    
    # Age - Numbers only
    Age = suppressWarnings(as.integer(Age)),
    Age = ifelse(is.na(Age), "N/A", as.character(Age)),
    
    # Ongoing_conditions - Specific parameters
    Ongoing_conditions = case_when(
      Ongoing_conditions == "N/A" ~ "N/A",
      str_detect(str_to_lower(Ongoing_conditions), "acid") ~ "Acid_reflux",
      str_detect(str_to_lower(Ongoing_conditions), "h[.]?\\s*pylori|gastritis") ~ "Helicobacter_pylori–associated_gastritis",
      str_detect(str_to_lower(Ongoing_conditions), "ibs") ~ "Irritable_bowel_syndrome_(IBS)",
      str_detect(str_to_lower(Ongoing_conditions), "type\\s*1") ~ "Diabetes(Type_1)",
      str_detect(str_to_lower(Ongoing_conditions), "diabetes") ~ "Diabetes_(Mellitus)",
      str_detect(str_to_lower(Ongoing_conditions), "metabolic") ~ "Metabolic_syndromes",
      str_detect(str_to_lower(Ongoing_conditions), "nafld") ~ "Non-alcoholic_fatty_liver_disease_(NAFLD)",
      str_detect(str_to_lower(Ongoing_conditions), "fld|fatty") ~ "Fatty_liver_disease_(FLD)",
      str_detect(str_to_lower(Ongoing_conditions), "ibd") ~ "Inflammatory_Bowel_Disease_(IBD)",
      str_detect(str_to_lower(Ongoing_conditions), "lupus|sle") ~ "Lupus_(SLE)",
      str_detect(str_to_lower(Ongoing_conditions), "crohn") ~ "Crohn's_disease",
      str_detect(str_to_lower(Ongoing_conditions), "rheumatoid") ~ "Rheumatoid_arthritis",
      str_detect(str_to_lower(Ongoing_conditions), "ache") ~ "Stomach_ache",
      TRUE ~ "N/A"
    ),
    
    # Appendix_removed - Yes/No parameters
    Appendix_removed = case_when(
      Appendix_removed == "N/A" ~ "N/A",
      str_to_lower(Appendix_removed) == "yes" ~ "Yes",
      str_to_lower(Appendix_removed) == "no" ~ "No",
      TRUE ~ "N/A"
    ),
    
    # UPDATED: Frequency_of_alcohol_consumption - Process BEFORE using in Alcohol_consumption
    Frequency_of_alcohol_consumption = case_when(
      Frequency_of_alcohol_consumption == "N/A" ~ "N/A",
      # First check for text-based entries
      str_detect(str_to_lower(Frequency_of_alcohol_consumption), "daily") ~ "Daily",
      str_detect(str_to_lower(Frequency_of_alcohol_consumption), "3[-–]5|3_5.*week") ~ "4_per_week",
      str_detect(str_to_lower(Frequency_of_alcohol_consumption), "1[-–]2|1_2.*week") ~ "2_per_week", 
      str_detect(str_to_lower(Frequency_of_alcohol_consumption), "1[-–]3|1_3.*month") ~ "1_per_week",
      # Then handle numeric values using the conversion function
      TRUE ~ convert_alcohol_frequency(Frequency_of_alcohol_consumption)
    ),
    
    # UPDATED: Alcohol_consumption - Now uses processed frequency data
    Alcohol_consumption = case_when(
      Frequency_of_alcohol_consumption == "N/A" ~ "N/A",
      Frequency_of_alcohol_consumption %in% c("Daily", "0_per_week") ~ case_when(
        Frequency_of_alcohol_consumption == "Daily" ~ "Yes",
        Frequency_of_alcohol_consumption == "0_per_week" ~ "No",
        TRUE ~ "N/A"
      ),
      str_detect(Frequency_of_alcohol_consumption, "_per_week$") ~ "Yes",
      TRUE ~ "N/A"
    ),
    
    # UPDATED: Alcohol_consumption - Now uses processed frequency data
    Alcohol_consumption = case_when(
      Frequency_of_alcohol_consumption == "N/A" ~ "N/A",
      Frequency_of_alcohol_consumption %in% c("Daily", "0_per_week") ~ case_when(
        Frequency_of_alcohol_consumption == "Daily" ~ "Yes",
        Frequency_of_alcohol_consumption == "0_per_week" ~ "No",
        TRUE ~ "N/A"
      ),
      str_detect(Frequency_of_alcohol_consumption, "_per_week$") ~ "Yes",
      TRUE ~ "N/A"
    ),
    
    # Allergies - Extract multiple allergy categories
    Allergies = str_to_lower(Allergies),
    Allergies = case_when(
      is.na(Allergies) | Allergies %in% c("na", "n/a", "none", "allergy free", "") ~ "Allergy_free",
      TRUE ~ map_chr(Allergies, function(x) {
        # Define allergy patterns and their corresponding labels
        allergy_patterns <- c(
          "Peanuts" = "peanut",
          "Shellfish" = "shellfish", 
          "Tree_nuts" = "tree.*nut|nut.*tree",
          "Eggs" = "\\begg\\b",
          "Milk" = "\\bmilk\\b|dairy",
          "Cats" = "\\bcat\\b|feline|cat dander",
          "Dogs" = "\\bdog\\b|canine",
          "Penicillin" = "penicillin|penecillian|pencillin|penicillian",
          "Amoxicillin" = "amoxicillin",
          "Zithromax" = "zithromax",
          "Ceclor" = "ceclor",
          "Accutane" = "accutane",
          "Sulfa_drugs" = "sulf|sulfa",
          "Latex" = "latex",
          "Benzoyl_peroxide" = "benzoyl\\s*peroxide",
          "Dust" = "dust",
          "Seasonal_allergies" = "pollen|seasonal|summer.*allerg|allerg.*summer",
          "Mold" = "\\bmold\\b|mould",
          "Asthma" = "asthma"
        )
        
        # Find matching allergies
        detected <- names(allergy_patterns)[map_lgl(allergy_patterns, ~ str_detect(x, .x))]
        
        # Return result
        if (length(detected) > 0) {
          paste(detected, collapse = ", ")
        } else {
          "Unspecified"
        }
      })
    ),
    
    # Dietary_Information - Specific parameters
    Dietary_Information = case_when(
      Dietary_Information == "N/A" ~ "N/A",
      str_detect(str_to_lower(Dietary_Information), "vegetarian.*seafood") ~ "Vegetarian_but_eat_seafood",
      str_detect(str_to_lower(Dietary_Information), "vegetarian") ~ "Vegetarian",
      str_detect(str_to_lower(Dietary_Information), "vegan") ~ "Vegan",
      str_detect(str_to_lower(Dietary_Information), "gluten") ~ "Gluten_free",
      str_detect(str_to_lower(Dietary_Information), "keto") ~ "Keto",
      str_detect(str_to_lower(Dietary_Information), "halal") ~ "Halal",
      str_detect(str_to_lower(Dietary_Information), "kosher") ~ "Kosher",
      str_detect(str_to_lower(Dietary_Information), "paleo") ~ "Paleo",
      str_detect(str_to_lower(Dietary_Information), "omnivore|avoid.*") ~ "Omnivore",
      TRUE ~ "N/A"
    ),
    
    # Bowel_movement_quality - Specific parameters
    Bowel_movement_quality = case_when(
      Bowel_movement_quality == "N/A" ~ "N/A",
      str_detect(str_to_lower(Bowel_movement_quality), "constipat") ~ "Constipated",
      str_detect(str_to_lower(Bowel_movement_quality), "normal") ~ "Normal",
      str_detect(str_to_lower(Bowel_movement_quality), "diarrh") ~ "Diarrhea",
      TRUE ~ "N/A"
    ),
    
    # Antibiotic_intake - Specific parameters
    Antibiotic_intake = case_when(
      Antibiotic_intake == "N/A" ~ "N/A",
      str_detect(str_to_lower(Antibiotic_intake), "week") ~ "Week",
      str_detect(str_to_lower(Antibiotic_intake), "month") & str_detect(str_to_lower(Antibiotic_intake), "6") ~ "6_months",
      str_detect(str_to_lower(Antibiotic_intake), "month") ~ "Month",
      str_detect(str_to_lower(Antibiotic_intake), "past.*year") ~ "Past_year",
      str_detect(str_to_lower(Antibiotic_intake), "year") ~ "Year",
      TRUE ~ "N/A"
    ),
    
    # Medications - Specific parameters
    Medications = case_when(
      Medications == "N/A" ~ "N/A",
      str_detect(str_to_lower(Medications), "antidiabet|diabet") ~ "Antidiabetics",
      str_detect(str_to_lower(Medications), "probiot") ~ "Probiotics",
      str_detect(str_to_lower(Medications), "prebiot") ~ "Prebiotics",
      str_detect(str_to_lower(Medications), "laxative") ~ "Laxatives",
      str_detect(str_to_lower(Medications), "antiemetic|mimetic") ~ "Antiemetic",
      str_detect(str_to_lower(Medications), "ppi|pump") ~ "PPIs_(proton-pump_inhibitors)",
      str_detect(str_to_lower(Medications), "immuno") ~ "Immunosupressors",
      str_detect(str_to_lower(Medications), "antidepress|antipsych|anxio") ~ "Antidepressors/Antipsicotics/Anxiolytics",
      str_detect(str_to_lower(Medications), "contracep") ~ "Contraceptives",
      str_detect(str_to_lower(Medications), "retinoid") ~ "Retinoids",
      str_detect(str_to_lower(Medications), "antihist") ~ "Antihistamines",
      str_detect(str_to_lower(Medications), "nsaid") ~ "NSAIDs",
      TRUE ~ "N/A"
    ),
    
    # Cancer - Yes/No parameters
    Cancer = case_when(
      Cancer == "N/A" ~ "N/A",
      str_to_lower(Cancer) == "yes" ~ "Yes",
      str_to_lower(Cancer) == "no" ~ "No",
      TRUE ~ "N/A"
    ),
    
    # Body_Mass_Index - Convert numeric to categories
    Body_Mass_Index = sapply(Body_Mass_Index, convert_bmi),
    
    # Exercise_frequency - Specific parameters
    Exercise_frequency = case_when(
      Exercise_frequency == "N/A" ~ "N/A",
      str_detect(str_to_lower(Exercise_frequency), "never") ~ "Never",
      str_detect(str_to_lower(Exercise_frequency), "rare") ~ "Rarely",
      str_detect(str_to_lower(Exercise_frequency), "daily") ~ "Daily",
      str_detect(str_to_lower(Exercise_frequency), "1-2|1_2") ~ "1-2_times_per_week",
      str_detect(str_to_lower(Exercise_frequency), "3-5|3_5") ~ "3-5_times_per_week",
      TRUE ~ "N/A"
    ),
    
    # Smoking_status - Specific parameters
    Smoking_status = case_when(
      Smoking_status == "N/A" ~ "N/A",
      str_detect(str_to_lower(Smoking_status), "non") ~ "Non-smoker",
      str_detect(str_to_lower(Smoking_status), "smok") ~ "Smoker",
      TRUE ~ "N/A"
    ),
    
    # Daily_cigarettes - Specific parameters
    Daily_cigarettes = case_when(
      Daily_cigarettes == "N/A" ~ "N/A",
      str_detect(Daily_cigarettes, "1[-–]5|1_5") ~ "1-5",
      str_detect(Daily_cigarettes, "6[-–]10|6_10") ~ "6-10",
      str_detect(Daily_cigarettes, "11[-–]15|11_15") ~ "11-15",
      str_detect(Daily_cigarettes, "16[-–]20|16_20") ~ "16-20",
      str_detect(Daily_cigarettes, "20\\+|\\+20") ~ "+20",
      TRUE ~ "N/A"
    ),
    
    # Notes_Samples - First letter uppercase with underscores
    Notes_Samples = clean_text(Notes_Samples)
  )
# ─────────────────────────────────────────────────────────────────────────────
# 6) Handle run-specific columns if they exist -------------------------------
# ─────────────────────────────────────────────────────────────────────────────

# Check if run-specific columns exist and format them
if ("Sequencing_Platform" %in% colnames(metadata)) {
  metadata <- metadata %>%
    mutate(
      Sequencing_Platform = case_when(
        Sequencing_Platform == "N/A" ~ "N/A",
        str_detect(str_to_lower(Sequencing_Platform), "miseq") ~ "MiSeq",
        str_detect(str_to_lower(Sequencing_Platform), "nextseq") ~ "NextSeq",
        str_detect(str_to_lower(Sequencing_Platform), "grindion") ~ "GrindION",
        str_detect(str_to_lower(Sequencing_Platform), "minion") ~ "MinION",
        str_detect(str_to_lower(Sequencing_Platform), "ion.*genestudio") ~ "Ion_GeneStudio_S5",
        TRUE ~ "N/A"
      ),
      
      Sequencing_Type = case_when(
        Sequencing_Type == "N/A" ~ "N/A",
        str_detect(str_to_lower(Sequencing_Type), "16s") ~ "16S_rRNA",
        str_detect(str_to_lower(Sequencing_Type), "wgs") ~ "WGS",
        str_detect(str_to_lower(Sequencing_Type), "shotgun") ~ "Shotgun_metagenomic",
        str_detect(str_to_lower(Sequencing_Type), "rna") ~ "RNA-Seq",
        str_detect(str_to_lower(Sequencing_Type), "long") ~ "Long-read",
        str_detect(str_to_lower(Sequencing_Type), "metatrans") ~ "Metatranscriptomics",
        TRUE ~ "N/A"
      ),
      
      Technician_name = clean_text(Technician_name),
      Notes_Runs = clean_text(Notes_Runs)
    )
}



# ─────────────────────────────────────────────────────────────────────────────
# 8) Write out the cleaned CSV ------------------------------------------------
# ─────────────────────────────────────────────────────────────────────────────
print("Writing cleaned metadata...")
write_csv(metadata, "controls/metadata_controls.csv", na = "")

print("✓ Successfully created metadata_cleaned.csv")
print(paste("Final dataset contains", nrow(metadata), "rows and", ncol(metadata), "columns"))

# Display summary
print("\n=== SUMMARY ===")
print(paste("Sample metadata entries:", nrow(metadata_sample)))
print(paste("Run metadata entries:", nrow(metadata_run))) 
print(paste("Merged entries:", nrow(metadata)))
print("\nColumn names in final dataset:")
print(colnames(metadata))

# Display sample of key formatted columns
print("\n=== SAMPLE OF FORMATTED DATA ===")
sample_cols <- c("Sample_Accession", "Gender", "Age", "Body_Mass_Index", "Dietary_Information", "Allergies")
available_cols <- intersect(sample_cols, colnames(metadata_clean))
if (length(available_cols) > 0) {
  print(head(metadata_clean[, available_cols, drop = FALSE], 10))
}



