
# Load necessary libraries
library(readr)
library(dplyr)
library(stringr)

##########################################
#Install libraries
##########################################
if (!require("tidyverse"))    install.packages("tidyverse")
if (!require("lubridate"))    install.packages("lubridate")
if (!require("janitor"))      install.packages("janitor")

library(tidyverse)
library(lubridate)
library(janitor)
# ─────────────────────────────────────────────────────────────────────────────
# 2) Read the template CSV ---------------------------------------------------
# ─────────────────────────────────────────────────────────────────────────────
# (Adjust the path if needed)
metadata <- read_csv("../metadata/template_metadata_merged/metadata_merged.csv")


# ─────────────────────────────────────────────────────────────────────────────
# 3) Clean up any new column names -------------------------------------------
#    - Trim whitespace
#    - Replace spaces or non-alphanumeric chars with "_"
#    - Preserve existing uppercase letters & digits
# ─────────────────────────────────────────────────────────────────────────────
clean_names_safe <- function(names_vec) {
  names_vec %>%
    str_trim() %>%
    str_replace_all("[^A-Za-z0-9]+", "_")
}

colnames(metadata) <- clean_names_safe(colnames(metadata))


# ─────────────────────────────────────────────────────────────────────────────
# 4) Fill only empty cells with "N/A" ----------------------------------------
#    - For character columns: if "" or NA, replace with "N/A"
#    - Numeric columns untouched
# ─────────────────────────────────────────────────────────────────────────────
metadata <- metadata %>%
  mutate(across(
    .cols = where(is.character),
    .fns  = ~ if_else(is.na(.x) | .x == "", "N/A", .x)
  ))

# ─────────────────────────────────────────────────────────────────────────────
# 4.2) Trim and replace spaces in all text fields ----------------------------
#    - Trim leading/trailing whitespace
#    - Replace any internal spaces with "_"
# ─────────────────────────────────────────────────────────────────────────────
metadata <- metadata %>%
  mutate(across(
    .cols = where(is.character),
    .fns  = ~ str_trim(.x) %>% str_replace_all("\\s+", "_")
  ))

# ─────────────────────────────────────────────────────────────────────────────
# 4.3) Apply Title Case to selected columns ----------------------------------
#    - Institution, Department, Analyst_Processor_Name, Notes_Samples
#    - Convert "_" back to space, title-case, then "_" again
# ─────────────────────────────────────────────────────────────────────────────
title_cols <- c("Institution",
                "Department",
                "Analyst_Processor_Name",
                "Notes_Samples")

metadata <- metadata %>%
  mutate(across(
    .cols = any_of(title_cols),
    .fns  = ~ .x %>%
      str_replace_all("_", " ") %>%
      str_to_title() %>%
      str_replace_all("\\s+", "_")
  ))

# ─────────────────────────────────────────────────────────────────────────────
# 5) Field-specific formatting ------------------------------------------------
#    - Collection_Date: DD/MM/YYYY
#    - Collection_Storage_Temperature: enforce allowed values
#    - Gender: Female, Male, Not_specified
#    - Age: integer or "N/A"
#    - (Extendable: other case_when for controlled vocabularies)
# ─────────────────────────────────────────────────────────────────────────────
metadata <- metadata %>%
  mutate(
    # Format dates as DD/MM/YYYY
    Collection_Date = dmy(Collection_Date) %>% format("%d/%m/%Y"),
    
    # Ensure storage temperature matches one of: -80, -20, 4, Room_Temperature
    Collection_Storage_Temperature = case_when(
      str_detect(Collection_Storage_Temperature, "^-?80$")   ~ "-80",
      str_detect(Collection_Storage_Temperature, "^-?20$")   ~ "-20",
      Collection_Storage_Temperature == "4"                 ~ "4",
      tolower(Collection_Storage_Temperature) %in% 
        c("room_temperature", "room temperature")           ~ "Room_Temperature",
      TRUE                                                  ~ Collection_Storage_Temperature
    ),
    
    # Normalize gender values
    Gender = case_when(
      tolower(Gender) == "female"        ~ "Female",
      tolower(Gender) == "male"          ~ "Male",
      TRUE                                ~ "Not_specified"
    ),
    
    # Force Age to integer; if invalid, set to "N/A"
    Age = suppressWarnings(as.integer(Age)),
    Age = if_else(is.na(Age), "N/A", as.character(Age))
    
    # Add more case_when() here for:
    # Ongoing_Conditions, Allergies, Dietary_Information, etc.
    # Example:
    # , Allergies = case_when(
    #     Allergies %in% c("Peanuts","Shellfish",...) ~ Allergies,
    #     TRUE                                       ~ "Other"
    #   )
  )

# ─────────────────────────────────────────────────────────────────────────────
# 6) Write out the cleaned CSV ------------------------------------------------
# ─────────────────────────────────────────────────────────────────────────────
write_csv(metadata, "metadata_cleaned.csv")

  
