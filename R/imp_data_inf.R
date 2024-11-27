# imp_data_inf.R
# import infauna data and convert to format for analysis

# Set up ####
### load packages ####
ld_pkgs <- c("tidyverse","readxl", "tictoc")
vapply(ld_pkgs, library, logical(1L),
       character.only = TRUE, logical.return = TRUE)
rm(ld_pkgs)
tictoc::tic.clearlog()

### load metadata ####
tic("load metadata")
source("R/set_meta.R")
toc(log = TRUE)

# copy & load data ####
tic("copy & load data")

# Set the source folder
source_folder <- paste0(fol, "Data/Infauna/")
# Set the destination folder
destination_folder <- "data_in/"
file_names <- c("2023 Infauna data for R_USE.xlsx",
                "2023_Infauna_Working_USE.xlsx")
# Iterate over each file name and copy if it doesn't exist in the destination
for (file_name in file_names) {
  source_file <- paste0(source_folder, file_name)
  destination_file <- paste0(destination_folder, file_name)
  
  # Check if the file exists in the destination folder
  if (!file.exists(destination_file)) {
    # Copy the file from the source folder to the destination folder
    if (file.copy(source_file, destination_file)) {
      cat("File", file_name, "copied successfully.\n")
    } else {
      cat("Failed to copy file:", file_name, "\n")
    }
  } else {
    cat("File", file_name, "already exists in the destination folder.\n")
  }
}

df0_2023 <- as_tibble(read_xlsx("data_in/2023 Infauna data for R_USE.xlsx",
                           sheet = "2023 whole site",
                           guess_max = 10000))
df0_2017 <- as_tibble(read_xlsx("data_in/2023 Infauna data for R_USE.xlsx",
                                sheet = "2017 whole site",
                                guess_max = 10000))
df0_2012 <- as_tibble(read_xlsx("data_in/2023 Infauna data for R_USE.xlsx",
                                sheet = "2012 Colne",
                                guess_max = 10000))

df_tax <- as_tibble(read_xlsx("data_in/2023_Infauna_Working_USE.xlsx",
                              sheet = "ALL WoRMS Match",
                              guess_max = 10000))

toc(log = TRUE)
