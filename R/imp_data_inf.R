# imp_data_inf.R
# import infauna data and convert to wide format for analysis

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

rm(source_file, source_folder,
   destination_file, destination_folder,
   file_name,file_names)
toc(log = TRUE)

# convert to long format & join taxon info ####
tic("convert to long format for joining")
# Define a function to pivot each dataframe
reshape_to_long <- function(df, id_cols, value_name_prefix) {
  df %>%
    pivot_longer(
      cols = -all_of(id_cols), # All columns except the first 14
      names_to = "taxon",   # Name of the new column containing the original column names
      values_to = "abundance" # Values column name with prefix
    )
}

# Specify the shared columns (first 14)
id_columns <- names(df0_2012)[1:14]

# Reshape each dataframe
df_long_2012 <- reshape_to_long(df0_2012, id_columns, "2012")
df_long_2017 <- reshape_to_long(df0_2017, id_columns, "2017")
df_long_2023 <- reshape_to_long(df0_2023, id_columns, "2023")

rm(reshape_to_long)

# Combine them into a single dataframe (optional)
df_long <- bind_rows(
  df_long_2012 %>% mutate(year = 2012),
  df_long_2017 %>% mutate(year = 2017),
  df_long_2023 %>% mutate(year = 2023)
)
rm(df0_2012,df0_2017,df0_2023,id_columns,
   df_long_2012,df_long_2017,df_long_2023)

## append taxon info to long data & format ####
df_long %>%
  ## remove '0' values
  filter(.,abundance !=0) %>% 
  left_join(., df_tax, by = c("taxon" = "AcceptedName_Qualifier")) %>% 
  ### remove variables flagged as 'remove' (and retain NA values)
  filter(.,!str_starts(Flag, "Remove")|is.na(Flag)) %>% 
  ### remove unnecessary columns
  dplyr::select(., c(Waterbody_Year:year,ScientificName_accepted)) %>% 
  dplyr::select(., -taxon) %>% 
  #rename 'ScientificName_accepted' to 'taxon'
  rename(taxon = ScientificName_accepted) %>% #names()
  ### group_by everything except abundance & sum 'duplicated' rows
  group_by(across(c(!abundance))) %>% 
  summarise(.,abundance=sum(abundance),.groups = "drop") %>%
  ### widen data (species as columns)
  pivot_wider(.,names_from = taxon,
              values_from = abundance,
              values_fill = 0) -> dfw

write.csv(dfw,file="data_out/inf_wide_all.csv", row.names = FALSE)
rm(df_tax)
toc(log=TRUE)

## convert to list ####
tic("convert to list")
dfwmeta <- dfw[,c(1:15)]
dfwspp <- dfw[,c(16:length(dfw))]
dfw <- list(dfwmeta, dfwspp)
rm(dfwmeta,dfwspp)
names(dfw) <- c("metadata", "abundance")

unlist(tictoc::tic.log())
