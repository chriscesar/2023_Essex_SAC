# NtaxSummary.R ####

# Set up ####
## load packages ####
ld_pkgs <- c("tidyverse", "tictoc")
vapply(ld_pkgs, library, logical(1L),
       character.only = TRUE, logical.return = TRUE)
rm(ld_pkgs)
tictoc::tic.clearlog()

## load data & prep for analysis ####
tic("load metadata")
source("R/set_meta.R")
source("R/imp_data_inf.R")

mt <- names(dfw$metadata %>% dplyr::select(.,Waterbody_Year, PSA))
dfw$metadata %>% dplyr::select(.,year, PSA) -> dfmt

### reduce to 1 dataframe for grouping
df_tmp <- as_tibble(cbind(dfmt,dfw$abundance))

## calculate means sp WB & year
df_tmp %>% 
  pivot_longer(., cols = -c(
    year,
    PSA)) %>% 
  group_by(year,PSA,name) %>% 
  summarise(.,mean = mean(value),.groups = "drop") %>% ungroup() %>% 
  pivot_wider(.,names_from = name, values_from = mean) -> df_out

write.csv(df_out, file = "data_out/inf_tmp.csv",row.names = FALSE)
