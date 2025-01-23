# psa_v2.R
# Import PSA data and generate ternary plots
# amended to generate plots by WB


## load required packages ####
ld_pkgs <- c("grid","tidyverse","tictoc","ggtern") # what packages do we need to load?
vapply(ld_pkgs, library, logical(1L), # load them and display TRUE/FALSE if loaded
       character.only = TRUE, logical.return = TRUE);rm(ld_pkgs)
tictoc::tic.clearlog()

## load metadata
tic("load metadata")
source("R/set_meta.R")
df_psa <- as_tibble(readxl::read_xlsx(paste0(fol,"Data/PSA/2023_Essex_PSA.xlsx"),sheet = "All PSA for R"))

# assign Water body values
df_psa %>% mutate(WB = case_when(
  startsWith(Station, "BLK") ~ "Blackwater",
  startsWith(Station, "COL") ~ "Colne",
  startsWith(Station, "CRO") ~ "Crouch",
  startsWith(Station, "ESX") ~ "Essex",
  startsWith(Station, "OBW") ~ "Blackwater Outer",
  TRUE ~ NA_character_ # Default value if no condition is met
  )) %>% relocate(WB) -> df_psa


###### TO DO: STANDARDISE CHART COLOURS TO MATCH THOSE IN MDS PLOTS
# Step 3: Get unique levels of Waterbody and PSA
waterbody_levels <- unique(df_psa$WB)
psa_levels <- unique(df_psa$Eunis_Code)

# Step 4: Create a list to store results
psa_results <- list()

for (waterbody in waterbody_levels) {
  for (psa in psa_levels) {
    bsh_nm <- BSH_codes$BSH_name[BSH_codes$BSH_code==psa]
    # Filter the data for the current combination
    subset_data <- df_psa %>%
      filter(WB == waterbody, Eunis_Code == psa)
    
    # Skip if there are not enough rows for ordination
    if (nrow(subset_data) < 1) {
      message(paste("Skipping Waterbody:", waterbody, "and PSA:", psa, "due to insufficient data."))
      next
    }
    
    tic(paste0(waterbody,"_",psa, " prep data"))
    ggplot(data=subset_data,
           aes(Silt,
               Gravel,
               Sand,
               colour=as.factor(Year))) +
      coord_tern()+
      theme_bw() +
      theme_showarrows() + custom_percent('%') + 
      geom_mask() +
      geom_text(aes(label=Station),show.legend = T,
              fontface=2) + 
      labs(Tarrow = "Gravel",
         Larrow = "Silt",
         Rarrow = "Sand")+
      #scale_colour_manual("Year", values=cbPalette)+
      scale_colour_manual(
      name = "Year",
      values = c("2012" = "#009E73", "2017" = "#0C7BDC", "2023" = "#d41159")
      ) +
      labs(title = paste0(waterbody,": ",psa," ",bsh_nm," BSH"),
         caption = "Sample identities coloured by sample year")+
      # geom_line(aes(group = MatchedSite),colour="darkgrey") +# Add this line to connect points
      theme(axis.title = element_text(face=2),
          legend.title = element_text(face=2),
          plot.title = element_text(hjust=0.5,vjust=0,face="bold",margin=margin(0,0,-20,0)),
          legend.position = c(.9, .5))
    ggsave(filename = paste0("figs/PSA_.",psa,"_",waterbody,".WB.png"),
         device = "png",
         width=10, height=10, units = "in")
  toc(log=TRUE)
  flush.console()
  }
}

toc(log=TRUE)

### house keeping

rm(list=ls(pattern = "^psa"))
rm(list=ls(pattern = "^waterbody"))
rm(list=ls(pattern = "^bsh"))
rm(df,BSH_codes,psa_results,subset_data, cbPalette,fol,nit,perm,ppi)
detach("package:ggtern",  unload=TRUE)
detach("package:tictoc",unload = TRUE)
detach("package:grid",unload = TRUE)
detach("package:tidyverse",unload = TRUE)
