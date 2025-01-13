# psa.R
# Import PSA data and generate ternary plots


## load required packages ####
ld_pkgs <- c("grid","tidyverse","tictoc","ggtern") # what packages do we need to load?
vapply(ld_pkgs, library, logical(1L), # load them and display TRUE/FALSE if loaded
       character.only = TRUE, logical.return = TRUE);rm(ld_pkgs)
tictoc::tic.clearlog()

## load metadata
tic("load metadata")
source("R/set_meta.R")
df_psa <- as_tibble(read_xlsx(paste0(fol,"Data/PSA/2023_Essex_PSA.xlsx"),sheet = "All PSA for R"))

###### TO DO: STANDARDISE CHART COLOURS TO MATCH THOSE IN MDS PLOTS



for (bshcode in unique(df_psa$Eunis_Code)) {
  # subset data by BSH
  #bsh_data <- subset(df_psa, Eunis_Code == "A5.4")
  bsh_data <- subset(df_psa, Eunis_Code == bshcode)
  
  
  tic(paste0(unique(bsh_data$Eunis_Code)[1], " prep data"))
  
  ggplot(data=bsh_data,
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
    labs(title = paste0(unique(bsh_data$Eunis_Code)[1]," BSH"),
         caption = "Sample identities coloured by sample year")+
    # geom_line(aes(group = MatchedSite),colour="darkgrey") +# Add this line to connect points
    theme(axis.title = element_text(face=2),
          legend.title = element_text(face=2),
          plot.title = element_text(hjust=0.5,vjust=0,face="bold",margin=margin(0,0,-20,0)),
          legend.position = c(.9, .5))
  
  ggsave(filename = paste0("figs/PSA_BSH_.",unique(bsh_data$Eunis_Code)[1],".png"),
         device = "png",
         width=10, height=10, units = "in")
  toc(log=TRUE)
  flush.console()
}

    