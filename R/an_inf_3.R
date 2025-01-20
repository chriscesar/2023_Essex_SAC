# an_inf_3.R ####
# analyses of infaunal data
# now splitting data into WB chunks and analysing by WB_BSH

# Set up ####
## load packages ####
ld_pkgs <- c("tidyverse", "vegan","tictoc","ggtext","ggpmisc","mvabund")
vapply(ld_pkgs, library, logical(1L),
       character.only = TRUE, logical.return = TRUE)
rm(ld_pkgs)
tictoc::tic.clearlog()

## load data & prep for analysis ####
tic("load metadata")
source("R/set_meta.R")
source("R/imp_data_inf.R")

# deal with Presence-only taxa ####
#create new list element to store tweaked data
dfw$abundance_raw <- dfw$abundance

## convert abundance values of -999 to 1
dfw$abundance[dfw$abundance <0] <- 1
toc(log = TRUE)

tic("Data analysis, plots & export")
# Ensure the abundance data matches metadata rows
dfw$metadata <- dfw$metadata %>% mutate(id = row_number())
dfw$abundance <- dfw$abundance %>% mutate(id = row_number())

# FOR LOOPS MDS ####
dfw$metadata <- dfw$metadata %>% mutate(id = row_number())
dfw$abundance <- dfw$abundance %>% mutate(id = row_number())

full_data <- dfw$metadata %>%
  inner_join(dfw$abundance, by = "id")

# Step 1: Merge metadata and abundance
full_data <- dfw$metadata %>%
  inner_join(dfw$abundance, by = "id")

# Step 2: Remove unneeded columns and filter levels of PSA
## remove uneeded cols
full_data %>%
  dplyr::select(.,-c(Waterbody_Year,Site,Replicate:Biotope_EUNIS,id)
  ) %>%
  #remove A5.1 and A5.2 due to limited presence across WBs
  filter(.,
         PSA != "A5.1",
         PSA != "A5.2"#,
         # Waterbody != "Blackwater",
         # Waterbody != "Colne"
         ) -> trim_data

# Step 3: Get unique levels of Waterbody and PSA
waterbody_levels <- unique(trim_data$Waterbody)
psa_levels <- unique(trim_data$PSA)

# Step 4: Create lists to store results
mds_results <- list()
mvabund_results <- list()
summary_list <- list()

# Step 5: Loop through Waterbody and PSA levels
for (waterbody in waterbody_levels) {
  for (psa in psa_levels) {
    # Filter the data for the current combination ####
    trim_data2 <- trim_data %>%
      filter(!(Waterbody=="Blackwater" & PSA == "A5.3")) %>% 
      filter(!(Waterbody=="Colne" & PSA == "A5.4"))
    subset_data <- trim_data2 %>%
      filter(Waterbody == waterbody, PSA == psa)
      # filter(Waterbody == "Crouch", PSA == "A5.3")
    
    # remove 'empty columns'
    subset_data %>%
      dplyr::select(1:4, where(~ is.numeric(.x) && sum(.x) != 0)) -> subset_data

    # Skip if there are not enough rows for ordination
    if (nrow(subset_data) < 2) {
      message(paste("Skipping Waterbody:", waterbody, "and PSA:", psa, "due to insufficient data."))
      next
    }
    
    # Convert to matrix for metaMDS (removing grouping variables)
    abundance_matrix <- as.matrix(subset_data %>%
                                    dplyr::select(-Waterbody, -PSA, -Site_ID,-year))
    
    ## Run metaMDS and store results ####
    mds_result <- metaMDS(abundance_matrix, distance = "bray", trymax = 500, try = 200)
    
    # Save the result in the list
    mds_results[[paste(waterbody, psa, sep = "_")]] <- list(
      Waterbody = waterbody,
      PSA = psa,
      mds = mds_result
    )
    
    ## output mds data ####
    mds_scores <- as_tibble(as.data.frame(scores(mds_result,"site")))
    mds_scores$PSA <- psa
    mds_scores$WB <- waterbody
    mds_scores$year <- subset_data$year
    mds_scores$stress <- mds_result$stress
    
    spp_scores <- as_tibble(as.data.frame(scores(mds_result,"species")))
    spp_scores$species <- colnames(subset_data)[-c(1:4)]
    spp_scores$species_sh <- vegan::make.cepnames(spp_scores$species)

    ## plot ####
    png(file=paste0("figs/",waterbody,"_",psa,".png"),
        width=12*ppi, height=6*ppi, res=ppi)
    
    mds_scores %>%
      ggplot(.,aes(x=NMDS1, y=NMDS2))+
      geom_text(data=spp_scores,
                aes(
                  x=NMDS1,
                  y= NMDS2,
                  label=species_sh
                ),
                col="darkgrey",
                size=4, show.legend = FALSE, alpha = 0.5)+
      geom_text(aes(
        label=year,
        colour=as.factor(year)
      ),
      show.legend = FALSE,
      size=7,
      fontface=2)+
      geom_text_npc(aes(npcx = .99, npcy = .99,
                        label=paste("Stress = ",
                                    round(mds_scores$stress[1], 3))))+
      labs(
        title = "MDS ordination of benthic infauna",
        subtitle = paste0(waterbody," ",psa," broadscale habitat"),
        colour = "Year")+
      coord_fixed()+
      scale_colour_manual(
        values = c("2012" = "#009E73", "2017" = "#0C7BDC", "2023" = "#d41159")
      ) +
      theme(
        axis.title = element_text(face=2)
      ) -> pl
    print(pl)
    dev.off()
    
    ## mvabund comparisons ####
    # generate offset for number of replicates
    subset_data %>% 
      group_by(year, Site_ID) %>% 
      summarise(reps = n(), .groups = "drop") -> offset
    subset_mvab <- left_join(subset_data, offset, by=c("Site_ID","year"))
    
    # Convert to mvabund object for analysis for metaMDS (removing grouping variables)
    abundance_mvab <- mvabund::mvabund(subset_mvab %>%
                                         dplyr::select(-reps,
                                                       -Waterbody,
                                                       -PSA,
                                                       -Site_ID,
                                                       -year
                                                       ))
    # Run manyglm and store results
    mvabund_result <- manyglm(abundance_mvab ~ subset_mvab$year,
                              offset = log(subset_mvab$reps),
                              family = "negative.binomial")
    # Save the result in the list
    mvabund_results[[paste(waterbody, psa, sep = "_")]] <- list(
      Waterbody = waterbody,
      PSA = psa,
      mvab = mvabund_result
    )
    rm(pl,mds_scores,mds_result,spp_scores,abundance_matrix)
  }
}

# save outputs
saveRDS(mvabund_results,file="data_out/mvabund_results_WB_BSH.Rdat")
saveRDS(mds_results,file="data_out/mds_results_WB_BSH.Rdat")
rm(abundance_mvab,offset,subset_mvab,subset_data,psa,waterbody)
toc(log=TRUE)

tic("summarise & export MVABUND results")
## summarise & export MVABUND results ####
# Loop through mvabund_results
for (key in names(mvabund_results)) {
  # Extract the mvabund_result object
  mvab_result <- mvabund_results[[key]]$mvab
  
  # Generate a summary of the result
  summary_output <- summary(mvab_result)
  
  # Save the summary in the list
  summary_list[[key]] <- summary_output
  
  # Define the file name for the text output
  output_file <- paste0("data_out/mvabSummary_", key, ".txt")
  
  # Save the summary to a text file
  capture.output(summary_output, file = output_file)
  
  message(paste("Summary for", key, "saved to", output_file))
}

# Save the summary list to an R object file
saveRDS(summary_list, file = "data_out/mvabund_summaries.rds")
summary_list <- readRDS("data_out/mvabund_summaries.rds")
toc(log=TRUE)

### generate plots to find significant taxa
tic("Anova on mvabund models")
# create holding pen
anova_out <- list()
# for(key in names(mvabund_results)){
#   # Extract the mvabund_result object
#   mvab_result <- mvabund_results[[key]]$mvab
#   fit.glm.out <- mvabund::anova.manyglm(mvab_result,p.uni = "adjusted", test="LR",show.time="all")
#   m2tmp1 <- t(as.data.frame(fit.glm.out$uni.p))[,2]
# 
#   nsig <- length(names(m2tmp1[m2tmp1<0.1]))
#   print(paste0(nsig," 'significant' taxa"))
# 
#   ##save anova outputs
#   anova_out[[key]] <- fit.glm.out
#   }
# 
# saveRDS(anova_out, file="data_out/mvabund_anova_out_WB_BSH.Rdat")
anova_out <- readRDS("data_out/mvabund_anova_out_WB_BSH.Rdat")
toc(log=TRUE)

# generate plots of significant taxa ####
tic("generate plots of significant taxa")

for(i in names(anova_out)){
  WB <-  stringr::str_split(i,"_",simplify = TRUE)[1]
  psa <- stringr::str_split(i,"_",simplify = TRUE)[2]
  x <- anova_out[[i]] ## anova outputs
  df_tmp <- trim_data %>% filter(.,Waterbody==WB & PSA==psa)
  m2tmp1 <- t(as.data.frame(x$uni.p))[,2]
  nm <- names(m2tmp1)
  nm <- nm %>% str_replace_all("\\."," ")
  names(m2tmp1) <- nm; rm(nm)
  nsig <- length(names(m2tmp1[m2tmp1<0.1]))
  print(paste0(i,": ",nsig," 'significant' taxa"))
  if(nsig == 0){
    print(paste0(i,": NO SIGNIFICANT TAXA"))
  } else{
  m2tx <- names(m2tmp1[m2tmp1 < 0.1]) # Which taxa are 'significantly' different?
  kptx <- names(df_tmp) %in% m2tx
  m2tx <- df_tmp[, kptx]
  m2tx$Year <- df_tmp$year
  # Make long
  m2txl <- m2tx %>%
    relocate(Year) %>% # Move Year to start
    pivot_longer(cols = c(2:ncol(.)))
  m2txl$Year <- as.factor(m2txl$Year)
    # plot
  ggplot(m2txl,aes(x=log(value+1),y=name,
                   fill=Year,
                   shape = Year))+
    geom_jitter(data=m2txl[,c(2:3)], inherit.aes = FALSE,
                aes(x=log(value+1),y=name,),
                height = 0.05,size=1, alpha = 0.5, colour = "grey") +
    geom_jitter(height = 0.05,size=3, alpha = 0.9) +
    scale_shape_manual(values = c(21:24))+
    labs(title = paste0(WB, " ", psa," BSH"),
         x="log(Taxon abundance (n+1))",
         caption=paste0("Displayed taxa are the ",paste0(length(names(m2tmp1[m2tmp1<0.1])))," taxa which showed differences in abundances between years at alpha = 0.1.\n",
                        "Taxon abundances across all years are presented in each facet, with abundances for a given year displayed by larger, coloured icons."))+
    scale_fill_manual(values = cbPalette)+ scale_colour_manual(values = cbPalette)+
    facet_wrap(.~Year)+
    scale_y_discrete(limits=rev)+
    theme(
      legend.position = "none",
      axis.title.y = element_blank(),
      axis.text.y = element_text(size=12,face="italic"),
      axis.title.x = element_text(face="bold"),
      strip.text.x = element_text(face="bold",size=12),
      plot.title.position = "plot",
      plot.title = element_text(face="bold",size=14)
      ) -> pl3
  ggsave(filename = paste0("figs/infauna_",i,"_relabund_WB_BSH.png"),
           width=12, height=6,plot=pl3);rm(pl3)
  rm(m2tx,kptx,m2txl)
  }
}
toc(log=TRUE)
