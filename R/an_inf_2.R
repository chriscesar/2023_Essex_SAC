# an_inf_2.R ####
# analyses of infaunal data
# now splitting data into WB chunks

# Set up ####
## load packages ####
ld_pkgs <- c("tidyverse", "vegan","lmerTest","rstatix", "mvabund","tictoc",
             "MASS","ggtext","ggpmisc", "gllvm","MASS","lme4")
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
toc(log = TRUE);toc(log = TRUE)

# split by WB
# BW_tmp <- dfw$metadata$Waterbody == "Blackwater"
# CO_tmp <- dfw$metadata$Waterbody == "Colne"
# CR_tmp <- dfw$metadata$Waterbody == "Crouch"
# ES_tmp <- dfw$metadata$Waterbody == "Essex"
# BO_tmp <- dfw$metadata$Waterbody == "OuterBlackwater"

# Ensure the abundance data matches metadata rows
dfw$metadata <- dfw$metadata %>% mutate(id = row_number())
dfw$abundance <- dfw$abundance %>% mutate(id = row_number())

# # Merge metadata and abundance
# full_data <- dfw$metadata %>%
#   inner_join(dfw$abundance, by = "id")
# 
# ## remove uneeded cols
# full_data %>% 
#   dplyr::select(.,-c(Waterbody_Year,Site,Replicate:Biotope_EUNIS,id)
#                 ) %>% 
#   #remove A5.1 and A5.2 due to limited presence across WBs
#   filter(.,
#          PSA != "A5.1",
#          PSA != "A5.2"
#   ) -> trim_data
# 
# # Nest data by WaterBody and BSH
# nested_data <- trim_data %>%
#   group_by(Waterbody, PSA) %>%
#   nest()
# 
# # Define a function for metaMDS
# run_metaMDS <- function(data) {
#   abundance_matrix <- as.matrix(data[ , !(names(data) %in% c("WaterBody", "BSH", "id"))])
#   metaMDS(abundance_matrix[,-c(1:4)], distance = "bray", trymax = 500, try=200)
# }
# 
# # Apply the function to each group
# nested_results <- nested_data %>%
#   mutate(mds_result = map(data, run_metaMDS))
# 

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
         PSA != "A5.2",
         Waterbody != "Blackwater",
         Waterbody != "Colne") -> trim_data


# Step 3: Get unique levels of Waterbody and PSA
waterbody_levels <- unique(trim_data$Waterbody)
psa_levels <- unique(trim_data$PSA)

# Step 4: Create a list to store results
mds_results <- list()

# Step 5: Loop through Waterbody and PSA levels
for (waterbody in waterbody_levels) {
  for (psa in psa_levels) {
    # Filter the data for the current combination
    subset_data <- trim_data %>%
      filter(Waterbody == waterbody, PSA == psa)
    
    # Skip if there are not enough rows for ordination
    if (nrow(subset_data) < 2) {
      message(paste("Skipping Waterbody:", waterbody, "and PSA:", psa, "due to insufficient data."))
      next
    }
    
    # Convert to matrix for metaMDS (removing grouping variables)
    abundance_matrix <- as.matrix(subset_data %>%
                                    dplyr::select(-Waterbody, -PSA, -Site_ID,-year))
    
    # Run metaMDS and store results
    mds_result <- metaMDS(abundance_matrix, distance = "bray", trymax = 500, try = 200)
    
    # Save the result in the list
    mds_results[[paste(waterbody, psa, sep = "_")]] <- list(
      Waterbody = waterbody,
      PSA = psa,
      mds = mds_result
    )
  }
}

# Step 6: Extract Results
# For example, stress values for all combinations
stress_values <- lapply(mds_results, function(x) x$mds$stress)

# Display stress values
stress_values



############################################
# MDS plots by BSH####
## A5.2 ####
tic("A5.2")
# which rows correspond to BSH?
tmp_row <- dfw$metadata$PSA == "A5.2"
# keep only those data
tmp_met <- dfw$metadata[tmp_row,]
tmp_abd <- dfw$abundance[tmp_row,]
rm(tmp_row)

### run ordination ####
tic("run ordination")
tmpord <- metaMDS(tmp_abd,try = 200,trymax = 500)
toc(log = TRUE)

tic("Extract ordination data for plotting using ggplot")
## extract Site scores
mds_scores <- as_tibble(as.data.frame(scores(tmpord,"site")))

mds_scores$Biotope <- tmp_met$Biotope
mds_scores$Biotope_EUNIS <- tmp_met$Biotope_EUNIS
mds_scores$Waterbody <- tmp_met$Waterbody
mds_scores$Year <- tmp_met$year
mds_scores$Site <- tmp_met$Site_ID
mds_scores$PSA <- tmp_met$PSA
mds_scores$Site_USE <- tmp_met$Site_USE

### extract species scores and groups
spp_scores <- as.data.frame(scores(tmpord, "species"))  #Using the scores function from vegan to extract the species scores and convert to a data.frame
spp_scores$species <- rownames(spp_scores)  # create a column of species, from the rownames of species.scores
spp_scores$species_sh <- make.cepnames(spp_scores$species)

### plot ####
pdf("figs/inf_mds_A5.2.pdf",
    # width=14,height = 14
    width=12,height = 12
)
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
    label=Site_USE,
    colour=as.factor(Year)
  ),
  size=7,
  fontface=2)+
  geom_text_npc(aes(npcx = .99, npcy = .99,
                    label=paste("Stress = ",
                                round(tmpord$stress, 3))))+
  labs(colour = "Year")+
  coord_fixed()+
  scale_colour_manual(
    values = c("2012" = "#009E73", "2017" = "#0C7BDC", "2023" = "#d41159")
  ) +
  theme(
    legend.title = element_text(face=2),
    legend.text = element_text(face=2),
    legend.position = c(.99, 0.01),
    legend.justification = c(1, 0),
    # legend.direction = "horizontal",
    legend.box = "vertical",
    axis.title = element_text(face=2)
  )
dev.off()
rm(tmp_row,tmp_abd,tmp_met,tmpord,mds_scores,spp_scores)

toc(log = TRUE);toc(log = TRUE);toc(log = TRUE)
##

## A5.3 ####
tic("A5.3")
# which rows correspond to BSH?
tmp_row <- dfw$metadata$PSA == "A5.3"
# keep only those data
tmp_met <- dfw$metadata[tmp_row,]
tmp_abd <- dfw$abundance[tmp_row,]
rm(tmp_row)

### run ordination ####
tic("run ordination")
tmpord <- metaMDS(tmp_abd,try = 200,trymax = 500)
toc(log = TRUE)

tic("Extract ordination data for plotting using ggplot")
## extract Site scores
mds_scores <- as_tibble(as.data.frame(scores(tmpord,"site")))

mds_scores$Biotope <- tmp_met$Biotope
mds_scores$Biotope_EUNIS <- tmp_met$Biotope_EUNIS
mds_scores$Waterbody <- tmp_met$Waterbody
mds_scores$Year <- tmp_met$year
mds_scores$Site <- tmp_met$Site_ID
mds_scores$PSA <- tmp_met$PSA
mds_scores$Site_USE <- tmp_met$Site_USE

### extract species scores and groups
spp_scores <- as.data.frame(scores(tmpord, "species"))  #Using the scores function from vegan to extract the species scores and convert to a data.frame
spp_scores$species <- rownames(spp_scores)  # create a column of species, from the rownames of species.scores
spp_scores$species_sh <- make.cepnames(spp_scores$species)

### plot ####
pdf("figs/inf_mds_A5.3.pdf",
    # width=14,height = 14
    width=12,height = 12
)
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
    label=Site_USE,
    colour=as.factor(Year)
  ),
  size=7,
  fontface=2)+
  geom_text_npc(aes(npcx = .99, npcy = .99,
                    label=paste("Stress = ",
                                round(tmpord$stress, 3))))+
  labs(colour = "Year")+
  coord_fixed()+
  scale_colour_manual(
    values = c("2012" = "#009E73", "2017" = "#0C7BDC", "2023" = "#d41159")
  ) +
  theme(
    legend.title = element_text(face=2),
    legend.text = element_text(face=2),
    legend.position = c(.99, 0.01),
    legend.justification = c(1, 0),
    # legend.direction = "horizontal",
    legend.box = "vertical",
    axis.title = element_text(face=2)
  )
dev.off()
rm(tmp_row,tmp_abd,tmp_met,tmpord,mds_scores,spp_scores)

toc(log = TRUE);toc(log = TRUE);toc(log = TRUE)
##

## A5.4 ####
tic("A5.4")
# which rows correspond to BSH?
tmp_row <- dfw$metadata$PSA == "A5.4"
# keep only those data
tmp_met <- dfw$metadata[tmp_row,]
tmp_abd <- dfw$abundance[tmp_row,]
rm(tmp_row)

### run ordination ####
tic("run ordination")
tmpord <- metaMDS(tmp_abd,try = 200,trymax = 500)
toc(log = TRUE)

tic("Extract ordination data for plotting using ggplot")
## extract Site scores
mds_scores <- as_tibble(as.data.frame(scores(tmpord,"site")))

mds_scores$Biotope <- tmp_met$Biotope
mds_scores$Biotope_EUNIS <- tmp_met$Biotope_EUNIS
mds_scores$Waterbody <- tmp_met$Waterbody
mds_scores$Year <- tmp_met$year
mds_scores$Site <- tmp_met$Site_ID
mds_scores$PSA <- tmp_met$PSA
mds_scores$Site_USE <- tmp_met$Site_USE

### extract species scores and groups
spp_scores <- as.data.frame(scores(tmpord, "species"))  #Using the scores function from vegan to extract the species scores and convert to a data.frame
spp_scores$species <- rownames(spp_scores)  # create a column of species, from the rownames of species.scores
spp_scores$species_sh <- make.cepnames(spp_scores$species)

### plot ####
pdf("figs/inf_mds_A5.4.pdf",
    # width=14,height = 14
    width=12,height = 12
)
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
    label=Site_USE,
    colour=as.factor(Year)
  ),
  size=7,
  fontface=2)+
  geom_text_npc(aes(npcx = .99, npcy = .99,
                    label=paste("Stress = ",
                                round(tmpord$stress, 3))))+
  labs(colour = "Year")+
  coord_fixed()+
  scale_colour_manual(
    values = c("2012" = "#009E73", "2017" = "#0C7BDC", "2023" = "#d41159")
  ) +
  theme(
    legend.title = element_text(face=2),
    legend.text = element_text(face=2),
    legend.position = c(.99, 0.01),
    legend.justification = c(1, 0),
    # legend.direction = "horizontal",
    legend.box = "vertical",
    axis.title = element_text(face=2)
  )
dev.off()
rm(tmp_row,tmp_abd,tmp_met,tmpord, mds_scores,spp_scores)

toc(log = TRUE);toc(log = TRUE);toc(log = TRUE)

# MVABUNDS ####
## ver 1: no offsets ####
tic("ANALYSES by BSH: Run mvabund models by BSH")

### remove A5.1 ####
remove <- which(dfw$metadata$PSA == "A5.1")
dfw_trim <- lapply(dfw, function(dfw) dfw[-remove,])
rm(remove)

# ## run MVABUND analysis
# for (bshcode in unique(dfw_trim$metadata$PSA)) {
#   ### which rows to keep
#   kp <- which(dfw$metadata$PSA == bshcode)
#   # subset data by BSH
#   bsh_data <- lapply(dfw, function(dfw) dfw[kp,])
#   
#   ## remove 'empty' columns
#   bsh_data$abundance <- bsh_data$abundance %>%
#     dplyr::select_if(where( ~ !is.numeric(.) || sum(.) !=0)) ## remove numeric variables that sum to 0
#   
#   ### produce meanvar plot ####
#   tic(paste0(unique(bsh_data$metadata$PSA)[1], " create meanvar plot"))
#   png(file = paste0("figs/infMeanVar.",unique(bsh_data$metadata$PSA)[1],".png"),
#       width=12*ppi, height=6*ppi, res=ppi)
#   mvpl <- mvabund::meanvar.plot(mvabund(bsh_data$abundance),
#                                 xlab="Mean",
#                                 ylab="Variance",
#                                 table=TRUE)
#   
#   # Step 1: Find the minimum and maximum values
#   min_value <- min(mvpl[,2])
#   max_value <- max(mvpl[,2])
#   
#   min_order <- floor(log10(min_value))
#   max_order <- floor(log10(max_value))
#   orders_of_magnitude_covered <- max_order - min_order
#   
#   ttl <- paste0("Very strong mean-variance relationship in invertebrate abundances in the ", unique(bsh_data$BSH)[1]," broadscale habitat")
#   sbtt <- paste0("Variance within the dataset covers *",orders_of_magnitude_covered," orders of magnitude*.")
#   
#   mtext(side=3, line = 1, at =-0.07, adj=0, cex = 1, ttl, font=1)
#   mtext(side=3, line = 0.25, at =-0.07, adj=0, cex = 0.7, sbtt)
#   
#   dev.off()
#   
#   rm(min_order,max_order,mvpl,min_value,max_value,orders_of_magnitude_covered,ttl,sbtt)
#   toc(log=TRUE)
#   
#   ### run model ####
#   tic(paste0(unique(bsh_data$metadata$PSA)[1], " fit manyglm"))
#   
#   ## untransformed
#   fit.glm <- manyglm(mvabund::as.mvabund(bsh_data$abundance) ~ bsh_data$metadata$year, family = "negative.binomial")
#   fit.glm.summary <- summary(fit.glm)
#   saveRDS(fit.glm, file = paste0("data_out/mvabund.inf.",unique(bsh_data$metadata$PSA)[1],".rdat"))
#   saveRDS(fit.glm.summary,file=paste0("data_out/mvabund.inf.",unique(bsh_data$metadata$PSA)[1],".summary.rdat"))
#   fit.glm.out <- mvabund::anova.manyglm(fit.glm,p.uni = "adjusted", test="LR",show.time="all")
#   saveRDS(fit.glm.out, file = paste0("data_out/mvabund.inf.",unique(bsh_data$metadata$PSA)[1],".pw.rdat"))
#   
#   m2tmp1 <- t(as.data.frame(fit.glm.out$uni.p))[,2]
#   names(m2tmp1) <- names(bsh_data$abundance)
#   
#   print(paste0(length(names(m2tmp1[m2tmp1<0.056]))," 'significant' taxa"))
#   
#   m2tx <- names(m2tmp1[m2tmp1<0.056])#which taxa are 'significantly' different?
#   kptx <- names(bsh_data$abundance) %in% m2tx
#   m2tx <- bsh_data$abundance[, kptx]
#   m2tx$Year <- bsh_data$metadata$year
#   ##make long
#   m2txl <- m2tx %>% 
#     relocate(Year) %>% #move Year to start
#     pivot_longer(cols = c(2:ncol(.)))
#   
#   toc(log=TRUE)
#   
#   tic(paste0(unique(bsh_data$metadata$PSA)[1], " produce ggplots"))
#   
#   ### produce plots ####
#   ggplot(m2txl)+
#     geom_boxplot(aes(x=name,y=value),varwidth = TRUE)+
#     facet_wrap(.~Year)+
#     coord_flip()+
#     labs(y="Taxon abundance",
#          caption=paste0(unique(bsh_data$metadata$PSA)[1]," BSH"))+
#     theme(axis.title.y = element_blank(),
#           # axis.text.x = element_blank(),
#           strip.text = element_text(face="bold")) -> pl2
#   
#   m2txl$Year <- as.factor(m2txl$Year)
#   ggplot(m2txl,aes(x=log(value+1),y=name,
#                    fill=Year,
#                    # colour=Year,stroke=1.5,
#                    shape = Year))+
#     geom_hline(yintercept = seq(from=1.5,
#                                 to=(length(unique(m2txl$name))-.5),
#                                 by=1),
#                colour="lightgrey", linetype=2) +
#     geom_jitter(data=m2txl[,c(2:3)], inherit.aes = FALSE,
#                 aes(x=log(value+1),y=name,),
#                 height = 0.05,size=1, alpha = 0.5, colour = "grey") +
#     geom_jitter(height = 0.05,size=3, alpha = 0.9) +
#     scale_shape_manual(values = c(21:24))+
#     # scale_shape_manual(values = c(1:3))+
#     labs(title = paste0(unique(bsh_data$metadata$PSA)[1]," BSH"),
#          x="log(Taxon abundance (n+1))",
#          caption=paste0("Displayed taxa are the ",paste0(length(names(m2tmp1[m2tmp1<0.056])))," taxa which showed significantly clear differences in abundances between years.\n",
#                         "Taxon abundances across all years are presented in each facet, with abundances for a given year displayed by larger, coloured icons."))+
#     scale_fill_manual(values = cbPalette)+
#     scale_colour_manual(values = cbPalette)+
#     facet_wrap(.~Year)+
#     scale_y_discrete(limits=rev)+
#     theme(
#       legend.position = "none",
#       axis.title.y = element_blank(),
#       axis.text.y = element_text(size=12,face="italic"),
#       axis.title.x = element_text(face="bold"),
#       strip.text.x = element_text(face="bold",size=12),
#       plot.title.position = "plot",
#       plot.title = element_text(face="bold",size=14)
#     ) -> pl3
#   
#   ggsave(filename = paste0("figs/infauna_",unique(bsh_data$metadata$PSA)[1],"_relabund.pdf"),
#          width = 14, height = 6, units="in",plot=pl2)
#   ggsave(filename = paste0("figs/infauna_",unique(bsh_data$metadata$PSA)[1],"_relabund_ver2.pdf"),
#          width = 14, height = 6, units="in",plot=pl3)
#   toc()
#   
#   ### ANOSIM ####
#   tic(paste0(unique(bsh_data$metadata$PSA)[1], " run ANOSIM model"))
#   # untransfomed
#   # fit.anosim <- vegan::anosim(bsh_dataord, group=bsh_data$Year, distance = "bray",permutations = perm)
#   # # log (n+1) transformed
#   # fit.anosim <- vegan::anosim(log(bsh_dataord+1), group=bsh_data$Year, distance = "bray",permutations = perm)
#   # untransformed
#   fit.anosim <- vegan::anosim(bsh_data$abundance,
#                               group=bsh_data$metadata$year,
#                               distance = "bray",
#                               permutations = perm)
#   saveRDS(fit.anosim, file = paste0("data_out/anosim.inf.",
#                                     unique(bsh_data$metadata$PSA)[1],".rdat"))
#   toc(log=TRUE)
#   
#   ### ADONIS2 ####
#   tic(paste0(unique(bsh_data$metadata$PSA)[1], " run PERMANOVA model"))
#   # untransformed
#   fit.adonis2 <- vegan::adonis2(bsh_data$abundance ~ bsh_data$metadata$year,
#                                 permutations = perm)
#   # # transformed
#   # fit.adonis2 <- vegan::adonis2(log(bsh_dataord+1) ~ bsh_data$Year,permutations = perm)
#   saveRDS(fit.adonis2, file = paste0("data_out/adonis2.inf.",
#                                      unique(bsh_data$BSH_CODE)[1],".rdat"))
#   toc(log=TRUE)
#   
#   ### SIMPER ####
#   tic(paste0(unique(bsh_data$metadata$PSA)[1], " run SIMPER"))
#   # untransformed
#   fit.simper <- vegan::simper(bsh_data$abundance, group=bsh_data$metadata$year,
#                               permutations = perm)
#   # # transformed
#   # fit.simper <- vegan::simper(log(bsh_dataord+1), group=bsh_data$Year,
#   #                             permutations = perm)
#   saveRDS(fit.simper, file = paste0("data_out/simper.inf.",
#                                     unique(bsh_data$BSH_CODE)[1],".rdat"))
#   toc(log = TRUE)
#   flush.console()
#   
# }
# 
# toc(log=TRUE)
# toc(log=TRUE)

## ver 2: with offsets ####
tic("ANALYSES by BSH: Run mvabund models by BSH with offset")

### calculate offsets ####
dfw$metadata %>% 
  group_by(year, Site_ID) %>% 
  summarise(reps = n(), .groups = "drop") -> offset
dfw$metadata <- left_join(dfw$metadata, offset, by=c("Site_ID","year"))
rm(offset)

### remove A5.1 ####
remove <- which(dfw$metadata$PSA == "A5.1")
dfw_trim <- lapply(dfw, function(dfw) dfw[-remove,])
rm(remove)

## run MVABUND analysis
for (bshcode in unique(dfw_trim$metadata$PSA)) {
  ### which rows to keep
  # test:
  #bshcode <- "A5.4"
  kp <- which(dfw$metadata$PSA == bshcode)
  # subset data by BSH
  bsh_data <- lapply(dfw, function(dfw) dfw[kp,])
  
  ## remove 'empty' columns
  bsh_data$abundance <- bsh_data$abundance %>%
    dplyr::select_if(where( ~ !is.numeric(.) || sum(.) !=0)) ## remove numeric variables that sum to 0
  
  ### produce meanvar plot ####
  tic(paste0(unique(bsh_data$metadata$PSA)[1], " create meanvar plot"))
  png(file = paste0("figs/infMeanVar.",unique(bsh_data$metadata$PSA)[1],".png"),
      width=12*ppi, height=6*ppi, res=ppi)
  mvpl <- mvabund::meanvar.plot(mvabund(bsh_data$abundance),
                                xlab="Mean",
                                ylab="Variance",
                                table=TRUE)
  
  # Step 1: Find the minimum and maximum values
  min_value <- min(mvpl[,2])
  max_value <- max(mvpl[,2])
  
  min_order <- floor(log10(min_value))
  max_order <- floor(log10(max_value))
  orders_of_magnitude_covered <- max_order - min_order
  
  ttl <- paste0("Very strong mean-variance relationship in invertebrate abundances in the ", unique(bsh_data$BSH)[1]," broadscale habitat")
  sbtt <- paste0("Variance within the dataset covers *",orders_of_magnitude_covered," orders of magnitude*.")
  
  mtext(side=3, line = 1, at =-0.07, adj=0, cex = 1, ttl, font=1)
  mtext(side=3, line = 0.25, at =-0.07, adj=0, cex = 0.7, sbtt)
  
  dev.off()
  
  rm(min_order,max_order,mvpl,min_value,max_value,orders_of_magnitude_covered,ttl,sbtt)
  toc(log=TRUE)
  
  ### run model ####
  tic(paste0(unique(bsh_data$metadata$PSA)[1], " fit manyglm"))
  
  ## untransformed
  fit.glm <- manyglm(mvabund::as.mvabund(bsh_data$abundance) ~ bsh_data$metadata$year,
                     family = "negative.binomial",
                     offset = log(bsh_data$metadata$reps))
  fit.glm.summary <- summary(fit.glm)
  saveRDS(fit.glm, file = paste0("data_out/mvabund.inf.",unique(bsh_data$metadata$PSA)[1],".offset.rdat"))
  saveRDS(fit.glm.summary,file=paste0("data_out/mvabund.inf.",unique(bsh_data$metadata$PSA)[1],".summary.offset.rdat"))
  fit.glm.out <- mvabund::anova.manyglm(fit.glm,p.uni = "adjusted", test="LR",show.time="all")
  saveRDS(fit.glm.out, file = paste0("data_out/mvabund.inf.",unique(bsh_data$metadata$PSA)[1],".pw.offset.rdat"))
  
  m2tmp1 <- t(as.data.frame(fit.glm.out$uni.p))[,2]
  names(m2tmp1) <- names(bsh_data$abundance)
  
  print(paste0(length(names(m2tmp1[m2tmp1<0.1]))," 'significant' taxa"))
  
  m2tx <- names(m2tmp1[m2tmp1<0.1])#which taxa are 'significantly' different?
  kptx <- names(bsh_data$abundance) %in% m2tx
  m2tx <- bsh_data$abundance[, kptx]
  m2tx$Year <- bsh_data$metadata$year
  ##make long
  m2txl <- m2tx %>% 
    relocate(Year) %>% #move Year to start
    pivot_longer(cols = c(2:ncol(.)))
  
  toc(log=TRUE)
  
  tic(paste0(unique(bsh_data$metadata$PSA)[1], " produce ggplots"))
  
  ### produce plots ####
  ggplot(m2txl)+
    geom_boxplot(aes(x=name,y=value),varwidth = TRUE)+
    facet_wrap(.~Year)+
    coord_flip()+
    labs(y="Taxon abundance",
         caption=paste0(unique(bsh_data$metadata$PSA)[1]," BSH"))+
    theme(axis.title.y = element_blank(),
          # axis.text.x = element_blank(),
          strip.text = element_text(face="bold")) -> pl2
  
  m2txl$Year <- as.factor(m2txl$Year)
  ggplot(m2txl,aes(x=log(value+1),y=name,
                   fill=Year,
                   # colour=Year,stroke=1.5,
                   shape = Year))+
    # geom_hline(yintercept = seq(from=1.5,
    #                             to=(length(unique(m2txl$name))-.5),
    #                             by=1),
    #            colour="lightgrey", linetype=2) +
    geom_jitter(data=m2txl[,c(2:3)], inherit.aes = FALSE,
                aes(x=log(value+1),y=name,),
                height = 0.05,size=1, alpha = 0.5, colour = "grey") +
    geom_jitter(height = 0.05,size=3, alpha = 0.9) +
    scale_shape_manual(values = c(21:24))+
    # scale_shape_manual(values = c(1:3))+
    labs(title = paste0(unique(bsh_data$metadata$PSA)[1]," BSH"),
         x="log(Taxon abundance (n+1))",
         caption=paste0("Displayed taxa are the ",paste0(length(names(m2tmp1[m2tmp1<0.1])))," taxa which showed differences in abundances between years at alpha = 0.1.\n",
                        "Taxon abundances across all years are presented in each facet, with abundances for a given year displayed by larger, coloured icons."))+
    scale_fill_manual(values = cbPalette)+
    scale_colour_manual(values = cbPalette)+
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
  
  ggsave(filename = paste0("figs/infauna_",unique(bsh_data$metadata$PSA)[1],"_relabund.pdf"),
         width = 14, height = 6, units="in",plot=pl2)
  ggsave(filename = paste0("figs/infauna_",unique(bsh_data$metadata$PSA)[1],"_relabund_ver2.pdf"),
         width = 14, height = 6, units="in",plot=pl3)
  toc()
  
  ### ANOSIM ####
  tic(paste0(unique(bsh_data$metadata$PSA)[1], " run ANOSIM model"))
  # untransfomed
  # fit.anosim <- vegan::anosim(bsh_dataord, group=bsh_data$Year, distance = "bray",permutations = perm)
  # # log (n+1) transformed
  # fit.anosim <- vegan::anosim(log(bsh_dataord+1), group=bsh_data$Year, distance = "bray",permutations = perm)
  # untransformed
  fit.anosim <- vegan::anosim(bsh_data$abundance,
                              group=bsh_data$metadata$year,
                              distance = "bray",
                              permutations = perm)
  saveRDS(fit.anosim, file = paste0("data_out/anosim.inf.",
                                    unique(bsh_data$metadata$PSA)[1],".rdat"))
  toc(log=TRUE)
  
  ### ADONIS2 ####
  tic(paste0(unique(bsh_data$metadata$PSA)[1], " run PERMANOVA model"))
  # untransformed
  fit.adonis2 <- vegan::adonis2(bsh_data$abundance ~ bsh_data$metadata$year,
                                permutations = perm)
  # # transformed
  # fit.adonis2 <- vegan::adonis2(log(bsh_dataord+1) ~ bsh_data$Year,permutations = perm)
  saveRDS(fit.adonis2, file = paste0("data_out/adonis2.inf.",
                                     unique(bsh_data$BSH_CODE)[1],".rdat"))
  toc(log=TRUE)
  
  ### SIMPER ####
  tic(paste0(unique(bsh_data$metadata$PSA)[1], " run SIMPER"))
  # untransformed
  fit.simper <- vegan::simper(bsh_data$abundance, group=bsh_data$metadata$year,
                              permutations = perm)
  # # transformed
  # fit.simper <- vegan::simper(log(bsh_dataord+1), group=bsh_data$Year,
  #                             permutations = perm)
  saveRDS(fit.simper, file = paste0("data_out/simper.inf.",
                                    unique(bsh_data$BSH_CODE)[1],".rdat"))
  toc(log = TRUE)
  flush.console()
  
}
toc(log = TRUE)
rm(fit.adonis2,fit.anosim,fit.glm,fit.glm.out,fit.glm.summary,fit.simper,
   BSH_codes,bsh_data,m2tx,m2txl,pl2,pl3,bshcode,kp,kptx,m2tmp1,df_long,dfw_trim)

# UNIVARIATE MODELS ####
## calculate taxon richness & Abundance values ####
S <- vegan::specnumber(dfw$abundance)

## convert abundance values of -999 to 0 for faunal density calculation
dfw$abundance_noP <- dfw$abundance_raw
dfw$abundance_noP[dfw$abundance_noP <0] <- 0
N <- rowSums(dfw$abundance_noP)

## create temporary DB
dftmp <- data.frame(Waterbody = dfw$metadata$Waterbody)
dftmp$S <- S;rm(S)
dftmp$N <- N;rm(N)
dftmp$year <- as.factor(dfw$metadata$year)
dftmp$PSA <- as.factor(dfw$metadata$PSA)
dftmp$year <- relevel(dftmp$year, ref="2023")


## TAXON RICHNESS (S) ####
### A5.2 ####
dftmpA5.2 <- dftmp %>% filter(PSA == "A5.2")

summary(fitSA52 <- glm.nb(S ~ year,data = dftmpA5.2))
saveRDS(fitSA52, file = "data_out/glmnb_A52_S.Rdat")
visreg::visreg(fitSA52)

summary(fitSA52.WB <- lme4::glmer.nb(S ~ year + (1|Waterbody),
                                     data=dftmpA5.2))
saveRDS(fitSA52.WB, file="data_out/glmnb.WB_A52_S.Rdat")
capture.output(summary(fitSA52.WB), file="data_out/glmnb.WB_A52_S.txt")
visreg::visreg(fitSA52.WB)
rm(dftmpA5.2,fitSA52,fitSA52.WB)

### A5.3  ####
dftmpA5.3 <- dftmp %>% filter(PSA == "A5.3")

summary(fitSA53 <- glm.nb(S ~ year,data = dftmpA5.3))
saveRDS(fitSA53, file = "data_out/glmnb_A53_S.Rdat")
visreg::visreg(fitSA53)

summary(fitSA53.WB <- lme4::glmer.nb(S ~ year + (1|Waterbody),
                                     data=dftmpA5.3))
saveRDS(fitSA53.WB, file="data_out/glmnb.WB_A53_S.Rdat")
capture.output(summary(fitSA53.WB), file="data_out/glmnb.WB_A53_S.txt")
visreg::visreg(fitSA53.WB)
rm(dftmpA5.3,fitSA53,fitSA53.WB)

### A5.4 ####
dftmpA5.4 <- dftmp %>% filter(PSA == "A5.4")

summary(fitSA54 <- glm.nb(S ~ year,data = dftmpA5.4))
saveRDS(fitSA54, file = "data_out/glmnb_A54_S.Rdat")
visreg::visreg(fitSA54)

summary(fitSA54.WB <- lme4::glmer.nb(S ~ year + (1|Waterbody),
                                     data=dftmpA5.4))
saveRDS(fitSA54.WB, file="data_out/glmnb.WB_A54_S.Rdat")
capture.output(summary(fitSA54.WB), file="data_out/glmnb.WB_A54_S.txt")
visreg::visreg(fitSA54.WB)
rm(dftmpA5.4,fitSA54,fitSA54.WB)

## TAXON Abundance (N) ####
### A5.2 ####
dftmpA5.2 <- dftmp %>% filter(PSA == "A5.2")

summary(fitNA52 <- glm.nb((N) ~ year,data = dftmpA5.2))
saveRDS(fitNA52, file = "data_out/glmnb_A52_N.Rdat")
visreg::visreg(fitNA52)

summary(fitNA52.WB <- lmer(log(N) ~ year + (1|Waterbody),
                           data=dftmpA5.2))
saveRDS(fitNA52.WB, file="data_out/glmnb.WB_A52_N.Rdat")
capture.output(summary(fitNA52.WB), file="data_out/glmnb.WB_A52_N.txt")
visreg::visreg(fitNA52.WB)
rm(dftmpA5.2,fitNA52,fitNA52.WB)

### A5.3  ####
dftmpA5.3 <- dftmp %>% filter(PSA == "A5.3")

summary(fitNA53 <- glm.nb(N ~ year,data = dftmpA5.3))
saveRDS(fitNA53, file = "data_out/glmnb_A53_N.Rdat")
visreg::visreg(fitNA53)

summary(fitNA53.WB <- lmer(log(N) ~ year + (1|Waterbody),
                           data=dftmpA5.3))
saveRDS(fitNA53.WB, file="data_out/glmnb.WB_A53_N.Rdat")
capture.output(summary(fitNA53.WB), file="data_out/glmnb.WB_A53_N.txt")
visreg::visreg(fitNA53.WB)
rm(dftmpA5.3,fitNA53,fitNA53.WB)

### A5.4 ####
dftmpA5.4 <- dftmp %>% filter(PSA == "A5.4")

summary(fitNA54 <- glm.nb(N ~ year,data = dftmpA5.4))
saveRDS(fitNA54, file = "data_out/glmnb_A54_N.Rdat")
visreg::visreg(fitNA54)

summary(fitNA54.WB <- lmer(log(N) ~ year + (1|Waterbody),
                           data=dftmpA5.4))
saveRDS(fitNA54.WB, file="data_out/glmnb.WB_A54_N.Rdat")
capture.output(summary(fitNA54.WB), file="data_out/glmnb.WB_A54_N.txt")
visreg::visreg(fitNA54.WB)
rm(dftmpA5.4,fitNA54,fitNA54.WB)

# MVABUNDS BY WB_BSH ####
tic("MVABUNDS BY WB_BSH")
full_data %>%
  dplyr::select(.,-c(Waterbody_Year,Site,Replicate:Biotope_EUNIS,id)
  ) %>%
  #remove A5.1 and A5.2 due to limited presence across WBs
  filter(.,
         PSA != "A5.1",
         PSA != "A5.2",
         Waterbody != "Blackwater",
         Waterbody != "Colne") -> trim_data

trim_data %>% 
  group_by(year, Site_ID) %>% 
  summarise(reps = n(), .groups = "drop") -> offset
trim_data <- left_join(trim_data, offset, by=c("Site_ID","year")) %>% 
  relocate(reps)
rm(offset)


# Step 3: Get unique levels of Waterbody and PSA
waterbody_levels <- unique(trim_data$Waterbody)
psa_levels <- unique(trim_data$PSA)

# Step 4: Create a list to store results
mvabund_results <- list()

for (waterbody in waterbody_levels) {
  for (psa in psa_levels) {
    # Filter the data for the current combination
    subset_data <- trim_data %>%
      filter(Waterbody == waterbody, PSA == psa)
    
    # Skip if there are not enough rows for ordination
    if (nrow(subset_data) < 2) {
      message(paste("Skipping Waterbody:", waterbody, "and PSA:", psa, "due to insufficient data."))
      next
    }
    
    # Convert to mvabund object for analysis for metaMDS (removing grouping variables)
    abundance_mvab <- mvabund::mvabund(subset_data %>%
                                    dplyr::select(-reps,-Waterbody, -PSA, -Site_ID,-year))
    
    # Run manyglm and store results
    mvabund_result <- manyglm(abundance_mvab ~ subset_data$year,
                              offset = log(subset_data$reps),
                              family = "negative.binomial")
    
    # Save the result in the list
    mvabund_results[[paste(waterbody, psa, sep = "_")]] <- list(
      Waterbody = waterbody,
      PSA = psa,
      mvab = mvabund_result
    )
  }
}

saveRDS(mvabund_results, file="data_out/mvabund_WB_BSH.Rdat")
toc(log=TRUE)

## summarise and export mvabunds ####
tic("summarise and export mvabunds")

# Initialize a list to store summary outputs
summary_list <- list()

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
saveRDS(summary_list, file = "mvabund_summaries.rds")
toc(log=TRUE)


## TO DO:
# produce plots of 'significant' taxa per WB_BSH between years
# update text to summarise analyses