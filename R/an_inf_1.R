# an_inf_1.R ####
# analyses of infaunal data

# Set up ####
## load packages ####
ld_pkgs <- c("tidyverse", "vegan","lmerTest","rstatix", "mvabund","tictoc",
             "MASS","ggtext","ggpmisc", "gllvm")
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

# analyse ALL data together ####
## run ordination ####
tic("run ordination")
ord <- metaMDS(dfw$abundance,try = 200,trymax = 500)
#plot(ord)
toc(log = TRUE)

## plot ####
tic("Extract ordination data for plotting using ggplot")
## extract Site scores
mds_scores <- as_tibble(as.data.frame(scores(ord,"site")))

mds_scores$Biotope <- dfw$metadata$Biotope
mds_scores$Biotope_EUNIS <- dfw$metadata$Biotope_EUNIS
mds_scores$Waterbody <- dfw$metadata$Waterbody
mds_scores$Year <- dfw$metadata$year
mds_scores$Site <- dfw$metadata$Site_ID
mds_scores$PSA <- dfw$metadata$PSA

### extract species scores and groups
spp_scores <- as.data.frame(scores(ord, "species"))  #Using the scores function from vegan to extract the species scores and convert to a data.frame
spp_scores$species <- rownames(spp_scores)  # create a column of species, from the rownames of species.scores
spp_scores$species_sh <- make.cepnames(spp_scores$species)

pdf("figs/inf_mds_all.pdf",
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
    label=PSA,
    colour=as.factor(Year)
  ),
  size=7,
  fontface=2)+
  geom_text_npc(aes(npcx = .99, npcy = .99,
                    label=paste("Stress = ",
                                round(ord$stress, 3))))+
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
rm(ord)
toc(log=TRUE)

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
    label=Site,
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
    label=Site,
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
    label=Site,
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
tic("ANALYSES by BSH: Run mvabund models by BSH")

## remove A5.1 ####
remove <- which(dfw$metadata$PSA == "A5.1")
dfw_trim <- lapply(dfw, function(dfw) dfw[-remove,])
rm(remove)

## run MVABUND analysis
for (bshcode in unique(dfw_trim$metadata$PSA)) {
  ### which rows to keep
  kp <- which(dfw$metadata$PSA == bshcode)
  # subset data by BSH
  bsh_data <- lapply(dfw, function(dfw) dfw[kp,])
  
  ## remove 'empty' columns
  bsh_data$abundance <- bsh_data$abundance %>%
    dplyr::select_if(where( ~ !is.numeric(.) || sum(.) !=0)) ## remove numeric variables that sum to 0
  
  ## produce meanvar plot ####
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
  
  # run model ####
  tic(paste0(unique(bsh_data$metadata$PSA)[1], " fit manyglm"))
  
  ## untransformed
  fit.glm <- manyglm(mvabund::as.mvabund(bsh_data$abundance) ~ bsh_data$metadata$year, family = "negative.binomial")
  fit.glm.summary <- summary(fit.glm)
  saveRDS(fit.glm, file = paste0("data_out/mvabund.inf.",unique(bsh_data$metadata$PSA)[1],".rdat"))
  saveRDS(fit.glm.summary,file=paste0("data_out/mvabund.inf.",unique(bsh_data$metadata$PSA)[1],".summary.rdat"))
  fit.glm.out <- mvabund::anova.manyglm(fit.glm,p.uni = "adjusted", test="LR",show.time="all")
  saveRDS(fit.glm.out, file = paste0("outputs/mvabund.inf.",unique(bsh_data$metadata$PSA)[1],".pw.rdat"))
  
  m2tmp1 <- t(as.data.frame(fit.glm.out$uni.p))[,2]
  names(m2tmp1) <- names(bsh_data$abundance)
  
  print(paste0(length(names(m2tmp1[m2tmp1<0.056]))," 'significant' taxa"))
  
  m2tx <- names(m2tmp1[m2tmp1<0.056])#which taxa are 'significantly' different?
  kptx <- names(bsh_data$abundance) %in% m2tx
  m2tx <- bsh_data$abundance[, kptx]
  m2tx$Year <- bsh_data$metadata$year
  ##make long
  m2txl <- m2tx %>% 
    relocate(Year) %>% #move Year to start
    pivot_longer(cols = c(2:ncol(.)))
  
  toc(log=TRUE)
  
  tic(paste0(unique(bsh_data$metadata$PSA)[1], " produce ggplots"))
  
  # produce plots ####
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
    geom_hline(yintercept = seq(from=1.5,
                                to=(length(unique(m2txl$name))-.5),
                                by=1),
               colour="lightgrey", linetype=2) +
    geom_jitter(data=m2txl[,c(2:3)], inherit.aes = FALSE,
                aes(x=log(value+1),y=name,),
                height = 0.05,size=1, alpha = 0.5, colour = "grey") +
    geom_jitter(height = 0.05,size=3, alpha = 0.9) +
    scale_shape_manual(values = c(21:24))+
    # scale_shape_manual(values = c(1:3))+
    labs(title = paste0(unique(bsh_data$metadata$PSA)[1]," BSH"),
         x="log(Taxon abundance (n+1))",
         caption=paste0("Displayed taxa are the ",paste0(length(names(m2tmp1[m2tmp1<0.056])))," taxa which showed significantly clear differences in abundances between years.\n",
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
  
  # ANOSIM ####
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
  
  # ADONIS2 ####
  tic(paste0(unique(bsh_data$metadata$PSA)[1], " run PERMANOVA model"))
  # untransformed
  fit.adonis2 <- vegan::adonis2(bsh_data$abundance ~ bsh_data$metadata$year,
                                permutations = perm)
  # # transformed
  # fit.adonis2 <- vegan::adonis2(log(bsh_dataord+1) ~ bsh_data$Year,permutations = perm)
  saveRDS(fit.adonis2, file = paste0("outputs/adonis2.inf.",
                                     unique(bsh_data$BSH_CODE)[1],".rdat"))
  toc(log=TRUE)
  
  # SIMPER ####
  tic(paste0(unique(bsh_data$metadata$PSA)[1], " run SIMPER"))
  untransformed
  fit.simper <- vegan::simper(bsh_data$abundance, group=bsh_data$metadata$year,
                              permutations = perm)
  # # transformed
  # fit.simper <- vegan::simper(log(bsh_dataord+1), group=bsh_data$Year,
  #                             permutations = perm)
  saveRDS(fit.simper, file = paste0("outputs/simper.inf.",
                                    unique(bsh_data$BSH_CODE)[1],".rdat"))
  toc(log = TRUE)
  flush.console()
  
}

toc(log=TRUE)
toc(log=TRUE)

###########################################################################
  ###########################################################################
  ###########################################################################
  ###########################################################################
  ###########################################################################
  ###########################################################################
  ###########################################################################
  ###########################################################################
  ###########################################################################
  ###########################################################################
  ################ FROM HERE ################################################
  ###########################################################################
  ###########################################################################
  ###########################################################################
  ###########################################################################
  ###########################################################################
  ###########################################################################
  ###########################################################################
  ###########################################################################
  ###########################################################################
  ###########################################################################
  

