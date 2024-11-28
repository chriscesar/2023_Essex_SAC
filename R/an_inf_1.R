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
plot(ord)
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

# split by BSH and reanalyse ####
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
##
unlist(tictoc::tic.log())
