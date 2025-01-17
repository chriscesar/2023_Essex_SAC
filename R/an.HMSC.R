# an.HMSC.R ###
# Testing use of HMSC method to deal with complex data structure

# Set up ####
## load packages ####
ld_pkgs <- c("tidyverse", "tictoc","Hmsc","ggpubr")
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

# prep data ####
## Y: Species occurrence ####
Y <- dfw$abundance
Y[Y>0] <- 1 #convert to presence-absence

## consider reducing number of taxa
# e.g. remove taxa that occur in n or fewer samples
n <- 5
Y <- Y %>%
  select(where(~ sum(.) >= n)) %>% 
  as.matrix(.)

ns <- ncol(Y)
###################
## Variable structuring ####

studyDesign0 <- dfw$metadata[,c(2,8,9,14,15)]
Waterbody <- as.factor(studyDesign$Waterbody)
PSA <- as.factor(studyDesign$PSA)
studyDesign <- data.frame(Waterbody,PSA)

## X: Environmental ####
Year <- as.factor(dfw$metadata$year)

X <- data.frame(Year)
X$Year <- as.factor(X$Year)
XFormula <- ~Year

## create model ####
# create model
simul <- Hmsc(Y=Y, XData = X,
              XFormula=XFormula,
              studyDesign = studyDesign,
              ranLevels = ranlevels,
              distr = "probit")

tic("Run model")
# Run model

test.run = TRUE ###'TRUE' runs 'shorter' version of full ('FALSE') model

if (test.run){
  thin = 10
  nChains = 1
  samples = 100
  transient = ceiling(thin*samples*.5)
}else{
  thin = 10
  nChains = 2
  samples = 1000
  transient = ceiling(thin*samples*.5)
}

## run model ####
mod_HMSC <- sampleMcmc(simul,
                       samples = samples,
                       thin = thin,
                       transient = transient,
                       nChains = nChains#,
                       # nParallel = nChains
                       )
saveRDS(mod_HMSC, file = paste0("data_out/mod_HMSC","_smp",samples,
                                "_thn",thin,"_trns",transient,"_chn",nChains,
                                ".Rdata"))
mod_HMSC <- readRDS(file = paste0("data_out/mod_HMSC","_smp",samples,
                                  "_thn",thin,"_trns",transient,"_chn",nChains,
                                  ".Rdata"))
print(paste0("MODEL RUN COMPLETE. NOTE THAT MODEL IS BASED ON PRESENCE-ABSENCE DATA OF SPECIES PRESENT IN ",n, "OR MORE SAMPLES"))
toc(log = TRUE)

# model diagnostics ####
## investigate model outputs ####
mpost <- convertToCodaObject(mod_HMSC) # model diagnostics/convergence
preds <- computePredictedValues(mod_HMSC) # model performance
MF <- evaluateModelFit(hM=mod_HMSC, predY = preds) # r2, etc
VP <- computeVariancePartitioning(mod_HMSC) # variance partitioning

# ### to check ####
# #VP warning:
# #In cor(lbeta[[i]][k, ], lmu[[i]][k, ]) : the standard deviation is zero
# 
ess.beta <- effectiveSize(mpost$Beta) %>%
  as_tibble() %>% dplyr::rename(ess_beta=value)
psrf.beta <- gelman.diag(mpost$Beta, multivariate = FALSE)$psrf %>%
  as_tibble() %>% dplyr::rename(psrf_beta = "Point est.")

diag_all <- ggarrange(ggplot(ess.beta,aes(x=ess_beta))+
                        geom_histogram()+
                        xlab("Effective sample size"),
                      ggplot(psrf.beta,aes(x=psrf_beta))+
                        geom_histogram()+
                        geom_vline(xintercept = 1.001, col=2)+
                        xlab("Gelman diagnostic"), align = "v")+
  ggtitle("All plots");diag_all

hist(ess.gamma <- effectiveSize(mpost$Gamma))
hist(psrf.gamma <- gelman.diag(mpost$Gamma, multivariate=FALSE)$psrf)

sppairs = matrix(sample(x = 1:ns^2, size = 100))
tmp = mpost$Omega[[1]]
for (chain in 1:length(tmp)){
  tmp[[chain]] = tmp[[chain]][,sppairs]
}
ess.omega = effectiveSize(tmp)
psrf.omega = gelman.diag(tmp, multivariate=FALSE)$psrf
hist(ess.omega)
hist(psrf.omega)

preds = computePredictedValues(simul)
MF = evaluateModelFit(hM=m, predY=preds)
hist(MF$R2, xlim = c(0,1), main=paste0("Mean = ", round(mean(MF$R2),2)))

MF$TjurR2 %>% mean(na.rm=TRUE)*100

# species niches
postBeta <- getPostEstimate(mod_HMSC, parName = "Beta")

plotVariancePartitioning(mod_HMSC, VP=VP, las=2, horiz=TRUE)
plotBeta(mod_HMSC,post=postBeta, param = "Support", #supportLevel=0.95,
         split = .4, spNamesNumbers = c(T,F))
