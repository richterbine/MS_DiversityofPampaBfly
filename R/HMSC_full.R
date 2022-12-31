# full script to make my thesis happen :)


# STEP 1 - model structure ------------------------------------------------

library(Hmsc)
library(ape)
library(ggplot2)

# read the data
da <- readRDS(here::here("data/processed/All_Data.rds"))
head(da)
str(da)

# Create the Plot variable to split out the data
da$Plot <- as.factor(paste(da$Local, da$SO, sep = ""))
da$Habitat <- as.factor(da$Habitat)
str(da)

colnames(da)[which(colnames(da) == "Forest Formation")] <- "Forest"

# Species mean traits data (T)
TrData <- as.data.frame(readRDS(here::here("data/processed/Mean_traits.rds")))
head(TrData)
rownames(TrData) <- TrData$Species

# phylogenetic relationship among 35 species
phyloTree <- readRDS(here::here("data/processed/tree_bfly_35spp.rds"))
phyloTree$tip.label <- gsub("_", " ", phyloTree$tip.label)
plot(phyloTree, cex = 0.5)


# defining the fixed effect part
# environmental covariates
XFormula <- ~ Temp + Forest + Grasslands + Crops

# As trait formula, we include mean forewing length and the mean wing-thorax ratio
TrFormula <- ~ mWTR

# create an empty list to store the data to be split
comm.sub <- Y <- Y.pa <- Y.hell <- xData <- xy <- Tr.sub <- phy.sub <- 
  studyDesign <- rL <- list()

## Split the data
for (i in 1:length(unique(da$Plot))) {
  comm.sub[[i]] <- subset(da, Plot == unique(da$Plot)[i])
}
names(comm.sub) <- unique(da$Plot)
comm.sub <- comm.sub[-1] # we will remove the data from EE1 (piloto)

for (i in 1:length(comm.sub)) {
  Y[[i]] <- as.matrix(comm.sub[[i]][,10:44])
  Y[[i]] <- Y[[i]][,-which(colSums(Y[[i]]) < 1)]
  Y.pa[[i]] <- vegan::decostand(Y[[i]], method = "pa")
  Y.hell[[i]] <- vegan::decostand(Y[[i]], method = "helli")
  
  # Environmental covariates (X)
  xData[[i]] <- comm.sub[[i]][,c(5, 8, 9, 45:47)]
  
  # spatial coordinates
  xy[[i]] <- as.matrix(comm.sub[[i]][, 6:7])
  rownames(xy[[i]]) <- comm.sub[[i]]$SiteID
  
  # prepare the traits data
  Tr.sub[[i]] <- TrData[match(colnames(Y[[i]]), rownames(TrData)),] # match
  Tr.sub[[i]]
  
  match(colnames(Y[[i]]), phyloTree$tip.label)
  
  phy.sub[[i]] <- geiger::treedata(phyloTree, setNames(colnames(Y[[i]]), 
                                                       colnames(Y[[i]])))$phy
  phy.sub[[i]]
  
  # We are now ready to define the HMSC model. 
  # define the random structure of the data
  studyDesign[[i]] <- data.frame(SU = as.factor(comm.sub[[i]]$SiteID))
  rL[[i]] <- HmscRandomLevel(sData = xy[[i]])
  
}

for (i in 1:length(Y)) {
  print(cbind(S = (rowSums(Y.pa[[i]])), N = (rowSums(Y[[i]], na.rm = T))))
  print(sum(Y[[i]]))
  print(dim(Y[[i]])[2])
}

# HMSC model
m.pa <- m.hell <- out.pa <- out.hell <- list()

nChains = 2  # number of independent MCMC chains to be run
thin = 10   # the number of MCMC steps between each recording of samples from the posterior
samples = 1000   # the number of MCMC samples to be obtained in each chain
transient <- 100*thin   # the number of MCMC steps that are executed before starting recording posterior samples
nParallel <- 2 # optional setting of nParallel

for (i in 1:length(xData)) {
  m.pa[[i]] <- Hmsc(Y = Y.pa[[i]], XData = xData[[i]], XFormula = XFormula, 
                    phyloTree = phy.sub[[i]], TrData = Tr.sub[[i]], TrFormula = TrFormula,
                    distr = "probit", studyDesign = studyDesign[[i]], 
                    ranLevels = list(SU = rL[[i]])) # Presence-absence data
  
  m.hell[[i]] <- Hmsc(Y = Y.hell[[i]], XData = xData[[i]], XFormula = XFormula, 
                      phyloTree = phy.sub[[i]], TrData = Tr.sub[[i]], TrFormula = TrFormula,
                      distr = "normal", studyDesign = studyDesign[[i]], 
                      ranLevels = list(SU = rL[[i]]), YScale = TRUE) # Normal

  out.pa[[i]] <- sampleMcmc(m.pa[[i]], thin = thin, samples = samples, 
                            transient = transient, nChains = nChains,
                            nParallel = nParallel)
  
  out.hell[[i]] <- sampleMcmc(m.hell[[i]], thin = thin, samples = samples, 
                              transient = transient, nChains = nChains,
                              nParallel = nParallel)
  }

saveRDS(out.pa, here::here("output/out_pa.2_1000samp.rds"))
saveRDS(out.hell, here::here("output/out_hell.2_1000samp.rds"))


# STEP 2 - Checking MCMC convergence -------------------------------------

library(Hmsc)
library(ggplot2)

# raw data
da <- readRDS(here::here("data/processed/All_Data.rds"))
da$Plot <- as.factor(paste(da$Local, da$SO, sep = ""))
str(da)

# read the models outputs
all.models <- list(Probit = NA, NormHell = NA)

all.models[[1]] <- readRDS(here::here("output/out_pa.2_1000samp.rds")) # poisson
all.models[[2]] <- readRDS(here::here("output/out_hell.2_1000samp.rds")) # normal

# we gonna take a loop for all model structure for examination of the MCMC convergence

for (m in 1:length(all.models)) {
  #m=1
  models <- all.models[[m]]
  
  # lets create a loop to examinate MCMC convergence for all local communities
  # for parameters related to fixed effects (beta, gamma and rho) and random effects 
  # (alpha - spatial) and omega (residual association matrix)
  
  mpost <- psrf.beta <- psrf.gamma <- psrf.rho <- psrf.omega <- psrf.alpha <- list()
  
  for (i in 1:length(models)) {
    mpost[[i]] <- convertToCodaObject(models[[i]], spNamesNumbers = c(T, F),
                                      covNamesNumbers = c(T, F))
    
    psrf.beta[[i]] <- gelman.diag(mpost[[i]]$Beta, multivariate = FALSE)$psrf
    
    psrf.gamma[[i]] <- gelman.diag(mpost[[i]]$Gamma, multivariate = FALSE)$psrf
    
    psrf.rho[[i]] <- gelman.diag(mpost[[i]]$Rho, multivariate = FALSE, autoburnin = FALSE)$psrf
    
    psrf.alpha[[i]] <- gelman.diag(mpost[[i]]$Alpha[[1]], multivariate = FALSE)$psrf
    
    psrf.omega[[i]] <- gelman.diag(mpost[[i]]$Omega[[1]], multivariate =  FALSE)$psrf
    
  }
  
  ls.beta <- ls.gamma <- ls.rho <- ls.alpha <- ls.omega <- ls.names <- list()
  
  for (l in 1:length(psrf.beta)) {
    ls.beta[[l]] <- as.vector(psrf.beta[[l]])
    ls.gamma[[l]] <- as.vector(psrf.gamma[[l]])
    ls.rho[[l]] <- as.vector(psrf.rho[[l]])
    ls.alpha[[l]] <- as.vector(psrf.alpha[[l]])
    ls.omega[[l]] <- as.vector(psrf.omega[[l]])
  }
  
  data.all <- data.frame(values = c(unlist(ls.beta), unlist(ls.gamma),
                                    unlist(ls.rho), unlist(ls.alpha), 
                                    unlist(ls.omega)),
                         Val.tipe = NA,
                         Comm = NA,
                         ID.param = NA)
  
  tmp <- list(list(), list(), list(), list(), list())
  for (l in 1:length(psrf.alpha)) {
    tmp[[1]][[l]] <- rep(c("PointEst", "UppCI"), each = length(ls.beta[[l]])/2)
    tmp[[2]][[l]] <- rep(c("PointEst", "UppCI"), each = length(ls.gamma[[l]])/2)
    tmp[[3]][[l]] <- rep(c("PointEst", "UppCI"), each = length(ls.rho[[l]])/2)
    tmp[[4]][[l]] <- rep(c("PointEst", "UppCI"), each = length(ls.alpha[[l]])/2)
    tmp[[5]][[l]] <- rep(c("PointEst", "UppCI"), each = length(ls.omega[[l]])/2)
  }
  data.all$Val.tipe <- unlist(tmp)
  
  # attaching the names to parameters values
  tmp <- list(list(), list(), list(), list(), list())
  for (l in 1:length(psrf.alpha)) {
    # l=1
    tmp[[1]][[l]] <- rep(levels(da$Plot)[-1][l], length(unlist(ls.beta[[l]])))
    tmp[[2]][[l]] <- rep(levels(da$Plot)[-1][l],length(unlist(ls.gamma[[l]])))
    tmp[[3]][[l]] <- rep(levels(da$Plot)[-1][l],length(unlist(ls.rho[[l]])))
    tmp[[4]][[l]] <- rep(levels(da$Plot)[-1][l],length(unlist(ls.alpha[[l]])))
    tmp[[5]][[l]] <- rep(levels(da$Plot)[-1][l],length(unlist(ls.omega[[l]])))
  }
  unlist(tmp)
  
  data.all$Comm <- unlist(tmp)
  
  data.all$ID.param <- c(rep("beta", length(unlist(psrf.beta))),
                         rep("gamma", length(unlist(psrf.gamma))),
                         rep("rho", length(unlist(psrf.rho))),
                         rep("alpha", length(unlist(psrf.alpha))),
                         rep("Omega", length(unlist(psrf.omega))))
  
  
  # density plot for the scale reduction factor values of important parameters estimated by
  # HMSC models
  data.all1 <- data.all
  data.all <- data.all[-which(data.all$ID.param == "rho"),]
  all.MCMC <- ggplot(data.all, 
                     aes(x = values, colour = Comm, fill = Comm)) +
    geom_density(alpha = .3) + facet_wrap(~ ID.param, scales = "free",
                                          labeller = label_parsed, ncol = 4) +
    geom_vline(aes(xintercept = 1.1), 
               linetype = "dashed") + #theme(legend.position = "none") +
    scale_color_viridis_d(option = "A", name = "Local \nComm") +
    scale_fill_viridis_d(option = "A", name = "Local \nComm") 
  all.MCMC
  
  if(m == 1){
    # verifying the psrf for phylogenetic signal rho
    write.table(subset(data.all1, ID.param == "rho"), here::here("output/Final_tables/rho_pa.txt"))
    mcmc.pa <- all.MCMC  
    # cowplot::save_plot(here::here("output/figures/psrf_Gelman_Probit.pdf"),
    #                    all.MCMC + labs(x = "Potential scale reduction factors", tag = "a)"), 
    #                    base_width = 12, base_height = 3)
  }
  if(m == 2){
    write.table(subset(data.all1, ID.param == "rho"), here::here("output/Final_tables/rho_normal.txt"))
    mcmc.hell <- all.MCMC
    # cowplot::save_plot(here::here("output/figures/psrf_Gelman_Normal.pdf"),
    #                    all.MCMC + labs(x = "Potential scale reduction factors", tag = "c)"), 
    #                    base_width = 12, base_height = 3)
  }
}

leg <- cowplot::get_legend(mcmc.pa)
p.mcmc <- cowplot::plot_grid(mcmc.pa + labs(x = NULL, tag = "a)") + theme(legend.position = "none"),
                   mcmc.hell + labs(x = "Potential scale reduction factors", tag = "b)") +
                     theme(legend.position = "none"),
                   nrow = 2)
p.mcmc

cowplot::plot_grid(p.mcmc, leg, ncol = 2, rel_widths = c(2,.5))
cowplot::save_plot(here::here("output/Final_figures/FigS2_psrfModels.png"),
                   cowplot::plot_grid(p.mcmc, leg, ncol = 2, rel_widths = c(2,.5)),
                   base_width = 10, base_height = 6)



  # The MCMC convergence diagnostics can be considered satisfactory,
  # as for most parameters the potential scale reduction factor is 
  # close to the ideal value of one.
  
  # for probit distributions the convergence can be considered satisfactory
  # for alpha and omega parameters.


# STEP 3 - Evaluating the model fit ---------------------------------------

# Extracting the explanatory power of the models 

library(Hmsc)
library(hablar)

# the explanatory power was accessed by R2 for normal model and Tjur R2 for non-normal models

out.pa <- readRDS(here::here("output/out_pa.2_1000samp.rds"))
out.hell <- readRDS(here::here("output/out_hell.2_1000samp.rds"))

nChains = 2  # number of independent MCMC chains to be run
thin = 10  # the number of MCMC steps between each recording of samples from the posterior
samples = 1000   # the number of MCMC samples to be obtained in each chain
transient = 100*thin   # the number of MCMC steps that are executed before starting recording posterior samples
nParallel <- 2 # optional setting of nParallel

result.pa <- result.hell <- matrix(NA, nrow = length(out.pa), ncol = 8)

set.seed(2022)

# for the poisson model
for (i in 1:length(out.pa)) {
  #i=1
  preds <- computePredictedValues(out.pa[[i]], expected = F)
  Mean.richness <- mean(apply(preds, 3, sum))
  MF <- evaluateModelFit(hM = out.pa[[i]], predY = preds)
  semNan <- s(MF$TjurR2)
  mean <- mean(semNan, na.rm = T)
  deviation <- sd(semNan, na.rm = T)
  quartile <- quantile(semNan, na.rm = T)
  result.pa[i,] <- round(c(Mean.richness, mean, deviation, quartile), 3)
  
  # for the normal model
  preds <- computePredictedValues(out.hell[[i]], expected = F)
  Mean.relAbund <- mean(apply(preds, 3, sum))
  R2.hmsc <- evaluateModelFit(hM = out.hell[[i]], predY = preds)
  semNan <- s(R2.hmsc$R2)
  mean <- mean(semNan, na.rm = T)
  deviation <- sd(semNan, na.rm = T)
  quartile <- quantile(semNan, na.rm = T)
  result.hell[i,] <- round(c(Mean.relAbund, mean, deviation, quartile), 3)
}

result.pa <- cbind(result.pa, unlist(lapply(Y.pa, sum)))

result.hell <- cbind(result.hell, unlist(lapply(Y.hell, sum)))

colnames(result.pa) <-c("Mean Est", "EP", "deviation", "qua_0", 
                        "qua_25", "qua_50", "qua_75", "qua_100", "Mean Obs")

rownames(result.pa) <- rownames(result.hell) <- names(comm.sub)

df.MF <- data.frame(Model = rep(c("PA", "Normal"), each = 7),
                    rbind(result.pa, result.hell))
df.MF

write.table(df.MF, here::here("output/Final_tables/ModelFit_all.txt"))


# STEP 4 - Exploring the parameters estimated from beta, alpha, gamma, and omega --------

library(Hmsc)
library(dplyr)
library(ggplot2)
library(cowplot)

# if (!require("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# BiocManager::install("ggtree")
# BiocManager::install("ggtreeExtra")

library(ggtree)
library(ggtreeExtra)
library(ggstance)


# read the models outputs
all.models <- list(Probit = NA, Normal = NA)

all.models[[1]] <- readRDS(here::here("output/out_pa.2_1000samp.rds")) # poisson
all.models[[2]] <- readRDS(here::here("output/out_hell.2_1000samp.rds")) # normal

# Run lines #349 to #371 to obtain the VP and beta parameters
# result <- list()
# 
# for (m in 1:length(all.models)) {
#   #m=2
#   models <- all.models[[m]]
#   
#   VP <- postBeta <- postGamma <- list()
#   
#   for (i in 1:length(models)) {
#     # variance partitioning
#     VP[[i]] <- computeVariancePartitioning(models[[i]])
#     
#     # estimated betas and gammas for each community
#     postBeta[[i]] <- getPostEstimate(models[[i]], parName = "Beta")
#     postGamma[[i]] <- getPostEstimate(models[[i]], parName = "Gamma")
#     
#   }
#   result[[m]] <- c(VP = VP, postBeta = postBeta, postGamma = postGamma)
#   
# }
# names(result) <- c("PA", "Normal")
# saveRDS(result, here::here("output/Final_tables/VP_and_fixParam.rds"))


result <- readRDS(here::here("output/Final_tables/VP_and_fixParam.rds"))
result[[1]]$postBeta1$mean[1:3,1:3]
result[[2]]$postBeta1$mean[1:3,1:3]
result$PA$VP1

# read the species data
bfly.species <- read.csv(here::here("data/raw/Species_data.csv"),
                         header = TRUE, sep = ";")
str(bfly.species)
bfly.species$Colors <- ifelse(bfly.species$Subfamily == "Biblidinae", "#0072B2",
                              ifelse(bfly.species$Subfamily == "Charaxinae", "#44AA99",
                                     ifelse(bfly.species$Subfamily == "Nymphalinae", "#CC79A7", "#E69F00")))
levels(as.factor(bfly.species$Colors))

tree <- list()
p.beta <- list(PA = list(), Normal = list())

for (m in 1:length(result)) {
  for (i in 8:14) { # only betas
    #m=1;i=8
    post <- result[[m]][[i]]
    
    betaP = post$support
    supportLevel = .9
    toPlot = 2 * betaP - 1
    toPlot = toPlot * ((betaP > supportLevel) + (betaP < 
                                                   (1 - supportLevel)) > 0)
    
    df.betas <- as.data.frame(t(toPlot))
    colnames(df.betas) <- all.models[[m]][[1]]$covNames
    
    phy <- all.models[[m]][[i-7]]$phyloTree
    
    df.betas <- df.betas[match(phy$tip.label, rownames(df.betas)),]
    #df.betas$label <- rownames(df.betas)
    
    df.betas$label <- paste(matrix(stringr::str_sub(unlist(strsplit(rownames(df.betas), " ")),1,3), byrow = T,
                                   ncol = 2)[,1], matrix(stringr::str_sub(unlist(strsplit(rownames(df.betas), " ")),1,3), byrow = T,
                                                         ncol = 2)[,2], sep = ".")
    print(paste("No doubles =",length(unique(df.betas$label)) == length(df.betas$label)))
    
    if(sum(match(df.betas$label, "Arc.dem"), na.rm = T) > 1){
      df.betas$label[which(df.betas$label == "Arc.dem")[1]] <- "Arc.dem1"
    }
    
    df.betas$Subfamily <- bfly.species[match(rownames(df.betas), bfly.species$Species), "Subfamily"]
    df.betas$colors <- bfly.species[match(rownames(df.betas), bfly.species$Species), "Colors"]
    
    t <- paste(matrix(stringr::str_sub(unlist(strsplit(phy$tip.label, " ")),1,3), byrow = T,
                      ncol = 2)[,1], matrix(stringr::str_sub(unlist(strsplit(phy$tip.label, " ")),1,3), byrow = T,
                                            ncol = 2)[,2], sep = ".")
    
    if(sum(match(t, "Arc.dem"), na.rm = T) > 1){
      t[which(t == "Arc.dem")[1]] <- "Arc.dem1"
    }
    
    phy$tip.label <- t
    
    tree[[i-7]] <- dplyr::full_join(phy, df.betas, by = "label")
    
    p1 <- ggtree(tree[[i-7]]) +
      geom_tiplab(size = 2, offset = 60, fontface = 3) +
      xlim(0,200)  + 
      geom_treescale(width = 10, x = 0, color = "black") +
      geom_tippoint(aes(color = Subfamily)) +
      scale_color_manual(values = levels(as.factor(df.betas$colors))) 
    
    rownames(df.betas) <- df.betas$label
    
    p.beta[[m]][[i-7]] <- gheatmap(p1, df.betas[,-match(c("label", "Subfamily", "colors"), 
                                                        colnames(df.betas))], offset = .5, 
                                   width = .6, font.size = 2, 
                                   custom_column_labels = c("Int", "Temp", "Forest", "Grass", "Crop"), 
                                   colnames_position = "top", 
                                   colnames_angle = 45, colnames_offset_y = 0, 
                                   hjust = 0, colnames_offset_x = 3) +
      scale_fill_gradient2() + theme(legend.position = "none") +
      coord_cartesian(clip = 'off')
  }
}


leg <- cowplot::get_legend(p.beta[[1]][[1]] + theme(legend.position = "right"))
p.beta$PA[[8]] <- leg
p.beta$Normal[[8]] <- leg

cowplot::plot_grid(plotlist = p.beta$PA, nrow = 2, labels = c("EEA2", "JG1", "JG2", "QR1", "QR2",
                                                              "SG1", "SG2"), label_size = 10,
                   label_fontface = "plain")

cowplot::save_plot(here::here("output/Final_figures/FigS3_PA_betas.png"), 
                   cowplot::plot_grid(plotlist = p.beta$PA, ncol = 4, 
                                      labels = c("EEA2", "JG1", "JG2", "QR1", "QR2", "SG1", "SG2"), 
                                      label_size = 10, label_fontface = "plain"),
                   base_height = 6, base_width = 12)


cowplot::plot_grid(plotlist = p.beta$Normal, nrow = 2, 
                   labels = c("EEA2", "JG1", "JG2", "QR1", "QR2", "SG1", "SG2"), 
                   label_size = 10, label_fontface = "plain")

cowplot::save_plot(here::here("output/Final_figures/FigS4_Normal_betas.png"), 
                   cowplot::plot_grid(plotlist = p.beta$Normal, ncol = 4, 
                                      labels = c("EEA2", "JG1", "JG2", "QR1", "QR2", "SG1", "SG2"), 
                                      label_size = 10, label_fontface = "plain"),
                   base_height = 6, base_width = 12)


# GAMMA - Traits and niche

post.gammas <- p.gamma <- list(PA = vector("list", 7), Normal = vector("list", 7))
names(post.gammas$PA) <- names(post.gammas$Normal) <- c("EEA1", "JG1", "JG2", "QR1", "QR2", "SG1", "SG2")
supportLevel = .9

for (m in 1:length(result)) {
  for (i in 1:7) { # only gammas
    #m=2; i=1
    post <- result[[m]][[i+14]]
    gammaP = post$support
    toPlot = 2 * gammaP - 1
    toPlot = toPlot * ((gammaP > supportLevel) + (gammaP < 
                                                   (1 - supportLevel)) > 0)
    colnames(toPlot) <- all.models[[1]][[1]]$trNames
    rownames(toPlot) <- all.models[[m]][[1]]$covNames
    
    post.gammas[[m]][[i]] <- reshape2::melt(toPlot)
    
    p.gamma[[m]][[i]] <- ggplot(data = post.gammas[[m]][[i]], aes(x = Var1, y = Var2, fill = value)) + 
      geom_tile(color = "white", size = 3, show.legend = F) + labs(x = NULL, y = NULL, fill = NULL,
                                                  title = paste(names(post.gammas[[m]])[i], " (", supportLevel,")", sep = "")) +
      scale_fill_gradient2(limit = c(-1, 1)) +
      theme_classic() #+  coord_fixed(ratio = 1)
  }
}

get.leg <- ggplot(data = post.gammas[[m]][[i]], aes(x = Var1, y = Var2, fill = value)) + 
  geom_tile(color = "white", size = 3) + labs(x = NULL, y = NULL, fill = NULL,
                                                               title = paste(names(post.gammas[[m]])[i], " (", supportLevel,")", sep = "")) +
  scale_fill_gradient2(limit = c(-1, 1)) +
  theme_classic()

legend <- get_legend(
  # create some space to the left of the legend
  get.leg + theme(legend.box.margin = margin(0, 0, 0, 12))
)

cowplot::plot_grid(plotlist = p.gamma$PA)
# none association at a support level of 0.9

# save_plot(here::here("output/figures/PA_gamma.png"), 
#           cowplot::plot_grid(p.gamma$PA[[4]], p.gamma$PA[[6]],
#                              legend, ncol = 3, rel_widths = c(2,2,1)),
#           base_width = 10, base_height = 4)


cowplot::plot_grid(plotlist = p.gamma$Normal)

save_plot(here::here("output/Final_figures/FigS4_Normal_gamma.png"), 
          cowplot::plot_grid(p.gamma$Normal[[3]], p.gamma$Normal[[4]], 
                             p.gamma$Normal[[5]], p.gamma$Normal[[6]],
                             p.gamma$Normal[[7]],
                             legend, ncol = 3, 
                             rel_widths = c(2,2,2,2,2,1)),
          base_width = 10, base_height = 6)


# RHO - phylogenetic signal and alpha - spatial signal
mpost <- ls.alpha <- list(PA = vector("list", 7), Normal = vector("list", 7))

ls.rho <- list(PA = matrix(NA, ncol = 3, nrow = 7), Normal = matrix(NA, ncol = 3, nrow = 7))

for (m in 1:length(all.models)) {
  for (i in 1:length(all.models[[1]])) {
    mpost[[m]][[i]] <- convertToCodaObject(all.models[[m]][[i]])
    ls.rho[[m]][i,] <- round(summary(mpost[[m]][[i]]$Rho, quantiles = c(0.025, 0.5, 0.975))[[2]],2)
    ls.alpha[[m]][[i]] <- round(summary(mpost[[m]][[i]]$Alpha[[1]], quantiles = c(0.025, 0.5, 0.975))[[2]],2)
  }
}

ls.rho
ls.alpha

df.rho.alpha <- cbind(PA = rbind(ls.alpha$PA[[1]][1:2,], ls.rho$PA[1,],
                            ls.alpha$PA[[2]][1:2,], ls.rho$PA[2,],
                            ls.alpha$PA[[3]][1:2,], ls.rho$PA[3,],
                            ls.alpha$PA[[4]][1:2,], ls.rho$PA[4,],
                            ls.alpha$PA[[5]][1:2,], ls.rho$PA[5,],
                            ls.alpha$PA[[6]][1:2,], ls.rho$PA[6,],
                            ls.alpha$PA[[7]][1:2,], ls.rho$PA[7,]), 
                      Normal = rbind(ls.alpha$Normal[[1]][1:2,], ls.rho$Normal[1,],
                            ls.alpha$Normal[[2]][1:2,], ls.rho$Normal[2,],
                            ls.alpha$Normal[[3]][1:2,], ls.rho$Normal[3,],
                            ls.alpha$Normal[[4]][1:2,], ls.rho$Normal[4,],
                            ls.alpha$Normal[[5]][1:2,], ls.rho$Normal[5,],
                            ls.alpha$Normal[[6]][1:2,], ls.rho$Normal[6,],
                            ls.alpha$Normal[[7]][1:2,], ls.rho$Normal[7,]))
write.table(df.rho.alpha, here::here("output/Final_tables/Rho_Alpha.txt"))

# OMEGA - residual correlation

library(corrplot)

p.omega <- list(PA = vector("list", 7), Normal = vector("list", 7))
comm.names <- c("EEA2", "JG1", "JG2", "QR1", "QR2", "SG1", "SG2")

for (m in 1:length(all.models)) {
  for (i in 1:length(all.models[[m]])) {
    #m=1;i=1
    OmegaCor <- computeAssociations(all.models[[m]][[i]])
    supportLevel <- 0.9
    toPlot <- ((OmegaCor[[1]]$support > supportLevel)
               + (OmegaCor[[1]]$support < (1 - supportLevel)) > 0) * OmegaCor[[1]]$mean
    toPlot <- sign(toPlot)
    plotOrder <- match(all.models[[m]][[i]]$phyloTree$tip.label, rownames(toPlot))
    toPlot <- toPlot[plotOrder,plotOrder]
    # corrplot(toPlot, method = "color", tl.cex = 0.4,
    #          col = colorRampPalette(c("blue", "white", "red"))(255))
    
    df.omega <- reshape2::melt(toPlot)
    
    p.omega[[m]][[i]] <- ggplot(data = df.omega, aes(x = Var1, y = Var2, fill = value)) + 
      geom_tile(show.legend = F, color = "gray90", size = .8) +
      scale_fill_gradient2(limit = c(-1, 1)) + labs(x = NULL, y = NULL) +
      labs(title = paste(comm.names[i])) +
      theme(axis.text.x = element_blank(), panel.background = element_blank(),
            title = element_text(size = 10))
    
  }
}

get.leg <- ggplot(data = df.omega, aes(x = Var1, y = Var2, fill = value)) + 
  geom_tile(color = "gray90", size = .8) +
  scale_fill_gradient2(limit = c(-1, 1)) + labs(x = NULL, y = NULL) +
  theme(axis.text.x = element_blank(), panel.background = element_blank(),
        title = element_text(size = 10))

legend <- get_legend(
  # create some space to the left of the legend
  get.leg + theme(legend.box.margin = margin(0, 0, 0, 12)))

cowplot::plot_grid(plotlist = p.omega$PA)

save_plot(here::here("output/Final_figures/FigS5_PA_omega.png"), 
          cowplot::plot_grid(p.omega$PA[[1]], p.omega$PA[[3]],
                             p.omega$PA[[5]], legend, ncol = 2,
                             rel_widths = c(2,2,2,1)),
          base_width = 10, base_height = 6)

cowplot::plot_grid(plotlist = p.omega$Normal)

save_plot(here::here("output/Final_figures/FigS6_Normal_omega.png"), 
          cowplot::plot_grid(p.omega$Normal[[1]], p.omega$Normal[[2]], 
                             p.omega$Normal[[3]], p.omega$Normal[[5]],
                             p.omega$Normal[[6]],
                             legend, ncol = 3, rel_widths = c(2,2,2,2,2,1)),
          base_width = 12, base_height = 6)


# Variance Partitioning
p.VarPart <- list(PA = vector("list", 7), Normal = vector("list", 7))

for (m in 1:length(all.models)) {
  for (i in 1:length(all.models[[m]])) {
    #m=2;i=1
    tmp <- result[[m]][[i]]
    tmp$vals <- tmp$vals[,match(all.models[[m]][[i]]$phyloTree$tip.label, colnames(tmp$vals))]
    
    p <- reshape2::melt(tmp$vals)
    
    prop.var <- round(rowMeans(result[[m]][[i]]$vals), 4)*100
    new.labs <- c(paste("Temperature (", prop.var[1], "%)", sep = ""),
                  paste("Forest (", prop.var[2], "%)", sep = ""),
                  paste("Grassland (", prop.var[3], "%)", sep = ""),
                  paste("Crops (", prop.var[4], "%)", sep = ""),
                  paste("Random (", prop.var[5], "%)", sep = ""))
    p$Var3 <- as.factor(paste(stringr::str_sub(matrix(unlist(strsplit(as.character(p$Var2), " ")), ncol = 2, byrow = T)[,1], 1,3),
                              stringr::str_sub(matrix(unlist(strsplit(as.character(p$Var2), " ")), ncol = 2, byrow = T)[,2], 1,3), sep = " "))
    
    levels(p$Var3) <- paste(stringr::str_sub(matrix(unlist(strsplit(levels(p$Var2), " ")), ncol = 2, byrow = T)[,1], 1,3),
                            stringr::str_sub(matrix(unlist(strsplit(levels(p$Var2), " ")), ncol = 2, byrow = T)[,2], 1,3), sep = " ")
    
    p.VarPart[[m]][[i]] <- ggplot(p, aes(fill = Var1, y = value, x = Var3)) + 
      geom_bar(position="fill", stat="identity") + labs(title = paste(comm.names[i])) +
      scale_fill_viridis_d(name = NULL, option = "A", direction = -1,
                           labels = new.labs) + coord_flip() + labs(x = NULL, y = NULL)+
      theme(panel.background = element_blank(), legend.position = "right",
            title = element_text(size = 10))
  }
}

p1 <- cowplot::plot_grid(plotlist = p.VarPart$PA, ncol = 4)
p1
save_plot(here::here("output/Final_figures/Fig2_PA_VP.png"), p1, 
          base_height = 10, base_width = 16)

p2 <- cowplot::plot_grid(plotlist = p.VarPart$Normal, ncol = 4)
p2
save_plot(here::here("output/Final_figures/Fig3_Normal_VP.png"), p2, 
          base_height = 10, base_width = 16)


# Accessing the proportion of explanation of traits on covariates and occurrence/abundance of species

df.R2T <- list(PA = matrix(NA, ncol = 6, nrow = 7), Normal = matrix(NA, ncol = 6, nrow = 7))
colnames(df.R2T$PA) <- names(unlist(result$PA$VP1$R2T))

for (m in 1:length(result)) {
  for (i in 1:7) {
    #m=1;i=1
    df.R2T[[m]][i,] <- unlist(result[[m]][[i]]$R2T)
    
  }
}
(df.R2T)
write.table(data.frame(Model = rep(names(df.R2T)[c(1,2)], each = 7), 
                       rbind(df.R2T$PA, df.R2T$Normal)),
            here::here("output/Final_tables/dfR2T_PA_Norm.txt"))
df.R2T <- read.table(here::here("output/Final_tables/dfR2T_PA_Norm.txt"), h = T)


# Community-level responses -----------------------------------------------

# We start by making gradient plots that visualize how the communities vary among the 
# environmental variables.
source(here::here("R/functions/plotGradient_modify_function.R"))

all.models[[1]][[1]]$XFormula
cbp <- c("#D16103", #EEA
         "#C4961A", #JG1
         "#E69F00", #JG2
         "#009E73", #QR1
         "#52854C", #QR2
         "#4E84C4", #SG1
         "#56B4E9") #SG2
scales::show_col(cbp)

# Focal: Mean Temperature
pred.Temp <- list(PA = vector("list", 14), 
                  Normal = vector("list", 14))

for (m in 1:length(all.models)) {
  for (i in 1:length(all.models[[m]])) {
    #m=1;i=1
    Gradient <- constructGradient(all.models[[m]][[i]], focalVariable = "Temp",
                                  non.focalVariables = list("Grasslands" = list(1), "Forest" = list(1),
                                                            "Crops" = list(1)))
    predY <- predict(all.models[[m]][[i]], Gradient=Gradient, expected = TRUE)
    pred.Temp[[m]][[i]] <- plotGradient_modify(hM = all.models[[m]][[i]], Gradient = Gradient, 
                                               predY = predY, measure = "S")
    pred.Temp[[m]][[i+7]] <- plotGradient_modify(hM = all.models[[m]][[i]], Gradient = Gradient, 
                                               predY = predY, measure = "T")
    
  }
}


# Focal: Forest proportion
pred.Forest <- list(PA = vector("list", 14), 
                  Normal = vector("list", 14))

for (m in 1:length(all.models)) {
  for (i in 1:length(all.models[[m]])) {
    #m=2;i=1
    Gradient <- constructGradient(all.models[[m]][[i]], focalVariable = "Forest",
                                  non.focalVariables = list("Grasslands" = list(1), "Temp" = list(1),
                                                            "Crops" = list(1)))
    predY <- predict(all.models[[m]][[i]], Gradient=Gradient, expected = TRUE)
    pred.Forest[[m]][[i]] <- plotGradient_modify(hM = all.models[[m]][[i]], Gradient = Gradient, 
                                                 predY = predY, measure = "S")
    pred.Forest[[m]][[i+7]] <- plotGradient_modify(hM = all.models[[m]][[i]], Gradient = Gradient, 
                                                 predY = predY, measure = "T")
    
  }
}


# Focal: Grassland proportion
pred.Grass <- list(PA = vector("list", 14),  
                   Normal = vector("list", 14))

for (m in 1:length(all.models)) {
  for (i in 1:length(all.models[[m]])) {
    #m=1;i=1
    Gradient <- constructGradient(all.models[[m]][[i]], focalVariable = "Grasslands",
                                  non.focalVariables = list("Forest" = list(1), "Temp" = list(1),
                                                            "Crops" = list(1)))
    predY <- predict(all.models[[m]][[i]], Gradient=Gradient, expected = TRUE)
    pred.Grass[[m]][[i]] <- plotGradient_modify(hM = all.models[[m]][[i]], Gradient = Gradient, 
                                                predY = predY, measure = "S")
    pred.Grass[[m]][[i+7]] <- plotGradient_modify(hM = all.models[[m]][[i]], Gradient = Gradient, 
                                                predY = predY, measure = "T")
  }
}


#Focal: Crop proportion
pred.Crops <- list(PA = vector("list", 14),
                   Normal = vector("list", 14))

for (m in 1:length(all.models)) {
  for (i in 1:length(all.models[[m]])) {
    #m=1;i=1
    Gradient <- constructGradient(all.models[[m]][[i]], focalVariable = "Crops",
                                  non.focalVariables = list("Grasslands" = list(1), "Forest" = list(1),
                                                            "Temp" = list(1)))
    predY <- predict(all.models[[m]][[i]], Gradient=Gradient, expected = TRUE)
    pred.Crops[[m]][[i]] <- plotGradient_modify(hM = all.models[[m]][[i]], Gradient = Gradient, 
                                                predY = predY, measure = "S")
    pred.Crops[[m]][[i+7]] <- plotGradient_modify(hM = all.models[[m]][[i]], Gradient = Gradient, 
                                                predY = predY, measure = "T")
    
  }
}

comm_level <- list(pred.Temp, pred.Forest, pred.Grass, pred.Crops)
saveRDS(comm_level, here::here("output/Comm_responses.rds"))

# plots in the script "S_Comm_level"

# Community size vs effect of environmental conditions --------------------

MF <- read.table(here::here("output/Final_tables/ModelFit_all.txt"))
MF
MF.norm <- subset(MF, Model == "Normal")$qua_50 # usando a mediana
comm.size <- vector()
for(i in 1:length(Y)){
  comm.size[i] <- sum(Y[[i]])
}

df.size <- data.frame(r2 = MF.norm, size = comm.size, comm = rownames(MF)[1:7])
lm.r2 <- summary(lm(r2 ~ size, data = df.size))


library(ggplot2)
R2.size <- ggplot(df.size, aes(x = size, y = r2, label = comm)) +
  geom_point() + geom_smooth(method = "lm") +
  geom_text(hjust = 'inward', vjust = 'inward') + labs(y = expression(paste("Global ", R^2, sep = "")), x = "Community size") +
  annotate("text", x = 220, y = 0.83, label = paste(bquote("R^2 == "), round(lm.r2$r.squared, 2), sep = " "),
           parse = TRUE, size = 3)
R2.size


# evaluating the relationship between VP and community size
result <- readRDS(here::here("output/Final_tables/VP_and_fixParam.rds"))

prop.var <- as.data.frame(matrix(NA, ncol = 5, nrow = 7))
colnames(prop.var) <- rownames(result$Normal[[1]]$vals)
rownames(prop.var) <- rownames(MF)[1:7]

for (i in 1:7) { # 7 because is the number of communities
    prop.var[i,] <- round(rowMeans(result$Normal[[i]]$vals), 4)*100
}
prop.var$size <- comm.size
prop.var$landscape <- rowSums(prop.var[,2:4])
hist(prop.var$Temp)

lm.tmp <- summary(lm(Temp ~ size, data = prop.var))
lm.for <- summary(lm(Forest ~ size, data = prop.var))
lm.gra <- summary(lm(Grasslands ~ size, data = prop.var))
lm.cro <- summary(lm(Crops ~ size, data = prop.var))
lm.rnd <- summary(lm(`Random: SU` ~ size, data = prop.var))
lm.land <- summary(lm(landscape ~ size, data = prop.var))
reshape2::melt(prop.var)

temp.size <- ggplot(prop.var, aes(x = size, y = Temp, label = rownames(prop.var))) +
  geom_point() + geom_smooth(method = "lm") +
  geom_text(hjust = 'inward', vjust = 'inward') + labs(y = "Temperature", x = "Community size") +
  annotate("text", x = 220, y = 13, 
           label = paste(bquote("R^2 =="), round(lm.tmp$r.squared, 2), sep = " "),
           parse = TRUE, size = 3) 
temp.size

forest.size <- ggplot(prop.var, aes(x = size, y = Forest, label = rownames(prop.var))) +
  geom_point() + geom_smooth(method = "lm") +
  geom_text(hjust = 'inward', vjust = 'inward') + labs(y = "Forest", x = "Community size") +
  annotate("text", x = 220, y = 41, label = paste(bquote("R^2 =="), round(lm.for$r.squared, 2), sep = " "),
           parse = TRUE, size = 3) 
forest.size

grass.size <- ggplot(prop.var, aes(x = size, y = Grasslands, label = rownames(prop.var))) +
  geom_point() + geom_smooth(method = "lm") +
  geom_text(hjust = 'inward', vjust = 'inward') + labs(y = "Grassland", x = "Community size") +
  annotate("text", x = 220, y = 33, label = paste(bquote("R^2 =="), round(lm.gra$r.squared, 2), sep = " "),
           parse = TRUE, size = 3) 
grass.size

crop.size <- ggplot(prop.var, aes(x = size, y = Crops, label = rownames(prop.var))) +
  geom_point() + geom_smooth(method = "lm") +
  geom_text(hjust = 'inward', vjust = 'inward') + labs(y = "Crop", x = "Community size") +
  annotate("text", x = 220, y = 1, label = paste(bquote("R^2 =="), round(lm.cro$r.squared, 2), sep = " "),
           parse = TRUE, size = 3) 
crop.size

random.size <- ggplot(prop.var, aes(x = size, y = `Random: SU`, label = rownames(prop.var))) +
  geom_point() + geom_smooth(method = "lm") +
  geom_text(hjust = "inward", vjust = 'inward') + labs(y = "Random effect", x = "Community size") +
  annotate("text", x = 220, y = 6, label = paste(bquote("R^2 =="), round(lm.rnd$r.squared, 2), sep = " "),
           parse = TRUE, size = 3)
random.size 

landscape.size <- ggplot(prop.var, aes(x = size, y = landscape, label = rownames(prop.var))) +
  geom_point() + geom_smooth(method = "lm") +
  geom_text(hjust = "inward", vjust = 'inward') + labs(y = "Landscape", x = "Community size") +
  annotate("text", x = 220, y = 71, label = paste(bquote("R^2 =="), round(lm.rnd$r.squared, 2), sep = " "),
           parse = TRUE, size = 3)
landscape.size 

table.lm <- round(rbind(lm.r2$coefficients, lm.tmp$coefficients,
                        lm.for$coefficients, lm.gra$coefficients,
                        lm.cro$coefficients, lm.rnd$coefficients), 3)
table.lm

cowplot::save_plot(here::here("output/Final_figures/Fig6_Comm_size.png"),
                   cowplot::plot_grid(R2.size, temp.size, forest.size,
                                      grass.size, crop.size, random.size,  
                                      labels = paste(letters[1:6], ")", sep = ""),
                                      label_fontface = "plain", label_size = 12),
                   base_height = 6, base_width = 12)
