#'---
#' title: "Data analysis"
#' author: "Larkin, Beck, & Bajer"
#' output: 
#'    html_document:
#'       toc: true
#'       theme: default
#'       toc_depth: 2
#'       toc_float:
#'           collapsed: false
#'--- 

#' ### Load libraries
#+ message=FALSE, warning=FALSE
library(vegan)
library(GGally)
library(mvabund)
library(dplyr)
library(mosaic)
library(ggplot2)
library(stringr)
library(Hmisc)
library(onewaytests)
library(FD)
library(tibble)
library(sjPlot)
library(ggeffects)
library(grid)
library(gridExtra)
options(width=160)

#' # Data prep
#' 
#' Read in and organize data
#' 
#' - "speciesTraits" = Native/non-native status, perennial/annual, growth form, 
#' and C-value for 206 taxa
speciesTraits <- read.csv("data/speciesTraits.csv")
#' 
#' - "allData_913" = environmental, carp, and macrophyte data for all 913 lakes
allData_913 <- read.csv("data/allData_913.csv")

#' - "commData_913" = macrophyte community data only (pres/abs) for all 913 lakes (206 taxa)
commData_913 <- allData_913
cols.keep <- which(colnames(commData_913) == "acorus_calamus_var._americanus"):dim(commData_913)[2] 
commData_913 <- commData_913[,cols.keep] 
cols.drop <- which(colSums(commData_913) == 0) 
commData_913 <- commData_913[,-cols.drop]
dim(commData_913)

#' Standardizing environmental predictor variables to mean 0, SD 1 
stdVar <- function(x){(x - mean(x)) / sd(x)}

allData_913$std.log.depth <- stdVar(log(allData_913$depthm + 1))
allData_913$std.log.carp <- stdVar(log(allData_913$cpue.carp + 1))
allData_913$std.phuman <- stdVar(allData_913$phuman)
allData_913$std.invSecchi <- stdVar(allData_913$invSecchi)
allData_913$std.sdi <- stdVar(allData_913$sdi)
allData_913$std.log.area <- stdVar(log(allData_913$aream2))

#' Reducing dataset by dropping rare taxa (occurring in <2.5% of lakes). This reduces dataset to 
#' 119 taxa, and 906 lakes w/ >=2 of the retained taxa
#' 
#' - "allData_906" = environmental, carp, and macrophyte data for 906 lakes
#' - "commData_906" = macrophyte community data for 906 lakes
commData_906 <- commData_913
numbRecs <- apply(commData_906, c(2), length)
lengthNotZero.fun <- function(x) length(which(x!=0))
numbNotZero <- apply(commData_906, c(2), lengthNotZero.fun)
sppFreq <- numbNotZero/numbRecs
commData_906 <- commData_906[, sppFreq >= 0.025]
drop.rows <- which(rowSums(commData_906) < 2)
commData_906 <- commData_906[-drop.rows,]
allData_906 <- allData_913
allData_906 <- allData_913[-drop.rows,]
dim(commData_906)

#' # NMDS ordination
#+ nmds, cache=TRUE, message=FALSE, results="hide"
set.seed(17)
commData_906.nms <- metaMDS(commData_906, k = 3, distance = "bray")
#'
commData_906.nms

#' ## Figure 1
#' 
#' See R file for code
#+ nmds.fig, cache=TRUE, echo=FALSE
tiff(file = "Fig.1.NMDS.tif", width = 6, height = 5.5, units = "in", res = 600,
     compression = "lzw")
par(mar=c(5, 4.5, 0.5, 1))
fig <- ordiplot(commData_906.nms, type="none", ylim = c(-0.8, 1.4))
points(fig, "sites", pch=21, bg="white",
       cex=1*as.numeric(allData_906$ecoregion == "eastern temperate forests")
       *(log(allData_906$cpue.carp + 1) + 0.5))
points(fig, "sites", pch=21, bg="gray60",
       cex=1*as.numeric(allData_906$ecoregion == "great plains")
       *(log(allData_906$cpue.carp + 1) + 0.5))
points(fig, "sites", pch=21, bg="black",
       cex=1*as.numeric(allData_906$ecoregion == "northern forests")
       *(log(allData_906$cpue.carp + 1) + 0.5))
legend.text <- c("N Forests","E Temp Forests","G Plains")
legend(0, 1.7, title = "Ecoregion", legend.text, bty = "n", pch = 21,
       pt.bg = c("black","white","gray60"), pt.cex = 1, cex = 0.8,
       title.adj = 0)
legend(0.85, 1.7, title = "Carp (CPUE)", c("0","1","10", "50"),
       bty = "n", pch = 21, pt.cex = c(0.5, log(1+1)+0.5, log(10+1)+0.5, log(50+1)+0.5),
       cex = 0.8, title.adj = 0, x.intersp = 1.5)
dev.off()

#' # Multivariate GLM
#' 
#' Predictor variables
#' 
#' - Ecoregion = 'ecoregion'
#' - Sampling year = 'fish.year'
#' - Sampling season = 'season'
#' - Lake area = 'aream2' (ln transformed)
#' - Maximum depth = 'depthm' (ln transformed)
#' - Inverse Secchi depth = 'invSecchi'
#' - Proportion anthropogenic land use = 'phuman'
#' - Lake morphometry (shoreline development index) = 'SDI'
#' - Carp CPUE = 'cpue.carp' (ln(x+1) transformed)
#' 

#' Checking for excess collinearity between predictor variables
#' 
#' All pairwise correlations <= 0.52. Reasonable to retain all predictors
#+ ggpairs, cache=TRUE
ggpairs(allData_906[, c(238:243)], lower = list(continuous = wrap("points", alpha = 0.1, size = 1)))

#' ## All ecoregions

#' Fitting full multivariate glm w/ interactions
#+ commData_906.glm.full, cache=TRUE
commData_906.mv <- mvabund(commData_906)
commData_906.glm.full <- manyglm(commData_906.mv ~ ecoregion*std.log.carp + std.invSecchi*std.log.carp
                            + std.log.area + std.log.depth + fish.year + std.phuman + std.sdi 
                            + season, data = allData_906, family=binomial("cloglog"))
plot(commData_906.glm.full, overlay = TRUE, cex = 0.1, pch = 1)

#' Analysis of deviance for mv glm
#+ mv.adev.full, cache=TRUE, warning=FALSE, dependson=commData_906.glm.full
nboot <- 1000
commData_906.glm.full.adev <- anova.manyglm(commData_906.glm.full, nBoot=nboot)
commData_906.glm.full.adev

#' Reducing model (backwards selection)
#' 
#' Dropping non-significant ecoregion x carp interaction
#+ mv.adev.reduced, cache=TRUE, warning=FALSE
commData_906.glm <- manyglm(commData_906.mv ~ ecoregion + std.log.carp + std.invSecchi*std.log.carp
                            + std.log.area + std.log.depth + fish.year + std.phuman + std.sdi 
                            + season, data = allData_906, family=binomial("cloglog"))
commData_906.glm.adev <- anova.manyglm(commData_906.glm, nBoot=nboot)
commData_906.glm.adev

#' ## Excluding northern forests
#' 
#' Dropping lakes from Northern Forests ecoregion where carp rare
# Proportion of NF lakes w/ carp
sum(as.numeric(allData_913$cpue.carp[allData_913$ecoregion == "northern forests"] > 0)) /
  length(allData_913$cpue.carp[allData_913$ecoregion == "northern forests"])
drop.rows <- which(allData_913$ecoregion == "northern forests")
allData_NoNF <- allData_913[-drop.rows,]
commData_NoNF <- commData_913[-drop.rows,]

#' Removing species found in <2.5% of lakes, and then lakes w/ <2 species
numbRecs <- apply(commData_NoNF, c(2), length)
numbNotZero <- apply(commData_NoNF, c(2), lengthNotZero.fun)
sppFreq <- numbNotZero/numbRecs
commData_NoNF <- commData_NoNF[, sppFreq >= 0.025]
drop.rows <- which(rowSums(commData_NoNF) < 2)
commData_NoNF <- commData_NoNF[-drop.rows,]
allData_NoNF <- allData_NoNF[-drop.rows,]
#' Numbers of lakes and taxa
dim(commData_NoNF)

#' Fitting multivariate glm
#+ commData_NoNF.glm, cache=TRUE
commData_NoNF.mv <- mvabund(commData_NoNF)
commData_NoNF.glm <- manyglm(commData_NoNF.mv ~  ecoregion + std.log.carp + std.invSecchi*std.log.carp
                             + std.log.area + std.log.depth + fish.year + std.phuman + std.sdi + season, 
                             data = allData_NoNF, family=binomial("cloglog"))
plot(commData_NoNF.glm, overlay = TRUE, cex = 0.1, pch = 1)

#' Analysis of deviance for mv glm
#+ mv.adev.nonf, cache=TRUE, warning=FALSE, dependson=commData_NoNF.glm
commData_NoNF.glm.adev <- anova.manyglm(commData_NoNF.glm, nBoot=nboot)
commData_NoNF.glm.adev

#' # Individual species responses
#' 
#' Excluding species w/ few records (<25) outside NF ecoregion (and thus 
#' minimal potential overlap w/ carp)
#+ commdata.uni, cache=TRUE
ecoregCounts <- aggregate(. ~ allData_906$ecoregion, commData_906, sum)
northernSpp <- names(which(colSums(ecoregCounts[-3,-1]) < 25))
commData_906.uni <- commData_906[ , !names(commData_906) %in% northernSpp]

#' Function for univariate GLM
#+ glm.fun, cache=TRUE
glm.fun <- function(x){
  glm(x ~ ecoregion + std.invSecchi*std.log.carp + std.log.area 
      + std.log.depth + fish.year + std.phuman + std.sdi + season, data = allData_906, 
      family=binomial("cloglog"))
}

#' Function to randomize carp CPUE data across lakes to generate expectations under null 
#' model of no carp effect
#+ glm.shuffle, cache=TRUE
glm.shuffle <- function(x){
  glm(x ~ ecoregion + std.invSecchi*shuffle(std.log.carp) + std.log.area 
      + std.log.depth + fish.year + std.phuman + std.sdi + season, data = allData_906, 
      family=binomial("cloglog"))
}

#' ## Example species
#+ glm.example, cache=TRUE, warning=FALSE, dependson=commdata.uni, dependson=glm.fun
j <- "myriophyllum_sibiricum"
test.glm <- glm.fun(commData_906.uni[,j])
#' Observed carp coefficient
test.stat <- c(coef(test.glm)["std.log.carp"], 
               coef(test.glm)["std.invSecchi:std.log.carp"]); test.stat

#' Distribution of expected carp coefficients under null hypothesis of no carp effect
#+ test.null, warning=FALSE, cache=TRUE, dependson=commdata.uni, dependson=glm.shuffle
onesided <- matrix(NA, nboot, 2)
for(i in 1:nboot){
  mod1 <- glm.shuffle(commData_906.uni[,j])
  onesided[i,] <- c(coef(mod1)["shuffle(std.log.carp)"], 
                    coef(mod1)["std.invSecchi:shuffle(std.log.carp)"])
}

#' One-sided test of probability that observed carp coefficient < null expectation
#' 
#' Carp main effect
sum(onesided[,1] <= test.stat[1])/nboot
#' Carp x turbidity interaction
sum(onesided[,2] <= test.stat[2])/nboot

#' ## All species
#' 
#' Observed carp coefficients for all species
#+ carp.coef, cache=TRUE, warning=FALSE, dependson=commdata.uni, dependson=glm.fun
species <- c(1:dim(commData_906.uni)[2])
carp.coef <- data.frame()
for(j in species){
  glm.j <- glm.fun(commData_906.uni[,j])
  carp.coef.j <- data.frame(row.names = names(commData_906.uni)[j], 
                            carp.coef = summary(glm.j)$coefficients["std.log.carp","Estimate"],
                            carp.se = summary(glm.j)$coefficients["std.log.carp","Std. Error"],
                            carp.turb.coef = summary(glm.j)$coefficients["std.invSecchi:std.log.carp","Estimate"],
                            carp.turb.se = summary(glm.j)$coefficients["std.invSecchi:std.log.carp","Std. Error"])
  carp.coef <- rbind(carp.coef, carp.coef.j)
}

#' Loop to generate null expectations for each species 
#+ carp.null, cache=TRUE, warning=FALSE, dependson=commdata.uni, dependson=glm.shuffle
set.seed(23)
carp.null <- vector(length(species), mode = "list")
names(carp.null) <- names(commData_906.uni)[species]
for(j in species){
  onesided <- matrix(NA, nboot, 2)
  for(i in 1:nboot){
    glm.shuffle.j <- glm.shuffle(commData_906.uni[,j])
    onesided[i,] <- cbind(coef(glm.shuffle.j)["shuffle(std.log.carp)"], 
                          coef(glm.shuffle.j)["std.invSecchi:shuffle(std.log.carp)"])
  }
  carp.null[[j]] <- data.frame(onesided)
  names(carp.null[[j]]) <- c("carp.main", "carp.interaction")
}

#' Calculating p-value for each species
#+ carp.coef.p, cache=TRUE, dependson=carp.coef, dependson=carp.null
carp.coef.p <- carp.coef
carp.coef.p <- transform(carp.coef.p, carp.p = 0, carp.turb.p = 0)
for(j in species){
  carp.coef.p[j,5] <- sum(unlist(carp.null[[j]][,"carp.main"]) <= carp.coef[j,"carp.coef"])/nboot
  carp.coef.p[j,6] <- sum(unlist(carp.null[[j]][,"carp.interaction"]) <= carp.coef[j,"carp.turb.coef"])/nboot
}
head(carp.coef.p)

#' ## Carp effect vs. turbidity effect
#' 
#' Loop for inverse Secchi coefficients for all species
#+ warning=FALSE
species <- c(1:dim(commData_906.uni)[2])
turb.coef <- data.frame()
for(j in species){
  glm.j <- glm.fun(commData_906.uni[,j])
  turb.coef.j <- data.frame(row.names = names(commData_906.uni)[j], 
                            turb.coef = summary(glm.j)$coefficients["std.invSecchi","Estimate"],
                            turb.se = summary(glm.j)$coefficients["std.invSecchi","Std. Error"])
  turb.coef <- rbind(turb.coef, turb.coef.j)
}
head(turb.coef)

#' Paired t-test
t.test(carp.coef$carp.coef, turb.coef$turb.coef, paired = T)

#' # Species characteristics
#' 
#' Ordering by significance and magnitude of carp effects
#+ carpresponse, cache=TRUE, dependson=carp.coef.p
carpResponse <- carp.coef.p[order(carp.coef.p$carp.p, carp.coef.p$carp.coef),]
sig <- ifelse(carpResponse$carp.coef < 0 & carpResponse$carp.p <= 0.05 | 
                carpResponse$carp.coef < 0 & carpResponse$carp.turb.p <= 0.05, 
              1, 0)
carpResponse <- transform(carpResponse, sig = sig, Carp_response = 0)
#' Taxa w/ sig. negative carp coefficients "sensitive", else "neutral"
carpResponse$Carp_response[carpResponse$sig == 1] <- "sensitive"
carpResponse$Carp_response[carpResponse$Carp_response==0] <- "neutral"
carpResponse$Species <- rownames(carpResponse)

#' Joining w/ species trait data
#+ warning=FALSE
responseTraits <- left_join(carpResponse, speciesTraits)
responseTraits %>%
  group_by(Carp_response) %>%
  summarise(n=n())

#' ## Figure 2
#' 
#' See R file for code
#+ echo=FALSE, warning=FALSE
SpeciesLabs <- as.character(responseTraits$Species)
SpeciesLabs <- str_split(SpeciesLabs, "_", simplify = TRUE)
SpeciesLabs <- str_trunc(SpeciesLabs, width = 6, ellipsis = "")
SpeciesLabs <- as.data.frame(SpeciesLabs)
SpeciesLabs <- paste(SpeciesLabs$V1, SpeciesLabs$V2, sep="_")
SpeciesLabs[SpeciesLabs=="persci_sp"] <- "persic_sp"
SpeciesLabs[SpeciesLabs=="fontin_"] <- "fontin_sp"
responseTraits$SpeciesLabs <- SpeciesLabs
responseTraits$SpeciesLabs <- capitalize(responseTraits$SpeciesLabs)
responseTraits$SpeciesLabs <- gsub("_", " ", responseTraits$SpeciesLabs)

responseTraits$Origins <- as.character(responseTraits$Origins)
responseTraits$Origins[is.na(responseTraits$Origins)] <- "Ambiguous"
responseTraits$Origins <- as.factor(responseTraits$Origins)
responseTraits$Origins <- factor(responseTraits$Origins, 
                                 levels(responseTraits$Origins)[c(3,2,1)])
levels(responseTraits$Aquatic.Form) <- c("B. Emergent", "C. Floating", "A. Submersed")
responseTraits$Aquatic.Form <- factor(responseTraits$Aquatic.Form, 
                                      levels(responseTraits$Aquatic.Form)[c(3,1,2)])
responseTraits$sigSymbol <- 0
responseTraits$sigSymbol[responseTraits$sig == 0] <- ""
responseTraits$sigSymbol[responseTraits$carp.p <= 0.10] <- ""
responseTraits$sigSymbol[responseTraits$carp.p <= 0.05] <- "*"
responseTraits$sigSymbol[responseTraits$carp.p <= 0.01] <- "**"
responseTraits$sigSymbol[responseTraits$carp.p <= 0.001] <- "***"
responseTraits$SpeciesLabs <- factor(responseTraits$SpeciesLabs, levels = responseTraits[order(-responseTraits$carp.coef), "SpeciesLabs"])

tiff(file = "Fig.2.CarpCoef.tif", width = 10, height = 5, units = "in", res = 600,
     compression = "lzw")
p <- ggplot(responseTraits) + 
  geom_segment(aes(x = carp.coef - carp.se, xend = carp.coef + carp.se, 
                   y = SpeciesLabs, yend = SpeciesLabs), lwd=0.3) +
  geom_point(size = 2, aes(x=carp.coef, y=SpeciesLabs, fill=Origins), shape = 21) +
  scale_fill_manual(values = c("white", "black", "gray60")) +
  labs(x = "Carp effect", y = "Species") +
  geom_vline(xintercept = 0, lwd=0.5) +
  geom_text(aes(y=SpeciesLabs, label=sigSymbol), x=-5.7, hjust = 0, vjust = 0.8, size = 4, 
            col = "grey30") +
  coord_cartesian(xlim = c(-5.5, 0.5)) +
  scale_x_continuous(breaks = seq(-5, 0, 1)) +
  theme(axis.ticks.y = element_blank(), axis.text.y = element_text(face = "italic"))
p + facet_wrap(~ Aquatic.Form, scales="free_y", labeller = )
dev.off()

#' ## Ecoregion prediction data
#' 
#' Ecoregional mean environmental covariates and carp values (0/mean/max)
#' 
#' Northern Forests
nForests <- droplevels(subset(allData_913, ecoregion == "northern forests"))
nForests.newdata <- data.frame(ecoregion = "northern forests", 
                               std.log.area = mean(nForests$std.log.area),
                               std.log.depth = mean(nForests$std.log.depth),
                               std.sdi = mean(nForests$std.sdi), 
                               fish.year = mean(nForests$fish.year), 
                               std.invSecchi = mean(nForests$std.invSecchi), 
                               std.phuman = mean(nForests$std.phuman),
                               std.log.carp = mean(nForests$std.log.carp),
                               season = "mid_sum")
nForests.newdata <- nForests.newdata[rep(seq_len(nrow(nForests.newdata[1])), each=3),]
nForests.newdata$std.log.carp[1] <- min(nForests$std.log.carp)
nForests.newdata$std.log.carp[3] <- max(nForests$std.log.carp)

#' Eastern Temperate Forests (same as above, see R file for code)
#+ echo=FALSE, warning=FALSE
eForests <- droplevels(subset(allData_913, ecoregion == "eastern temperate forests"))
eForests.newdata <- data.frame(ecoregion = "eastern temperate forests", 
                               std.log.area = mean(eForests$std.log.area),
                               std.log.depth = mean(eForests$std.log.depth),
                               std.sdi = mean(eForests$std.sdi), 
                               fish.year = mean(eForests$fish.year), 
                               std.invSecchi = mean(eForests$std.invSecchi), 
                               std.phuman = mean(eForests$std.phuman),
                               std.log.carp = mean(eForests$std.log.carp),
                               season = "mid_sum")
eForests.newdata <- eForests.newdata[rep(seq_len(nrow(eForests.newdata[1])), each=3),]
eForests.newdata$std.log.carp[1] <- min(eForests$std.log.carp)
eForests.newdata$std.log.carp[3] <- max(eForests$std.log.carp)

#' Great Plains (same as above, see R file for code)
#+ echo=FALSE, warning=FALSE
gPlains <- droplevels(subset(allData_913, ecoregion == "great plains"))
gPlains.newdata <- data.frame(ecoregion = "great plains", 
                              std.log.area = mean(gPlains$std.log.area),
                              std.log.depth = mean(gPlains$std.log.depth),
                              std.sdi = mean(gPlains$std.sdi), 
                              fish.year = mean(gPlains$fish.year), 
                              std.invSecchi = mean(gPlains$std.invSecchi), 
                              std.phuman = mean(gPlains$std.phuman),
                              std.log.carp = mean(gPlains$std.log.carp),
                              season = "mid_sum")
gPlains.newdata <- gPlains.newdata[rep(seq_len(nrow(gPlains.newdata[1])), each=3),]
gPlains.newdata$std.log.carp[1] <- min(gPlains$std.log.carp)
gPlains.newdata$std.log.carp[3] <- max(gPlains$std.log.carp)

#' Untransformed carp CPUE mean and max. values 
gPlains.carp <- droplevels(subset(allData_913, ecoregion == "great plains"))
c(mean(gPlains.carp$cpue.carp), max(gPlains.carp$cpue.carp))
eForests.carp <- droplevels(subset(allData_913, ecoregion == "eastern temperate forests"))
c(mean(eForests.carp$cpue.carp), max(eForests.carp$cpue.carp))
nForests.carp <- droplevels(subset(allData_913, ecoregion == "northern forests"))
c(mean(nForests.carp$cpue.carp), max(nForests.carp$cpue.carp))

#' Generating prediction data for exemplar spp
#' 
#' All introduced species and example native species
plotSpecies <- c("myriophyllum_spicatum", "myriophyllum_sibiricum", "utricularia_intermedia", 
                 "bidens_beckii", "stuckenia_pectinata", "potamogeton_crispus", 
                 "potamogeton_amplifolius","potamogeton_compressus", "potamogeton_gramineus", 
                 "potamogeton_natans", "typha_angustifolia", "typha_latifolia", 
                 "zizania_palustris", "sagittaria_rigida", "bolboschoenus_fluviatilis", 
                 "phalaris_arundinacea", "lythrum_salicaria", "lemna_trisulca", "wolffia_sp", 
                 "nymphaea_odorata")

#' Eastern temperate forests, all covariates set to mean values while varying carp
eForests.log.carp.max <- log(max(allData_906$cpue.carp[allData_906$ecoregion=="eastern temperate forests"]))
eForests.newdata.100 <- eForests.newdata
eForests.newdata.100 <- eForests.newdata.100[rep(seq_len(nrow(eForests.newdata.100[1])), each=100),]
eForests.newdata.100$log.carp <- seq(0, eForests.log.carp.max, len=100)

#' Predicting probability of occurrence on new carp data
#+ warning=FALSE
allData_906$log.carp <- log(allData_906$cpue.carp + 1)
glm.fun2 <- function(x){
  glm(x ~ ecoregion + std.invSecchi*log.carp + std.log.area 
      + std.log.depth + fish.year + std.phuman + std.sdi + season, data = allData_906, 
      family=binomial("cloglog"))
}
glm.logcarp <- vector(length(species), mode = "list")
names(glm.logcarp) <- names(commData_906.uni)[species]
for(i in species){
  glm.logcarp[[i]] <- glm.fun2(commData_906.uni[,i])
}

eForests.pred <- vector(length(plotSpecies), mode = "list")
names(eForests.pred) <- plotSpecies
for(i in plotSpecies){
  eForests.pred[[i]] <- predict(glm.logcarp[[i]], newdata = eForests.newdata.100, type = "response", se.fit = T)
  eForests.pred[[i]] <- cbind(eForests.newdata.100, data.frame(eForests.pred[[i]])[,1:2], species = i)
}
eForests.pred <- do.call("rbind", eForests.pred)

#' ## Figure 3
#' 
#' See R file for code
#+ fig3, cache=TRUE, echo=FALSE, warning=FALSE, fig.height=6
# Cleaning up and ordering species names
species <- (as.character(eForests.pred$species))
facetLabs <- capitalize(as.character(eForests.pred$species))
facetLabs <- gsub("_", " ", facetLabs)
levels <- unique(facetLabs)
eForests.pred$species.new <- facetLabs
eForests.pred$species.new <- factor(eForests.pred$species.new, levels = levels)

dat_text <- data.frame(
  species.new = unique(eForests.pred$species.new),
  Species = unique(eForests.pred$species),
  x     = 2.5,
  y     = Inf)
dat_text <- left_join(dat_text, responseTraits[,c(8,9,11,13)])
dat_text$Carp_response <- capitalize(dat_text$Carp_response)
dat_text$Aquatic.Form <- gsub("A. ", "", dat_text$Aquatic.Form)
dat_text$Aquatic.Form <- gsub("B. ", "", dat_text$Aquatic.Form)
dat_text$Aquatic.Form <- gsub("C. ", "", dat_text$Aquatic.Form)
dat_text$fontface <- as.numeric(dat_text$Origins)
dat_text$color <- dat_text$fontface
dat_text$fontface[dat_text$fontface==2] <- "bold"
dat_text$fontface[dat_text$fontface==1] <- "plain"
dat_text$color[dat_text$color==2] <- "red"
dat_text$color[dat_text$color==1] <- "black"

dat_text$fontface2 <- as.numeric(as.factor(dat_text$Carp_response))
dat_text$color2 <- dat_text$fontface2
dat_text$fontface2[dat_text$fontface2==2] <- "bold"
dat_text$fontface2[dat_text$fontface2==1] <- "plain"
dat_text$color2[dat_text$color2==2] <- "blue"
dat_text$color2[dat_text$color2==1] <- "black"
dat_text$x[dat_text$Species=="potamogeton_crispus"] <- 0.25

tiff(file = "Fig.3.SppGlm.tif", width = 10, height = 6, units = "in", res = 600,
     compression = "lzw")
ticks <- c(1, 2, 5, 13, 41)
p <- ggplot(eForests.pred, aes(x = log.carp, y = fit)) +
  geom_ribbon(aes(ymin = fit - se.fit, ymax = fit + se.fit), fill = 'grey', alpha = 0.6) +
  geom_line() +
  labs(y = "Probability of occurrence") +
  scale_x_continuous(name="Carp CPUE", limits=c(0, max(eForests.pred$log.carp)),
                     breaks=log(ticks), labels = ticks-1) +
  theme(strip.text = element_text(face = "italic"), 
        strip.background = element_rect(fill="white")) +
  geom_text(data = dat_text, mapping = aes(x = x, y = y, label = Origins, 
                                           fontface = fontface), 
            hjust=0, vjust=1.5, size=2.5, color = dat_text$color) +
  geom_text(data = dat_text, mapping = aes(x = x, y = y, label = Aquatic.Form), 
            hjust=0, vjust=3.25, size=2.5) +
  geom_text(data = dat_text, mapping = aes(x = x, y = y, label = Carp_response, 
                                           fontface = fontface2), 
            hjust=0, vjust=5, size=2.5, color = dat_text$color2)
p + facet_wrap( ~ species.new, nrow = 4, scales = "free_y")
dev.off()

#' ## Differences in sensitivity
#' 
#' Native/Introduced (significant)
Origins <- with(responseTraits, table(Carp_response, Origins)); Origins
fisher.test(Origins[,-3])

#' Growth form (significant)
Aquatic.Form <- with(responseTraits, table(Carp_response, Aquatic.Form)); Aquatic.Form
fisher.test(Aquatic.Form)

#' Annual/Perennial (non-significant)
Duration <- with(responseTraits, table(Carp_response, Duration)); Duration
fisher.test(Duration)

#' C-values (significant)
C.sens <- responseTraits$C.value[responseTraits$Carp_response == "sensitive" 
                                 & is.na(responseTraits$C.value) == FALSE]
C.neutr <- responseTraits$C.value[responseTraits$Carp_response == "neutral" 
                                  & is.na(responseTraits$C.value) == FALSE]
t.test(C.sens, C.neutr)
# Mean +/- S.E.
c(mean(C.neutr), sd(C.neutr)/sqrt(length(C.neutr)))
c(mean(C.sens), sd(C.sens)/sqrt(length(C.sens)))

#' ## Table S1
#'
TableA1 <- responseTraits
TableA1 <- TableA1[, c("Carp_response", "Species", "Origins", "Aquatic.Form", "Duration",
                       "C.value", "carp.coef", "carp.turb.coef")]
names(TableA1)[names(TableA1) %in% c("Aquatic.Form", "C.value")] <- c("Growth_form", "C_value")
TableA1$Carp_response<- capitalize(as.character(TableA1$Carp_response))
TableA1$Species <- capitalize(as.character(TableA1$Species))
TableA1$Growth_form <- capitalize(as.character(TableA1$Growth_form))
TableA1$Species <- gsub("_", " ", TableA1$Species)
TableA1 <- arrange(TableA1, desc(Carp_response), Species)
write.csv(TableA1, file = "TableA1.csv", row.names = FALSE)

#' # Community-level responses
#' 
#' ## Lake-level community metrics
#' 
#' Species richness
SR <- rowSums(commData_913)
summary(SR)

#' Mean conservatism, proportions of introduced and submersed species, calculated 
#' as community-weighted mean trait value using
#+ cwm, cache=TRUE
rownames(speciesTraits) <- speciesTraits$Species
speciesTraits <- speciesTraits[,-1]
speciesTraits <- speciesTraits[order(rownames(speciesTraits)),]
CWM <- functcomp(speciesTraits, as.matrix(commData_913), CWM.type = "all")
summary(CWM)[,c(1,2,8)]

#' ## Linear models
#' Species richness
#+ sr.lm, cache=TRUE, warning=FALSE, message=FALSE
SR.lm <- lm(log(SR) ~  ecoregion + std.invSecchi*std.log.carp + std.log.area 
            + std.log.depth + fish.year + std.phuman + std.sdi + season, 
            data = allData_913)
# R-squared
SR.sum <- summary(SR.lm)
SR.sum$adj.r.squared
# Organizing output for table
SR.sum <- as.data.frame(coef(SR.sum))
SR.sum <- data.frame(SR.sum[-c(1:3, 11:12), ])
SR.anova <- anova(SR.lm)
SR.anova <- as.data.frame(SR.anova)
SR.anova <- SR.anova[,-3]
SR.table <- left_join(rownames_to_column(SR.anova), rownames_to_column(SR.sum))
SR.table <- SR.table[, -c(7:9)]
SR.table <- SR.table[, c(1, 2, 6, 3:5)]
SR.table$Estimate <- round(SR.table$Estimate, digits=3)
SR.table$'Sum Sq' <- round(SR.table$'Sum Sq', digits=3)
SR.table$'F value' <- round(SR.table$'F value', digits=2)
SR.table$'Pr(>F)' <- round(SR.table$'Pr(>F)', digits=4)
write.csv(SR.table, "Table.SR.csv", row.names = F)
SR.table

#' Conservatism
#+ cv.lm, cache=TRUE, warning=FALSE, message=FALSE
CV.lm <- lm(CWM$C.value ~  ecoregion + std.invSecchi*std.log.carp + std.log.area 
            + std.log.depth + fish.year + std.phuman + std.sdi + season, 
            data = allData_913)
# R-squared
CV.sum <- summary(CV.lm)
CV.sum$adj.r.squared
# Organizing output for table
CV.sum <- as.data.frame(coef(CV.sum))
CV.sum <- data.frame(CV.sum[-c(1:3, 11:12), ])
CV.anova <- anova(CV.lm)
CV.anova <- as.data.frame(CV.anova)
CV.anova <- CV.anova[,-3]
CV.table <- left_join(rownames_to_column(CV.anova), rownames_to_column(CV.sum))
CV.table <- CV.table[, -c(7:9)]
CV.table <- CV.table[, c(1, 2, 6, 3:5)]
CV.table$Estimate <- round(CV.table$Estimate, digits=3)
CV.table$'Sum Sq' <- round(CV.table$'Sum Sq', digits=3)
CV.table$'F value' <- round(CV.table$'F value', digits=2)
CV.table$'Pr(>F)' <- round(CV.table$'Pr(>F)', digits=4)
write.csv(CV.table, "Table.Cvalue.csv", row.names = F)
CV.table

#' ## GLMs 
#' Proportion of species introduced
#+ intro.glm, cache=TRUE, warning=FALSE, message=FALSE
Intro.glm <- glm(CWM$Origins_Introduced ~  ecoregion + std.invSecchi*std.log.carp + std.log.area 
                 + std.log.depth + fish.year + std.phuman + std.sdi + season, 
                 data = allData_913, family = binomial, weights = SR)
# Organizing output for table
Intro.sum <- summary(Intro.glm)
Intro.sum <- as.data.frame(coef(Intro.sum))
Intro.sum <- data.frame(Intro.sum[-c(1:3, 11:12), ])
Intro.anova <- anova(Intro.glm, test = "Chisq") # significance test
Intro.anova <- as.data.frame(Intro.anova)
Intro.table <- left_join(rownames_to_column(Intro.anova), rownames_to_column(Intro.sum))
Intro.table <- Intro.table[, -c(8:10)]
Intro.table <- Intro.table[-1, c(1, 2, 7, 3:6)]
Intro.table$Estimate <- round(Intro.table$Estimate, digits=3)
Intro.table$Deviance <- round(Intro.table$Deviance, digits=2)
Intro.table$'Resid. Dev' <- round(Intro.table$'Resid. Dev', digits=2)
Intro.table$'Pr(>Chi)' <- round(Intro.table$'Pr(>Chi)', digits=4)
write.csv(Intro.table, "Table.PropIntro.csv", row.names = F)
Intro.table

#' Proportion of species submersed
#+ subm.glm, cache=TRUE, warning=FALSE, message=FALSE
Subm.glm <- glm(CWM$Aquatic.Form_submersed ~  ecoregion + std.invSecchi*std.log.carp + std.log.area 
                + std.log.depth + fish.year + std.phuman + std.sdi + season, 
                data = allData_913, family = binomial, weights = SR)
# Organizing output for table
Subm.sum <- summary(Subm.glm)
Subm.sum <- as.data.frame(coef(Subm.sum))
Subm.sum <- data.frame(Subm.sum[-c(1:3, 11:12), ])
Subm.anova <- anova(Subm.glm, test = "Chisq") # significance test
Subm.anova <- as.data.frame(Subm.anova)
Subm.table <- left_join(rownames_to_column(Subm.anova), rownames_to_column(Subm.sum))
Subm.table <- Subm.table[, -c(8:10)]
Subm.table <- Subm.table[-1, c(1, 2, 7, 3:6)]
Subm.table$Estimate <- round(Subm.table$Estimate, digits=3)
Subm.table$Deviance <- round(Subm.table$Deviance, digits=2)
Subm.table$'Resid. Dev' <- round(Subm.table$'Resid. Dev', digits=2)
Subm.table$'Pr(>Chi)' <- round(Subm.table$'Pr(>Chi)', digits=4)
write.csv(Subm.table, "Table.PropSubm.csv", row.names = F)
Subm.table

#' Probability lake invaded
#+ invaded.glm, cache=TRUE, warning=FALSE, message=FALSE
Invaded <- ceiling(CWM$Origins_Introduced) #does lake have at least one introduced spp?
Invaded.glm <- glm(Invaded ~  ecoregion + std.invSecchi*std.log.carp + std.log.area 
                   + std.log.depth + fish.year + std.phuman + std.sdi + season, 
                   data = allData_913, family = binomial)
# Organizing output for table
Invaded.sum <- summary(Invaded.glm)
Invaded.sum <- as.data.frame(coef(Invaded.sum))
Invaded.sum <- data.frame(Invaded.sum[-c(1:3, 11:12), ])
Invaded.anova <- anova(Invaded.glm, test = "Chisq") # significance test
Invaded.anova <- as.data.frame(Invaded.anova)
Invaded.table <- left_join(rownames_to_column(Invaded.anova), rownames_to_column(Invaded.sum))
Invaded.table <- Invaded.table[, -c(8:10)]
Invaded.table <- Invaded.table[-1, c(1, 2, 7, 3:6)]
Invaded.table$Estimate <- round(Invaded.table$Estimate, digits=3)
Invaded.table$Deviance <- round(Invaded.table$Deviance, digits=2)
Invaded.table$'Resid. Dev' <- round(Invaded.table$'Resid. Dev', digits=2)
Invaded.table$'Pr(>Chi)' <- round(Invaded.table$'Pr(>Chi)', digits=4)
write.csv(Invaded.table, "Table.ProbInvaded.csv", row.names = F)
Invaded.table

#' # Predicted community effects 
#' 
#' ## Linear models 
#' Following Vittinghoff et al. (2011), pp. 128-129
#' 
#' Species richness
# log-transformed response (SR) and predictor (carp) 
pct <- 1.10 #10% change in predictor
SR.coef <- SR.table[3,3]
SR.resp <- 100 * (exp(SR.coef * log(pct)) - 1)
SR.ci <- confint(SR.lm)
SR.ci <- as.vector(SR.ci[which(rownames(SR.ci)=="std.log.carp"),1:2])
SR.lo <- min(SR.ci)
SR.lo <- 100 * (exp(SR.lo * log(pct)) - 1)
SR.hi <- max(SR.ci)
SR.hi <- 100 * (exp(SR.hi * log(pct)) - 1)
c(SR.resp, SR.lo, SR.hi)

#' C-value
# untransformed response (C-value) and log-transformed predictor (carp) 
pct <- 1.10 #10% change in predictor
CV.coef <- CV.table[3,3]
# decrease in C-value for each 'pct' increase in carp
CV.resp <- CV.coef * log(pct)
CV.ci <- confint(CV.lm)
CV.ci <- as.vector(CV.ci[which(rownames(CV.ci)=="std.log.carp"),1:2])
CV.lo <- min(CV.ci)
CV.lo <- CV.lo * log(pct)
CV.hi <- max(CV.ci)
CV.hi <- CV.hi * log(pct)
c(CV.resp, CV.lo, CV.hi)

#' ## GLMs 
#' 
#' Order of results: 1. carp = 0, 1.1 carp = ecoregion mean, 1.2. carp = ecoregion max.
#'
#' Invaded.glm
#' 
#' northern forests
100*(predict.glm(Invaded.glm, nForests.newdata, type = "response"))
#' eastern temperate forests
100*(predict.glm(Invaded.glm, eForests.newdata, type = "response"))
#' great plains
100*(predict.glm(Invaded.glm, gPlains.newdata, type = "response"))

#' Subm.glm
#' 
#' northern forests
100*(predict.glm(Subm.glm, nForests.newdata, type = "response"))
#' eastern temperate forests
100*(predict.glm(Subm.glm, eForests.newdata, type = "response"))
#' great plains
100*(predict.glm(Subm.glm, gPlains.newdata, type = "response"))

#' Intro.glm
#' 
#' northern forests
100*(predict.glm(Intro.glm, nForests.newdata, type = "response"))
#' eastern temperate forests
100*(predict.glm(Intro.glm, eForests.newdata, type = "response"))
#' great plains
100*(predict.glm(Intro.glm, gPlains.newdata, type = "response"))

#' ## Figure 4
#' 
#' See R file for code
#+ fig4, cache=TRUE, message=FALSE, echo=FALSE
# Changing ecoregion orders for plotting
allData_913$ecoregion <- capitalize(as.character(allData_913$ecoregion))
allData_913$ecoregion[allData_913$ecoregion=="Eastern temperate forests"] <- "Eastern temp. forests"
allData_913$ecoregion <- factor(allData_913$ecoregion, 
                              levels = c("Northern forests", "Eastern temp. forests", "Great plains"))
# Models updated with log.carp as predictor 
SR.lm2 <- lm(log(SR) ~  ecoregion + std.invSecchi*log(cpue.carp+1) + std.log.area 
             + std.log.depth + fish.year + std.phuman + std.sdi + season, 
            data = allData_913)
CV.lm2 <- lm(CWM$C.value ~  ecoregion + std.invSecchi*log(cpue.carp+1) + std.log.area 
             + std.log.depth + fish.year + std.phuman + std.sdi + season, 
            data = allData_913)
Invaded.glm2 <- glm(Invaded ~  ecoregion + std.invSecchi*log(cpue.carp+1) + std.log.area 
                    + std.log.depth + fish.year + std.phuman + std.sdi + season, 
                   data = allData_913, family = binomial)
Intro.glm2 <- glm(CWM$Origins_Introduced ~  ecoregion + std.invSecchi*log(cpue.carp+1) + std.log.area 
                  + std.log.depth + fish.year + std.phuman + std.sdi + season, 
                 data = allData_913, family = binomial, weights = SR)
Subm.glm2 <- glm(CWM$Aquatic.Form_submersed ~  ecoregion + std.invSecchi*log(cpue.carp+1) + std.log.area 
                 + std.log.depth + fish.year + std.phuman + std.sdi + season, 
                data = allData_913, family = binomial, weights = SR)
# Creating prediction plots 
set_theme(plot.margin = unit(c(-0.5, 0.5, 0, 0.5), "cm"))
Fig.SR <- plot_model(SR.lm2, type = "pred", show.se = TRUE, 
                     terms = c("cpue.carp [exp]", "ecoregion"), title = "", 
                     grid=TRUE, transform = "exp", colors="black",
                     axis.title = c("", "\nSpecies richness"))
Fig.CV <- plot_model(CV.lm2, type = "pred", show.se = TRUE, 
                     terms = c("cpue.carp [exp]", "ecoregion"), title = "", 
                     grid=TRUE, transform = "exp", colors="black",
                     axis.title = c("", "\nMean C-value"))
Fig.Invaded <- plot_model(Invaded.glm2, digits=0, type = "pred", show.se = TRUE, 
                          terms = c("cpue.carp [exp]", "ecoregion"), title = "", 
                          grid=TRUE, show.intercept = FALSE, transform = "exp", colors="black", 
                          axis.title = c("", "Probability\nlake invaded"))
Fig.Intro <- plot_model(Intro.glm2, type = "pred", show.se = TRUE, 
                          terms = c("cpue.carp [exp]", "ecoregion"), title = "", 
                          grid=TRUE, transform = "exp", colors="black",
                          axis.title = c("", "Proportion\nintroduced species"))
Fig.Subm <- plot_model(Subm.glm2, type = "pred", show.se = TRUE, 
                        terms = c("cpue.carp [exp]", "ecoregion"), title = "", 
                        digits = 0, grid=TRUE, transform = "exp", colors="black",
                        axis.title = c("Carp CPUE", "Proportion\nsubmersed species"))
# Combining plots into single figure 
Fig.SR <- ggplotGrob(Fig.SR)
Fig.CV <- ggplotGrob(Fig.CV)
Fig.Invaded <- ggplotGrob(Fig.Invaded)
Fig.Intro <- ggplotGrob(Fig.Intro)
Fig.Subm <- ggplotGrob(Fig.Subm)
maxWidth = grid::unit.pmax(Fig.SR$widths[2:5], Fig.CV$widths[2:5], Fig.Invaded$widths[2:5],
                           Fig.Intro$widths[2:5], Fig.Subm$widths[2:5])
Fig.SR$widths[2:5] <- as.list(maxWidth)
Fig.CV$widths[2:5] <- as.list(maxWidth)
Fig.Invaded$widths[2:5] <- as.list(maxWidth)
Fig.Intro$widths[2:5] <- as.list(maxWidth)
Fig.Subm$widths[2:5] <- as.list(maxWidth)
tiff(file = "Fig.4.CommConseq.tif", width = 7, height = 10, units = "in", res = 600,
     compression = "lzw")
grid.arrange(Fig.SR, Fig.CV, Fig.Invaded, Fig.Intro, Fig.Subm, nrow = 5, 
             heights = c(1, 1, 1, 1, 1))
grid.text(c("A.", "B.", "C.", "D.", "E."), x = 0.018, 
          y = c(0.99, 0.79, 0.59, 0.39, 0.19), vjust=1, hjust=0, 
          gp=gpar(fontface=1, fontsize=13))
dev.off()

#' # Additional analyses 
#' 
#' Is carp-invasive plant relationship an artifact of both increasing concurrently
#' over the 20-year span of the monitoring data?
#' 
#' Proportion introduced macrophytes did not increase over time
InvProp.year <- lm(CWM$Origins_Introduced ~  fish.year, data = allData_913)
summary(InvProp.year)
#' Carp CPUE did not increase over time
CarpCpue.year <- lm(log(cpue.carp + 1) ~  fish.year, data = allData_913)
summary(CarpCpue.year)

#' # Saving file outputs
save(list=ls(), file = "carp.RData")

#' ### R Information (version and packages)
devtools::session_info()