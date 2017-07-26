#Cassin's Sparrows

library(tidyverse)

CASP_behavior <- read.csv(file = "CASP/20170207_CASP_playback_responses_combined_longform.csv")
CASP_veg_daub <- read.csv(file = "CASP/combined_vegetation_percent_cover_CASP.csv")
CASP_veg_sp   <- read.csv(file = "CASP/combined_vegetation_surveys_CASP.csv")

#get everything under same column name to assign points
CASP_veg_daub$Pointnames <- CASP_veg_daub$Point
CASP_veg_sp$Pointnames <- CASP_veg_sp$Point

CASP_behavior <- unite(CASP_behavior,
      Pointnames,
      Transect, Point, 
      sep = "", 
      remove = TRUE)

#Summarize daub

CASP_veg_daub_sum <- CASP_veg_daub %>%
  group_by(Location, Pointnames) %>%
  summarize(mForb = mean(Forb),
            mGrass = mean(Grass),
            mDead = mean(Dead),
            mBare = mean(Bare),
            mLitter = mean(Litter),
            mOther = mean(OtherFromNotes),
            mSum = mean(Sum)
            ) %>%
  filter(mSum >=95 & mSum<=105) #filtering out any more than 5% off

behavior.veg1 <- left_join(CASP_behavior,
                           CASP_veg_daub_sum,
                           by = c("Location",
                                  "Pointnames"))

#summarize shrub/tree counts:
CASP_veg_sp$Sum.counts <- CASP_veg_sp$Yucca +
  CASP_veg_sp$Sage+
  CASP_veg_sp$Sandplum+
  CASP_veg_sp$Cholla+ 
  CASP_veg_sp$Tree+ 
  CASP_veg_sp$Other.shrub
CASP_veg_sp_sum <- CASP_veg_sp %>%
  group_by(Location, Pointnames, Height) %>%
  summarize(sumYucca = sum(Yucca,
                           na.rm = TRUE),
            sumSage = sum(Sage,
                          na.rm = TRUE),
            sumSandplum = sum(Sandplum,
                              na.rm = TRUE),
            sumCholla = sum(Cholla,
                            na.rm = TRUE),
            sumTree = sum(Tree,
                          na.rm = TRUE),
            sumOther.shrub = sum(Other.shrub,
                                 na.rm = TRUE),
            sumCounts = sum(Sum.counts,
                            na.rm = TRUE)
  )

CASP_veg_sp_sum$Height <- gsub("<1m ",
                               "below1",
                               CASP_veg_sp_sum$Height,
                               fixed=TRUE)
CASP_veg_sp_sum$Height <- gsub(">1m",
                               "above1",
                               CASP_veg_sp_sum$Height,
                               fixed=TRUE)

CASP_veg_sp_sum$Height <- as.factor(CASP_veg_sp_sum$Height)
levels(CASP_veg_sp_sum$Height)

CASP_veg_sp_sum_spread <- CASP_veg_sp_sum %>% 
  gather(temp, score, starts_with("sum")) %>% 
  unite(temp1, Height, temp, sep = ".") %>% 
  spread(temp1, score)

behavior.veg <- left_join(behavior.veg1,
                          CASP_veg_sp_sum_spread,
                          by = c("Location",
                                 "Pointnames"))



#filter by points that have veg
behavior.veg <- behavior.veg %>%
  filter(!is.na(mSum))

veg.daub <- c("mForb",
              "mGrass",
              "mDead",
              "mBare",
              "mLitter",
              "mOther")

pca.veg.daub <- prcomp(behavior.veg[,
                                    veg.daub],
                       scale=TRUE, #correlation matrix used instead of covariance matrix, which is only appropriate if everything in same units.
                       retx=TRUE) #required to get PC scores for each individual.

summary(pca.veg.daub)
(pca.eigenvalues<-pca.veg.daub$sdev^2) #Get eigenvalues.
screeplot(pca.veg.daub) #Plot eigenvalues.
biplot(pca.veg.daub)

pca.veg.daub$rotation #eigenvectors.  Again the signs are arbitrary so don't worry
#if they differ but absolute values are the same between different programs or versions of R.
(loadings.pca<-cor(pca.veg.daub$x,
                   behavior.veg[,
                                veg.daub],
                   method="pearson"))
#Pearson's correlation of components with original variables.  Easier to interpret.
#Eigenvectors are how you get PCs, so also a sort of weight, just harder to think about.



pscores<-data.frame(pca.veg.daub$x) #puts PCA scores in a data frame

#Keeping those with eigenvalues above 1
behavior.veg$daubPC1 <- pscores$PC1
behavior.veg$daubPC2 <- pscores$PC2
behavior.veg$daubPC3 <- pscores$PC3


#counts of shrubs etc >1m tall
veg.above <- c("above1.sumYucca",
              "above1.sumSage",
              "above1.sumSandplum",
              "above1.sumCholla",
              "above1.sumTree",
              "above1.sumOther.shrub")

pca.veg.above <- prcomp(behavior.veg[,
                                    veg.above],
                       scale=TRUE, #correlation matrix used instead of covariance matrix, which is only appropriate if everything in same units.
                       retx=TRUE) #required to get PC scores for each individual.

summary(pca.veg.above)
(pca.eigenvalues.above<-pca.veg.above$sdev^2) #Get eigenvalues.
screeplot(pca.veg.above) #Plot eigenvalues.
biplot(pca.veg.above)

pca.veg.above$rotation #eigenvectors.  Again the signs are arbitrary so don't worry
#if they differ but absolute values are the same between different programs or versions of R.
(loadings.pca.above<-cor(pca.veg.above$x,
                   behavior.veg[,
                                veg.above],
                   method="pearson"))
#Pearson's correlation of components with original variables.  Easier to interpret.
#Eigenvectors are how you get PCs, so also a sort of weight, just harder to think about.



pscores.above<-data.frame(pca.veg.above$x) #puts PCA scores in a data frame

#Keeping those with eigenvalues above 1
behavior.veg$abovePC1 <- pscores.above$PC1
behavior.veg$abovePC2 <- pscores.above$PC2
behavior.veg$abovePC3 <- pscores.above$PC3

#shrubs below 1m
#counts of shrubs etc >1m tall
veg.below <- c("below1.sumYucca",
               "below1.sumSage",
               "below1.sumSandplum",
               "below1.sumCholla",
               "below1.sumTree",
               "below1.sumOther.shrub")

pca.veg.below <- prcomp(behavior.veg[,
                                     veg.below],
                        scale=TRUE, #correlation matrix used instead of covariance matrix, which is only appropriate if everything in same units.
                        retx=TRUE) #required to get PC scores for each individual.

summary(pca.veg.below)
(pca.eigenvalues.below<-pca.veg.below$sdev^2) #Get eigenvalues.
screeplot(pca.veg.below) #Plot eigenvalues.
biplot(pca.veg.below)

pca.veg.below$rotation #eigenvectors.  Again the signs are arbitrary so don't worry
#if they differ but absolute values are the same between different programs or versions of R.
(loadings.pca.below<-cor(pca.veg.below$x,
                         behavior.veg[,
                                      veg.below],
                         method="pearson"))
#Pearson's correlation of components with original variables.  Easier to interpret.
#Eigenvectors are how you get PCs, so also a sort of weight, just harder to think about.



pscores.below<-data.frame(pca.veg.below$x) #puts PCA scores in a data frame

#Keeping those with eigenvalues below 1
behavior.veg$belowPC1 <- pscores.below$PC1
behavior.veg$belowPC2 <- pscores.below$PC2
behavior.veg$belowPC3 <- pscores.below$PC3

#load libraries
library(lme4)
library(lmerTest)
#Strongest reaction by veg where present
lm.distance.veg <- lmer(ClosestDistance ~ daubPC1 + daubPC2 + daubPC3 +
                        abovePC1 + abovePC2 + abovePC3 +
                        abovePC1 + abovePC2 + abovePC3 + (1|Location),
                      data = behavior.veg[behavior.veg$Response==1,])

summary(lm.distance.veg)

#Presence/absence by veg

glm.presence.veg <- glmer(Response ~ daubPC1 + daubPC2 + daubPC3 +
                          abovePC1 + abovePC2 + abovePC3 +
                          belowPC1 + belowPC2 + belowPC3+ (1|Location),
                      data = behavior.veg,
                      family = "binomial")

summary(glm.presence.veg)

lm.strong.veg <- lmer(Strongest_behavior ~ daubPC1 + daubPC2 + daubPC3 +
                        abovePC1 + abovePC2 + abovePC3 +
                      belowPC1 + belowPC2 + belowPC3 + (1|Location),
                      data = behavior.veg)
summary(lm.strong.veg)


#Figures

fig2pred <- predict(glm.presence.veg)

ndFig2<- data.frame("latitude"=seq(min(Depredations$latitude),
                                max(Depredations$latitude),
                                length.out=length(Depredations$latitude)))
#plot the prediction with the new data (otherwise it uses rownumber and stretches the line out uselessly).
lines(ndW$latitude,
      predict.glm(glm.wolf,
                  newdata=ndW,
                  type="response"),
      lty="solid",lwd=2)

plot(Response ~ abovePC3,
     data = behavior.veg,
     xlab = "AbovePC3: decreasing yucca (-0.76), increasing trees (0.53)")

plot(Response ~ belowPC1,
     data = behavior.veg,
     xlab = "BelowPC1: decreasing yucca (-0.38), increasing sagebrush (0.73),
     increasing sandplum (0.59), decreasing cholla (-0.38), and
     decreasing other shrub sp (-0.42)")
