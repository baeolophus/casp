#Cassin's Sparrows

library(tidyverse)

CASP_behavior <- read.csv(file = "CASP/20170207_CASP_playback_responses_combined_longform.csv")
CASP_veg_daub <- read.csv(file = "CASP/combined_vegetation_percent_cover_CASP.csv")
CASP_veg_sp   <- read.csv(file = "CASP/combined_vegetation_surveys_CASP.csv")

#get everything under same column name to assign points
CASP_veg_daub$Pointnames <- CASP_veg_daub$Point
CASP_veg_sp$Pointnames <- CASP_veg_sp$Point

#Change $Response to be 0-2 as 0 (no response, could just be a neighbor) and 3-7 as 1 (yes).
CASP_behavior$Response[CASP_behavior$Strongest_behavior <= 2] <- 0
CASP_behavior$Response[CASP_behavior$Strongest_behavior > 2] <- 1

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


check.for.cor.below <-cor(
  behavior.veg[,
               veg.below],
  method="pearson")

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

check.for.cor.above <-cor(
                         behavior.veg[,
                                      veg.above],
                         method="pearson")

#Check for above/below cors
check.for.cor.both <-cor(  behavior.veg[,
                                        veg.below],
  behavior.veg[,
               veg.above],
  method="pearson")
#Yucca above/below correlated at 0.89194720
#Cholla above/below correlated at 0.78313087
#Sage above/below correlated at 0.55894587
#However, none are exact corerlations and all provide information that we want,
#so I don't think we are
#overemphasizing their contribution here.  Everything else is below |0.5|.
#Similarly low correlations among types of vegetation in above and below categories on their own.

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
                          belowPC1 + belowPC2 + belowPC3+abratio+(1|Location),
                      data = behavior.veg[behavior.veg$Response==1,])

summary(lm.distance.veg)

lm.distance.abovebelow <- lmer(ClosestDistance ~ above1.sumCounts +below1.sumCounts+(1|Location),
                        data = behavior.veg[behavior.veg$Response==1,])

summary(lm.distance.abovebelow)

#Presence/absence by veg

glm.presence.veg <- glmer(Response ~ daubPC1 + daubPC2 + daubPC3 +
                          abovePC1 + abovePC2 + abovePC3 +
                          belowPC1 + belowPC2 + belowPC3+ abratio+(1|Location),
                      data = behavior.veg,
                      family = "binomial")

summary(glm.presence.veg)

#Simply using counts alone wastes a lot of info.
glm.presence.abovebelow <- glmer(Response ~ above1.sumCounts +below1.sumCounts+abratio+ (1|Location),
                          data = behavior.veg,
                          family = "binomial")

summary(glm.presence.abovebelow)

#Add ratio of above/below.

behavior.veg$abratio <- behavior.veg$above1.sumCounts/behavior.veg$below1.sumCounts

#Per behavior tests vs absence
library(mixcat)
attach(behavior.veg)
behavior.veg$fLocation<- as.factor(behavior.veg$Location)
test <-npmlt(Strongest_behavior+1 ~ daubPC1 + daubPC2 + daubPC3 +
        abovePC1 + abovePC2 + abovePC3 +
        belowPC1 + belowPC2 + belowPC3,

      random=~1,
      id=fLocation,
      k=2,
      EB=FALSE)

summary(test)
z <- test$coefficients/test$SE.coefficients
z
##          (Intercept) sesmiddle seshigh  write
## general        2.445   -1.2018  -2.261 -2.706
## vocation       4.485    0.6117  -1.650 -5.113
# 2-tailed z test
p <- (1 - pnorm(abs(z), 0, 1)) * 2
p
round(p,2)
#possibility of surface map-
#CASP presence vs vegtation to account for spatial autocorrelation (some birds are same)


#Figures

ndFig2<- data.frame("daubPC1" = mean(behavior.veg$daubPC1),
                    "daubPC2" = mean(behavior.veg$daubPC2),
                    "daubPC3" = mean(behavior.veg$daubPC3),
                    "abovePC3"=seq(min(behavior.veg$abovePC3),
                                max(behavior.veg$abovePC3),
                                length.out=length(behavior.veg$abovePC3)),
                    "abovePC2" = mean(behavior.veg$abovePC2),
                    "abovePC1" = mean(behavior.veg$abovePC1),
                    "belowPC1" = mean(behavior.veg$belowPC1),
                    "belowPC2" = mean(behavior.veg$belowPC2),
                    "belowPC3" = mean(behavior.veg$belowPC3))
#plot the prediction with the new data (otherwise it uses rownumber and stretches the line out uselessly).


lines(ndFig2$abovePC3,
      predict(glm.presence.veg,
                  newdata=ndFig2,
                  type="response",
              re.form = NA),
      lty="solid",lwd=2)

plot(Response ~ abovePC3,
     data = behavior.veg,
     xlab = "AbovePC3: decreasing yucca (-0.76), increasing trees (0.53)")


ndFig3<- data.frame("daubPC1" = mean(behavior.veg$daubPC1),
                    "daubPC2" = mean(behavior.veg$daubPC2),
                    "daubPC3" = mean(behavior.veg$daubPC3),
                    "abovePC1"= mean(behavior.veg$abovePC1),
                    "abovePC2" = mean(behavior.veg$abovePC2),
                    "abovePC3" = mean(behavior.veg$abovePC3),
                    "belowPC1" = seq(min(behavior.veg$belowPC1),
                                     max(behavior.veg$belowPC1),
                                     length.out=length(behavior.veg$belowPC1)),
                    "belowPC2" = mean(behavior.veg$belowPC2),
                    "belowPC3" = mean(behavior.veg$belowPC3))
#plot the prediction with the new data (otherwise it uses rownumber and stretches the line out uselessly).



par(mar=c(7,5,5,4))
plot(Response ~ belowPC1,
     data = behavior.veg,
     xlab = "")
mtext(
      "BelowPC1: decreasing yucca (-0.38), increasing sagebrush (0.73),
     increasing sandplum (0.59), decreasing cholla (-0.38), and
     decreasing other shrub sp (-0.42)",
      side=1, line=5)
lines(ndFig3$belowPC1,
      predict(glm.presence.veg,
              newdata=ndFig2,
              type="response",
              re.form = NA),
      lty="solid",lwd=2)

#Interpretation for marginal PC


ndnone<- data.frame("daubPC1" = mean(behavior.veg$daubPC1),
                    "daubPC2" = mean(behavior.veg$daubPC2),
                    "daubPC3" = mean(behavior.veg$daubPC3),
                    "abovePC1"= seq(min(behavior.veg$abovePC1),
                                    max(behavior.veg$abovePC1),
                                    length.out=length(behavior.veg$abovePC1)),
                    "abovePC2" = mean(behavior.veg$abovePC2),
                    "abovePC3" = mean(behavior.veg$abovePC3),
                    "belowPC1" = mean(behavior.veg$belowPC1),
                    "belowPC2" = mean(behavior.veg$belowPC2),
                    "belowPC3" = mean(behavior.veg$belowPC3))
#plot the prediction with the new data (otherwise it uses rownumber and stretches the line out uselessly).



par(mar=c(7,5,5,4))
plot(Response ~ abovePC1,
     data = behavior.veg,
     xlab = "")
mtext(
  "abovePC1: increasing sagebrush (0.65), increasing sandplum (0.59), 
  decreasing cholla (-0.44), and decreasing other shrubs (-0.56)",
  side=1, line=4)
lines(ndnone$abovePC1,
      predict(glm.presence.veg,
              newdata=ndnone,
              type="response",
              re.form = NA),
      lty="solid",lwd=2)

par(mar=c(7,5,5,4))
plot(abovePC3~
     belowPC1,
     data = behavior.veg,
     pch = as.numeric(as.factor(Location)),
     bg = as.factor(Response),
     col = as.factor(Response),
     xlab = "",
     ylab = "AbovePC3: decreasing yucca (-0.76), increasing trees (0.53)")
mtext(
  "BelowPC1: decreasing yucca (-0.38), increasing sagebrush (0.73),
     increasing sandplum (0.59), decreasing cholla (-0.38), and
  decreasing other shrub sp (-0.42)",
  side=1, line=5)
