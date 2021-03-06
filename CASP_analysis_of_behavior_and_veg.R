#Cassin's Sparrows

###Load libraries
library(lme4)
library(lmerTest)
#library(raster)
#library(rasterVis)
library(tidyr)
library(tidyverse)

###Load data files
CASP_behavior <- read.csv(file = "20190417_CASP_playback_responses_combined_longform_corrected_BMSP_GPS.csv")
CASP_veg_daub <- read.csv(file = "combined_vegetation_percent_cover_CASP.csv")
CASP_veg_sp   <- read.csv(file = "combined_vegetation_surveys_CASP.csv")

###Taking raw data to summarized formats

#get everything under same column name to assign points
CASP_veg_daub$Pointnames <- CASP_veg_daub$Point
CASP_veg_sp$Pointnames <- CASP_veg_sp$Point

#Change $Response to be 0-2 as 0 (no response, could just be a neighbor) and 3-7 as 1 (yes).
CASP_behavior$Response[CASP_behavior$Strongest_behavior <= 2] <- 0
CASP_behavior$Response[CASP_behavior$Strongest_behavior > 2] <- 1

CASP_behavior <- tidyr::unite(CASP_behavior,
                       Pointnames,
                       Transect, Point, 
                       sep = "", 
                       remove = TRUE)

#Summarize daubenmire measurements

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
            sumOtherWoody = sum(Other.shrub + Tree,
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

#Create dataset that will be used in analyses.
behavior.veg <- left_join(behavior.veg1,
                          CASP_veg_sp_sum_spread,
                          by = c("Location",
                                 "Pointnames"))



#filter by points that have veg
behavior.veg <- behavior.veg %>%
  filter(!is.na(mSum))

#Add ratio of above/below.
behavior.veg$abratio <- behavior.veg$above1.sumCounts/behavior.veg$below1.sumCounts

###Reduce variables using PCA.

#First for Daubenmire.
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
               "above1.sumOtherWoody")

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
               "below1.sumOtherWoody")

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


check.for.cor.below <-cor(
  behavior.veg[,
               veg.below],
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

###ANALYSES

#Strongest reaction (as measured by distance of closest approach) by veg where present
lm.distance.veg <- lmer(ClosestDistance ~ daubPC1 + daubPC2 + daubPC3 +
                          abovePC1 + abovePC2 + abovePC3 +
                          belowPC1 + belowPC2 + belowPC3+(1|Location),
                        data = behavior.veg[behavior.veg$Response==1,])

summary(lm.distance.veg)

lm.distance.veg.no_RE <- lm(ClosestDistance ~ daubPC1 + daubPC2 + daubPC3 +
                          abovePC1 + abovePC2 + abovePC3 +
                          belowPC1 + belowPC2 + belowPC3,
                        data = behavior.veg[behavior.veg$Response==1,])

summary(lm.distance.veg.no_RE)

#Presence/absence of defense by veg

glm.presence.veg <- glmer(Response ~ daubPC1 + daubPC2 + daubPC3 +
                            abovePC1 + abovePC2 + abovePC3 +
                            belowPC1 + belowPC2 + belowPC3+(1|Location),
                          data = behavior.veg,
                          family = "binomial")

summary(glm.presence.veg)

glm.presence.veg.no_RE <- glm(Response ~ daubPC1 + daubPC2 + daubPC3 +
                            abovePC1 + abovePC2 + abovePC3 +
                            belowPC1 + belowPC2 + belowPC3,
                          data = behavior.veg,
                          family = "binomial")

summary(glm.presence.veg.no_RE)

#Table 3 (loadings for pc axes that were significant)
loadings.pca.above.df <- data.frame(loadings.pca.above)
loadings.pca.below.df <- data.frame(loadings.pca.below)

abovepc1 <- t(loadings.pca.above.df[1,])
abovepc3 <- t(loadings.pca.above.df[3,])
belowpc1 <- t(loadings.pca.below.df[1,])
loadings <- data.frame(cbind(abovepc1, belowpc1, abovepc3))

###FIGURES

#For manuscript, figures of the two significantly correlated PC axes.

# #Figure 2, abovePC1
# ndFig2<- data.frame("daubPC1" = mean(behavior.veg$daubPC1),
#                     "daubPC2" = mean(behavior.veg$daubPC2),
#                     "daubPC3" = mean(behavior.veg$daubPC3),
#                     "abovePC1"= seq(min(behavior.veg$abovePC1),
#                                     max(behavior.veg$abovePC1),
#                                     length.out=length(behavior.veg$abovePC1)),
#                     "abovePC2" = mean(behavior.veg$abovePC2),
#                     "abovePC3" = mean(behavior.veg$abovePC3),
#                     "belowPC1" = mean(behavior.veg$belowPC1),
#                     "belowPC2" = mean(behavior.veg$belowPC2),
#                     "belowPC3" = mean(behavior.veg$belowPC3))
# #plot the prediction with the new data (otherwise it uses rownumber and stretches the line out uselessly).
# 
# fig2_predict <- predict(glm.presence.veg.no_RE,
#         newdata=ndFig2,
#         type="response",
#         se.fit = TRUE)
# 
# svg("Fig2.svg",
#     width = 7,
#     height = 5)
# par(mar=c(7,5,5,4))
# plot(Response ~ abovePC1,
#      data = behavior.veg,
#      xlab = "",
#      ylab = "Presence of agonistic behavior in response to playback")
# mtext(
#   "Tall woody vegetation PC1: increasing sagebrush (0.65), increasing sandplum (0.59), 
#   decreasing cholla (-0.44), and decreasing other shrubs (-0.56)",
#   side=1, line=4)
# lines(ndFig2$abovePC1,
#       fig2_predict[[1]],
#       lty="solid",
#       lwd=2)
# lines(ndFig2$abovePC1,
#       fig2_predict[[1]] + fig2_predict[[2]],
#       lty="dotted",
#       lwd=2)
# lines(ndFig2$abovePC1,
#       fig2_predict[[1]] - fig2_predict[[2]],
#       lty="dotted",
#       lwd=2)
# dev.off()

#Figure 2, belowPC1
ndFig2<- data.frame("daubPC1" = mean(behavior.veg$daubPC1),
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
fig2_predict <- predict(glm.presence.veg.no_RE,
        newdata=ndFig2,
        type="response",
        se.fit = TRUE)

svg("Fig2.svg",
    width = 7,
    height = 5)
par(mar=c(7,5,5,4))
plot(Response ~ belowPC1,
     data = behavior.veg,
     xlab = "",
     ylab = "Presence of agonistic behavior in response to playback")
mtext(
  "Short woody vegetation PC1: decreasing yucca (-0.37), increasing sagebrush (0.73),
  increasing sandplum (0.60), decreasing cholla (-0.37), and
  decreasing other woody vegetation sp (-0.43)",
  side=1, line=5)
lines(x = ndFig2$belowPC1,
      y = fig2_predict[[1]],
      lty="solid",
      lwd=2)
lines(x = ndFig2$belowPC1,
      y = fig2_predict[[1]] + fig2_predict[[2]],
      lty="dotted",
      lwd=2)
lines(x = ndFig2$belowPC1,
      y = fig2_predict[[1]] - fig2_predict[[2]],
      lty="dotted",
      lwd=2)

dev.off()

# 
# #Figure 4, both PCA
# svg("Fig4.svg",
#     width = 7,
#     height = 7)
# par(mar=c(7,5,5,4))
# default.palette <- palette()
# palette(c("gray", "black"))
# plot(abovePC1 ~ belowPC1,
#      data = behavior.veg,
#      pch = as.numeric(as.factor(Location)),
#      col = Response+1,
#      bg = Response+1,
#      xlab = "",
#      ylab = "")
# mtext(
#   "BelowPC1: decreasing yucca (-0.38), increasing sagebrush (0.73),
#   increasing sandplum (0.59), decreasing cholla (-0.38), and
#   decreasing other shrub sp (-0.42)",
#   side=1, line=5)
# mtext(
#   "abovePC1: increasing sagebrush (0.65), increasing sandplum (0.59),
#   decreasing cholla (-0.44), and decreasing other shrubs (-0.56)",
#   side=2, line=2)
# legend("bottomright",
#        legend = c("Black Mesa State Park",
#                   "Cimarron Hills WMA",
#                   "Optima WMA Location 1",
#                   "Optima WMA Location 2",
#                   "Packsaddle WMA",
#                   "Rita Blanca WMA",
#                   "Selman Ranch",
#                   "Undefended point",
#                   "Defended point"),
#        pch = c(seq(1:7), 
#                1,
#                1),
#        col = c(rep(2, 7),
#                1,
#                2),
#        cex = 0.7)
# dev.off()
# 

###Map of study sites (Figure 1)
# 
# #Load ecoregions raster
# #create temporary raster files on large drive because they occupy 10-30 GB
# rasterOptions()$tmpdir
# rasterOptions(tmpdir=paste0(getwd(),
#                             "/CASP/rastertemp"))
# 
# 
# ecoregions <- raster(x = paste0(getwd(),
#                                 "/CASP/Raster/ok_vegetation.img"))
# 
# #check CRS
# crs(ecoregions)
# 
# #Add spatial data to points and transform (quicker than transforming whole raster)
# behavior.veg$lon <- behavior.veg$Longitude
# behavior.veg$lat <- behavior.veg$Latitude
# behavior.veg.sp <- behavior.veg
# coordinates(behavior.veg.sp) <- c("lat", 
#                                   "lon") #They are named backwards... lat is actually longitude and vice versa.
# proj4string(behavior.veg.sp) <- CRS("+init=epsg:4326 +proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs +towgs84=0,0,0")
# summary(behavior.veg.sp)
# behavior.veg.ecoregion <- spTransform(behavior.veg.sp,
#                                       crs(ecoregions))
# crs(behavior.veg.ecoregion)
# 
# #Assign metadata to raster
# library(XML)
# ecoregions.test <- xmlTreeParse(paste0(getwd(),
#                                        "/CASP/Raster/oklahoma_vegetation_raster_metadata.xml"))
# xml_data <- xmlToList(ecoregions.test)
# str(xml_data)
# types <- unlist(xml_data[[4]][[3]][[4]])
# types.m <-data.frame(matrix(data = types,
#                             ncol = 3,
#                             byrow = TRUE),
#                      stringsAsFactors = FALSE)
# colnames(types.m) <- c("ID",
#                        "regionname",
#                        "delete")
# types.m$ID <- as.numeric(types.m$ID)
# types.m$delete <- NULL
# 
# 
# #Extract ecoregion data from raster to study points.
# 
# behavior.veg$study_region_values <- raster::extract(ecoregions,
#                                                     behavior.veg.ecoregion)
# #Summarize using group_by to see what ecoregion each site is in
# #and if each site has more than one ecoregion.
# 
# #Table 1
# ecoregion.summary.sites <- behavior.veg %>% 
#   group_by(Location,
#            study_region_values,
#            Response) %>%
#   summarize("points" = n())%>%
#   left_join(.,
#             types.m,
#             by = c("study_region_values"="ID"))%>%
#   dplyr::select(Location, regionname, Response, points) %>%
#   arrange(Location, regionname, Response) %>%
#   print()

#Sample sizes given in results
sample.sizes <- behavior.veg %>% 
  group_by(Response) %>%
  summarize("points" = n())%>%
  print()


# #Figure 1
# 
# map <- extent(behavior.veg.ecoregion)
# small.eco <- crop(ecoregions,
#                   map+1000)
# small.eco.f <- as.factor(small.eco)
# ## Add a landcover column to the Raster Attribute Table
# rat <- levels(small.eco.f)[[1]]
# rat2 <- left_join(rat,
#                   types.m,
#                   by = c("ID"="ID"))
# rat[["landcover"]] <- rat2$regionname
# levels(small.eco.f) <- rat
# rat2$my_col <- rev(terrain.colors(n = nrow(rat2)))
# 
# my_habitats <- left_join(ecoregion.summary.sites,
#                                  rat2,
#                                  by = c("study_region_values"="ID")) %>%
#   distinct(regionname.x, my_col)%>%
#   print()
# ## Plot
# # levelplot(small.eco.f, 
# #           col.regions=rev(terrain.colors(nrow(rat2))),
# #           xlab="",
# #           ylab="")
# 
# 
# 
# plot(small.eco.f,
#      legend = FALSE, 
#      col = rat2$my_col)
# legend(x='bottomleft', 
#        legend = my_habitats$regionname.x, 
#        fill = my_habitats$my_col,
#        cex = 0.6,
#        ncol = 2)
# plot(behavior.veg.ecoregion, 
#      col = "black",
#      add = TRUE)


#New Fig. 1 without raster but with site labels
# https://gis.stackexchange.com/questions/222799/create-an-inset-map-in-r 
tidyverse <- "package:tidyverse"
detach(tidyverse, unload=TRUE, character.only = TRUE) #interferes with something in maps

library(maps)
library(GISTools)  


dev.off() #reset par
svg(file = "Fig1.svg",
    height = 5, 
    width = 8)

maps::map(database = 'state',
          regions = c('oklahoma'),
    fill = TRUE,
    col = "gray",
    mar=c(5,8,4,2)+0.1) # https://stackoverflow.com/questions/44806661/preventing-the-y-axis-label-being-chopped-off-by-r-maps
map.axes(cex.axis=0.8)
title(x = "Longitude", 
      y = "Latitude")
points(y = behavior.veg$Longitude,
       x = behavior.veg$Latitude,
       pch = as.numeric(as.factor(behavior.veg$Location)))

maps::map.scale(x=-102.75, y=34.2, 
                ratio=FALSE,
                relwidth=0.2)
north.arrow(xb=-101.5,
            yb=35,
            len=0.1,
            lab="N") 

# Inmap
par(usr=c(-300, -63, 22, 144))
rect(xleft =-126.2,
     ybottom = 23.8,
     xright = -65.5,
     ytop = 50.6,
     col = "white")
maps::map("usa", 
    xlim=c(-126.2,-65.5),
    ylim=c(23.8,50.6),add=T)
maps::map("state", 
    xlim=c(-126.2,-65.5),
    ylim=c(23.8,50.6),add=T,
    boundary = F,
    interior = T,
    lty=1)
maps::map("state", 
    region="oklahoma",
    fill=T, 
    add=T,
    col = "gray")
legend("topright",
       legend = c("Black Mesa State Park",
                  "Cimarron Hills WMA",
                  "Optima WMA Location 1",
                  "Optima WMA Location 2",
                  "Packsaddle WMA",
                  "Rita Blanca WMA",
                  "Selman Ranch"),
       pch = c(seq(1:7), 
               1,
               1),
       cex = 0.7)
dev.off()
