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
behavior.veg <- left_join(behavior.veg1,
                          CASP_veg_sp,
                          by = c("Location",
                                 "Pointnames"))

#Strongest reaction by veg where present
lm.presence.veg <- lm(ClosestDistance ~ Yucca +Sage + Sandplum + Cholla + Tree + Other.shrub,
                      data = behavior.veg[behavior.veg$Response==1,])

summary(lm.presence.veg)

#Presence/absence by veg

lm.presence.veg <- lm(c(Yucca +Sage + Sandplum + Cholla + Tree + Other.shrub) ~ Response,
                      data = behavior.veg[behavior.veg$Height==">1m",])

summary(lm.presence.veg)
