setwd("C:/Users/scro3122/Documents/ModelTest")
library(dplyr)
library(ggplot2)
library(lubridate)
inc <- read.csv("../indi.data.for.treatmentseeking.final.csv")

years <- dplyr::select(inc, year)[[1]]
hist(years)
print(table(years))

inc2015 <- filter(inc, year == 2015)

##do Mali in 2015
ngr15 <- filter(inc2015, Country == "Nigeria")
locs <- cbind(ngr15$LONG, ngr15$LAT)
plot(locs)

# 
# #covaraites
# source("getCovs.R")
# source("C:/Users/scro3122/Documents/ModelTest/getTempCovs3Par.R")
# #temporal covs
# startpathsT <- list("Z:/mastergrids/MODIS_Global/MOD11A2_v6_LST/LST_Day/1km/8-Daily/",
#                     "Z:/mastergrids/MODIS_Global/MCD43B4_BRDF_Reflectance/EVI/1km/8-Daily/",
#                     "Z:/mastergrids/Other_Global_Covariates/Rainfall/CHIRPS/1k/Monthly/")
# 
# 
# endpathsT <- list(".Data.1km.Data.tif",
#                   "_EVI_Filled_Data.tif",
#                   ".1km.tif")
# 
# obsdates <- as.Date(paste0("1/", ngr15$month, "/", ngr15$year), format="%d/%m/%Y")
# byMonth=c(FALSE, FALSE, TRUE)
# 
# tempcovnames <- list("LST", "EVI", "rain")
# timelags <- list(c(2,4,6,8), c(2,4,6,8), c(2,4,6,8))
# 
# 
# 
# tempCovs <- getTempCovsDate3Par(locs[,2:1], obsdates, startpathsT, endpathsT, tempcovnames, byMonth,
#                                 timelags = list(c(2,4,6,8), c(2,4,6,8), c(2,4,6,8), c(2,4,6,8)))
# #static covs
# covpaths <- list("Z:/mastergrids/Other_Global_Covariates/UrbanAreas/Global_Urban_Footprint/From_12m/1km/Global_Urban_Footprint_1km_12m-Prop-Urban-Unclipped.tif",
#                  "Z:/mastergrids/Other_Global_Covariates/Population/Worldpop_GPWv4_Hybrid_201708/1km/Global_Hybrid_Pop_v2_1km_UNAdj_2017-Interp.tif",
#                  "Z:/mastergrids/Other_Global_Covariates/NightTimeLights/DMSP_F18_nighttime_lights_2010_1km_global.tif",
#                  "Z:/mastergrids/Other_Global_Covariates/Elevation/Ferranti-Elevation/1km/Synoptic/ferranti_30sec_elev_max_clip.tif",
#                  "Z:/mastergrids/Other_Global_Covariates/Accessibility/Weiss/friction_surface_2015_v1.tif")
# staticCovs <- getStaticCovs(locs[,2:1], covpaths, covnames=list("Urban", "FrankenPop", "NightLights", "Elevation", "Friction"))
# 
# save(list = c("staticCovs", "tempCovs", "startpathsT", "endpathsT", "covpaths"), file = "Nigeria_covs.RData")

load("Nigeria_covs.RData")


###
obs.raw <- c(tempCovs, staticCovs)
obs.raw$"pr" <- ngr15$rdt.result
obs.raw$"net" <- ngr15$net

edu <- ngr15$mother.education
edu.use <- rep(0, length(ngr15$mother.education))
edu.use[edu == "No education"] <- 1
edu.use[edu == "Primary"] <- 2
edu.use[edu == "Secondary"] <- 3
edu.use[edu == "Higher"] <- 5

obs.raw$"edu" <- edu.use
obs.raw$"wealth" <- ngr15$wealth.index.modelling
obs <- as.list(na.omit(as.data.frame(obs.raw)))


nvar <- length(obs)
G_0 <- matrix(1, nrow=nvar, ncol=nvar) - diag(nvar)

temp.index <- c(1,6,11)
temp.index.m2 <- temp.index + 1
temp.index.m4 <- temp.index.m2 + 1
temp.index.m6 <- temp.index.m4 + 1
temp.index.m8 <- temp.index.m6 + 1

G_0[temp.index, c(temp.index.m2, temp.index.m4, temp.index.m6, temp.index.m8)] <- 0
G_0[temp.index.m2, c(temp.index.m4, temp.index.m6, temp.index.m8)] <- 0
G_0[temp.index.m4, c(temp.index.m6, temp.index.m8)] <- 0
G_0[temp.index.m6, temp.index.m8] <- 0

#malaria doesn't cause things
G_0[21, ] <- 0

#nets can only cause malaria
G_0[22, ] <- 0
G_0[22, 21] <- 1

#education can only cause nets or malaria
G_0[23, ] <- 0
G_0[23, 21:22] <- 1

#wealth can only cause education, nets and malaria
G_0[24, ] <- 0
G_0[24, 21:23] <- 0


save(list=c("obs", "ngr15", "locs", "G_0"), file = "Nigeria_dat_wwealth.RData")

