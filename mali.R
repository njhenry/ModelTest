setwd("C:/Users/scro3122/Documents/ModelTest")
library(dplyr)
library(ggplot2)
library(lubridate)
inc <- read.csv("../indi.data.for.treatmentseeking.final.csv")

years <- dplyr::select(inc, year)[[1]]
hist(years)
print(table(years))

inc2015 <- filter(inc, year == 2015)
country2015 <- dplyr::select(inc2015, Country)
print(table(country2015))

print(sum(dplyr::select(filter(inc2015, Country == "Nigeria"), rdt.result)[[1]]))
print(sum(dplyr::select(filter(inc2015, Country == "Mali"), rdt.result)[[1]]))
print(sum(dplyr::select(filter(inc2015, Country == "Rwanda"), rdt.result)[[1]]))
print(sum(dplyr::select(filter(inc2015, Country == "Kenya"), rdt.result)[[1]]))
print(table(dplyr::select(filter(inc2015, Country == "Mali"), urb.rural)))


##do Mali in 2015
mali15 <- filter(inc2015, Country == "Mali")
locs <- cbind(mali15$LONG, mali15$LAT)
plot(locs)

########################
# #to be the same as botswana covs
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
# obsdates <- as.Date(paste0("1/", mali15$month, "/", mali15$year), format="%d/%m/%Y")
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
#save(list = c("staticCovs", "tempCovs", "startpathsT", "endpathsT", "covpaths"), file = "MaliCovs_sameasBWA.RData")
load("MaliCovs_sameasBWA.RData")
######################


# source("getCovs.R")
# source("C:/Users/scro3122/Documents/ModelTest/getTempCovs3Par.R")
# #temporal covs
# startpathsT <- list("Z:/mastergrids/MODIS_Global/MOD11A2_LST/LST_Day/1km/Monthly/",
#                     "Z:/mastergrids/Other_Global_Covariates/TemperatureSuitability/TSI_Pf_Dynamic/1km/Monthly/",
#                     "Z:/mastergrids/Other_Global_Covariates/Rainfall/CHIRPS/1k/Monthly/")
# 
# 
# endpathsT <- list(".Mean.1km.Data.tif",
#                   ".Data.1km.Data.tif",
#                   ".1km.tif")
# 
# obsdates <- as.Date(paste0("1/", mali15$month, "/", mali15$year), format="%d/%m/%Y")
# byMonth=c(TRUE, TRUE, TRUE)
# 
# tempcovnames <- list("LST", "TempSuit", "rain")
# timelags <- list(c(2,4,6,8), c(2,4,6,8), c(2,4,6,8))
# 
# 
# 
# tempCovs <- getTempCovsDate3Par(locs[,2:1], obsdates, startpathsT, endpathsT, tempcovnames, byMonth,
#                                 timelags = list(c(2,4,6,8), c(2,4,6,8), c(2,4,6,8), c(2,4,6,8)))

##old temporal covs and current static covs I guess
# 
# # 
# # 
# # tempCovs <- getTempCovsMonth(locs[,2:1], mali15$year, mali15$month, startpathsT, endpathsT, tempcovnames, timelags)
# # 
# # 
# # #static covs
# # covpaths <- list("Z:/mastergrids/Other_Global_Covariates/UrbanAreas/Global_Urban_Footprint/From_12m/1km/Global_Urban_Footprint_1km_12m-Prop-Urban-Unclipped.tif",
# #                  "Z:/mastergrids/Other_Global_Covariates/Population/Worldpop_GPWv4_Hybrid_201708/1km/Global_Hybrid_Pop_v2_1km_UNAdj_2017-Interp.tif",
# #                  "Z:/mastergrids/Other_Global_Covariates/NightTimeLights/DMSP_F18_nighttime_lights_2010_1km_global.tif",
# #                  "Z:/mastergrids/Other_Global_Covariates/Elevation/Ferranti-Elevation/1km/Synoptic/ferranti_30sec_elev_max_clip.tif",
# #                  "Z:/mastergrids/Other_Global_Covariates/Accessibility/Weiss/friction_surface_2015_v1.tif")
# # staticCovs <- getStaticCovs(locs[,2:1], covpaths, covnames=list("Urban", "FrankenPop", "NightLights", "Elevation", "Friction"))
# # 

# save(list = c("staticCovs", "tempCovs", "startpathsT", "endpathsT", "covpaths"), file = "MaliCovs.RData")
# load("MaliCovs.RData")




# df <- data.frame(x=locs[,1], y=locs[,2], LST=tempCovs$LST, month=as.numeric(mali15$month) < 10, TempSuit=tempCovs$TempSuit,
#                  urban=staticCovs$Urban)
# p <- ggplot(df, aes(x=x,y=y,color=LST)) + geom_point()
# print(p)
# p <- ggplot(df, aes(x=x,y=y,color=TempSuit)) + geom_point()
# print(p)
# p <- ggplot(df, aes(x=x,y=y,color=month)) + geom_point()
# print(p)
# p <- ggplot(df, aes(x=x,y=y,color=urban)) + geom_point()
# print(p)
# fullpath <- "Z:/mastergrids/MODIS_Global/MOD11A2_LST/LST_Day/1km/Monthly/LST_Day.2015.09.Mean.1km.Data.tif"
# r <- crop(raster(fullpath), extent(c(-12.5, -2.5, 10, 16)))
# r2 <- crop(raster(covpaths[[1]]), extent(c(-12.5, -2.5, 10, 16)))
# plot(r)
# sum(is.na(values(r)))
# 
# library(FSIC)
# pr <- gpReg(c(locs, mali15$month), mali15$rdt.result)
# temp <- gpReg(c(locs, mali15$month), tempCovs$LST)
# ts <- gpReg(c(locs, mali15$month), tempCovs$TempSuit)
# 
# plot(temp, pr)
# plot(ts, pr)
# 
# df <- data.frame(x=locs[,1], y=locs[,2], smoothtemp = temp + tempCovs$LST)
# p <- ggplot(df, aes(x=x,y=y, color=smoothtemp)) + geom_point()
# print(p)

obs <- c(staticCovs, tempCovs)
obs$"net" <- mali15$net
obs$"net"[is.na(obs$"net")] <- 0
obs$"pr" <- mali15$rdt.result

#save(list=c("obs", "mali15", "locs"), file = "MaliObs.RData")
setwd("C:/Users/scro3122/Documents/ModelTest")
source("runpc.R")
setwd("C:/Users/scro3122/Documents/ModelTest")

# 
##set graph
nvar <- length(obs)
G_0 <- matrix(1, nrow=nvar, ncol=nvar) - diag(nvar)
G_0[length(obs), ] <- 0
temp.index <- c(6, 11, 16)
temp.index.m2 <- temp.index + 1
temp.index.m4 <- temp.index.m2 + 1
temp.index.m6 <- temp.index.m4 + 1
temp.index.m8 <- temp.index.m6 + 1

G_0[temp.index, c(temp.index.m2, temp.index.m4, temp.index.m6, temp.index.m8)] <- 0
G_0[temp.index.m2, c(temp.index.m4, temp.index.m6, temp.index.m8)] <- 0
G_0[temp.index.m4, c(temp.index.m6, temp.index.m8)] <- 0
G_0[temp.index.m6, temp.index.m8] <- 0


#save(list=c("obs", "mali15", "locs", "G_0"), file = "Mali_pc_dat_same_BWA.RData")
ptm <- proc.time()
pc <- runpc(obs, locs, mali15$month, 0.1, nSample=250, G_0=G_0)
print(proc.time() - ptm)

plot.minimal(pc[[1]], targetIndex = length(obs), names(obs))

# #######aggregate
# lat.min <- min(locs[,1])
# lat.max <- max(locs[,1])
# lon.min <- min(locs[,2])
# lon.max <- max(locs[,2])
# 
# N.blocks <- 25
# obsAgg <- list()
# point.block <- c()
# 
# #work out which points are in which block
# for(i in 1:(dim(locs)[1])){
#   loc <- locs[i, ]
#   lat <- loc[1]
#   lon <- loc[2]
#   
#   lat.n <- round(N.blocks * (lat - lat.min) / (lat.max - lat.min))
#   lon.n <- round(N.blocks * (lon - lon.min) / (lon.max - lon.min))
#   point.block[i] <- lat.n * N.blocks + lon.n + 1
# }
# 
# centres <- c()
# blocks <- unique(point.block)
# #agregate for each block for each variable
# for(i in 1:length(obs)){
#   var.agg <- c()
#   for(j in 1:length(blocks)){
#     block <- blocks[j]
#     var.agg <- c(var.agg, mean(obs[[i]][which(point.block==block)]))
#   }
#   obsAgg[[i]] <- var.agg
# }
# 
# for(j in 1:length(blocks)){
#   block <- blocks[j]
#   centres <- rbind(centres, colMeans(locs[which(point.block==block), ]))
# }
# 
# names(obsAgg) <- names(obs)
# 
# ptm <- proc.time()
# pc <- runpc(obsAgg, locs, mali15$month, 0.2)
# print(proc.time() - ptm)






# 
# 
# ##########ICP
# source("C:/Users/scro3122/Documents/ModelTest/ICP.R")
# obs$'cluster' <- mali15$cluster.x

# 
# icp.test <- icp(obs, locs, mali15$month, 0.5, 13, 14, nSample=600)
# 
# #(obs, locs, times , alpha, target.index, environ.index, nSample = length(obsDat[[1]]))
# print(names(obs)[unique(unlist(icp.test[[3]]))])
  