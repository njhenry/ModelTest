setwd("C:/Users/scro3122/Documents")
library(dplyr)
library(ggplot2)
inc <- read.csv("indi.data.for.treatmentseeking.final.csv")

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
source("getCovs.R")

#temporal covs
startpathsT <- list("Z:/mastergrids/MODIS_Global/MOD11A2_LST/LST_Day/1km/Monthly/LST_Day.",
                   "Z:/mastergrids/Other_Global_Covariates/TemperatureSuitability/TSI_Pf_Dynamic/1km/Monthly/TSI-Martens2-Pf."
                   )

endpathsT <- list(".Mean.1km.Data.tif",
                 ".Data.1km.Data.tif")

tempCovs <- getTempCovsMonth(locs[,2:1], mali15$year, mali15$month, startpathsT, endpathsT, list("LST", "TempSuit"))


#static covs
covpaths <- list("Z:/mastergrids/Other_Global_Covariates/UrbanAreas/Global_Urban_Footprint/From_12m/1km/Global_Urban_Footprint_1km_12m-Prop-Urban-Unclipped.tif")
staticCovs <- getStaticCovs(locs[,2:1], covpaths, covnames=list("Urban"))







df <- data.frame(x=locs[,1], y=locs[,2], LST=tempCovs$LST, month=as.numeric(mali15$month) < 10, TempSuit=tempCovs$TempSuit,
                 urban=staticCovs$Urban)
p <- ggplot(df, aes(x=x,y=y,color=LST)) + geom_point()
print(p)
p <- ggplot(df, aes(x=x,y=y,color=TempSuit)) + geom_point()
print(p)
p <- ggplot(df, aes(x=x,y=y,color=month)) + geom_point()
print(p)
p <- ggplot(df, aes(x=x,y=y,color=urban)) + geom_point()
print(p)
fullpath <- "Z:/mastergrids/MODIS_Global/MOD11A2_LST/LST_Day/1km/Monthly/LST_Day.2015.09.Mean.1km.Data.tif"
r <- crop(raster(fullpath), extent(c(-12.5, -2.5, 10, 16)))
r2 <- crop(raster(covpaths[[1]]), extent(c(-12.5, -2.5, 10, 16)))
plot(r)
sum(is.na(values(r)))

library(FSIC)
pr <- gpReg(c(locs, mali15$month), mali15$rdt.result)
temp <- gpReg(c(locs, mali15$month), tempCovs$LST)
ts <- gpReg(c(locs, mali15$month), tempCovs$TempSuit)

plot(temp, pr)
plot(ts, pr)

df <- data.frame(x=locs[,1], y=locs[,2], smoothtemp = temp + tempCovs$LST)
p <- ggplot(df, aes(x=x,y=y, color=smoothtemp)) + geom_point()
print(p)
