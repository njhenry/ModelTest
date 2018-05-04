# obsdates <- c(rep(as.Date("02/12/2014", format="%d/%m/%Y"), 2),
#               rep(as.Date("12/03/2015", format="%d/%m/%Y"), 2)
# )
# locs <- cbind(runif(4, -15, -12.5), runif(4, 30, 32))
# folderpaths <- list("Z:/mastergrids/MODIS_Global/MOD11A2_LST/LST_Day/1km/Monthly/",
#                     "Z:/mastergrids/MODIS_Global/MCD43B4_BRDF_Reflectance/TCB/1km/8-Daily/",
#                    "Z:/mastergrids/MODIS_Global/MCD43B4_BRDF_Reflectance/EVI/1km/8-Daily/")
# endpaths <- list(".Mean.1km.Data.tif",
#                  "_TCB_Filled_Data.tif",
#                  "_EVI_Filled_Data.tif")
# covnames <- NULL
# byMonth <- c(TRUE, FALSE, FALSE)

library(raster)
library(ggplot2)
library(lubridate)

getTempCovsDate2 <- function(locs, obsdates, folderpaths, endpaths, covnames=NULL, byMonth=rep(TRUE, length(folderpaths))){
  outVals <- list()
  if(is.null(covnames)) namesSupplied <- FALSE else namesSupplied <- TRUE
  
  ##for each covariate]
  for(i in 1:length(folderpaths)){
    if(namesSupplied) print(covnames[[i]]) else print(paste0("cov ", i))
    outCovVals <- rep(NA, length(obsdates))
    folderpath <- folderpaths[[i]]
    endpath <- endpaths[[i]]
    obsmonths <- as.numeric(month(obsdates))
    obsyears <- as.numeric(year(obsdates))
    
    #get list of dates available
    file.names <- list.files(folderpath, pattern=paste0(endpath, "$") )
    if(byMonth[i]){
      file.dates <- regmatches(file.names, regexpr("[0-9]{4}\\.[0-9]{2}", file.names))
      file.years <- as.numeric(substr(file.dates, 1, 4))
      file.months <- as.numeric(substr(file.dates, 6, 7))
      # print(file.months)
      
      nearestFile <- c()
      for(k in 1:length(obsdates)){
        nearestFile[k] <- which.min(abs(file.years*12+file.months - (obsyears[k]*12+obsmonths[k])))
      }
      # print(nearestFile)
      
    }
    else{
      file.dates <- as.Date(regmatches(file.names, regexpr("[0-9]{7}", file.names)), format="%Y%j")
      
      #sort file names by date in case not already in the right order
      file.names <- file.names[order(file.dates)]
      file.dates <- sort(file.dates)
      
      #for each obsdate find which of the files are the nearest
      nearestFile <- c()
      for(k in 1:length(obsdates)){
        nearestFile[k] <- which.min(abs(file.dates - obsdates[k]))
      }
    }
    
    nearestIndex <- sapply(1:length(file.names), function(k) which(nearestFile==k))
    
    #find number of files that need to be opened
    nonempty <- c()
    for(j in 1:length(nearestIndex)){
      if(length(nearestIndex[[j]] > 0)) nonempty <- c(nonempty, j)
    }
    #for each file that has a non-zero number of nearest obsdates, open that covariate and
    #extract data at locations
    for(l in 1:length(nonempty)){
      j <- nonempty[l]
      print(paste0(l, " of ", length(nonempty)))
      nearIndex <- nearestIndex[[j]]
      if(length(nearIndex) > 0){
        covpath <- paste0(folderpath, file.names[j])
        # print(covpath)
        cov <- raster(covpath)
        #extract values
        cov.vals <- extract(cov, SpatialPoints(locs[nearIndex, 2:1, drop=FALSE]))
        outCovVals[nearIndex] <- cov.vals
        rm(cov)
      }
    }

    if(namesSupplied) outVals[[covnames[[i]]]] <- outCovVals else outVals[[i]] <- outCovVals
  }
  return(outVals)
}



#v <- getTempCovsDate2(locs, obsdates, folderpaths, endpaths, byMonth=byMonth)