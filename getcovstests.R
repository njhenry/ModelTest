obsdates <- c(rep(as.Date("02/12/2014", format="%d/%m/%Y"), 2),
         rep(as.Date("12/03/2015", format="%d/%m/%Y"), 2)
)
locs <- cbind(runif(4, -15, -12.5), runif(4, 30, 32))
startpaths <- list("Z:/mastergrids/MODIS_Global/MCD43B4_BRDF_Reflectance/TCB/1km/8-Daily/",
                   "Z:/mastergrids/MODIS_Global/MCD43B4_BRDF_Reflectance/EVI/1km/8-Daily/")
endpaths <- list("_TCB_Filled_Data.tif",
                 "_EVI_Filled_Data.tif")
namesSupplied <- FALSE

outVals <- list()
#for each covariate
for(i in 1:length(startpaths)){
  startpath <- startpaths[[i]]
  endpath <- endpaths[[i]]
  outCovVals <- rep(NA, dim(locs)[1])
  #get list of times available
  files.l <- list.files(startpath, pattern=paste0(endpath, "$"))
  files.dat <- regmatches(files.l, regexpr(pattern="[0-9]{7}", files.l))
  files.dates <- as.Date(files.dat, "%Y%j")
  files.date.df <- data.frame(str=files.dat, dat=files.dates, fl=files.l)
  files.date.df <- files.date.df[order(files.date.df$dat), ]
  
  
  #for each date, find which date it is next to
  nearestTime <- c()
  for(i3 in 1:length(obsdates)){
    obsdate <- obsdates[i3]
    nearestindex <- max(which(files.date.df$dat < obsdate))
    ratio <- as.numeric(obsdate - files.date.df$dat[nearestindex]) / 
      as.numeric(files.date.df$dat[nearestindex + 1] - files.date.df$dat[nearestindex]) 
    nearestTime <- rbind(nearestTime, c(nearestindex, ratio))
  }
  
  #get unique indices and which dates are nearest to them
  unqNearestTimes <- unique(nearestTime[, 1])
  timesIndex <- list()
  for(i2 in 1:length(unqNearestTimes)){
    timesIndex[[i2]] <- which(nearestTime[, 1] == unqNearestTimes[[i2]])
  }
  
  #for each unique time open rasters and get relevant values
  for(j in 1:length(timesIndex)){
    # print(j)
    timeIndex <- timesIndex[[j]]
    unqNearestTime <- unqNearestTimes[j]
    # print(timeIndex)
    
    locsT <- locs[timeIndex, ]
    latmin <- floor(min(locsT[,1]))
    latmax <- ceiling(max(locsT[,1]))
    lonmin <- floor(min(locsT[,2]))
    lonmax <- ceiling(max(locsT[,2]))
    
    covpath <- paste0(startpath, grep(files.date.df$str[unqNearestTime], files.l, value=TRUE))
    # print(files.date.df$str[unqNearestTime])
    # print(grep(paste0("*", files.date.df$str[unqNearestTime], "*"), files.l, value=TRUE))
    # print(covpath)
    covfull <- raster(covpath)
    cov <- crop(covfull, extent(c(lonmin, lonmax, latmin, latmax)))
    if(length(timeIndex) == 1){
      locsT <- rbind(locs[timeIndex, ], locs[timeIndex, ])
      covVals <- extract(cov, SpatialPoints(locsT[,2:1]))[1]
      # outCovVals[timeIndex[[j]]] <- covVals
    }else{
      covVals <- extract(cov, SpatialPoints(locsT[,2:1]))
      # outCovVals[timeIndex[[j]]] <- covVals
    }
    
    covpathnext <- paste0(startpath, grep(files.date.df$str[unqNearestTime + 1], files.l, value=TRUE))
    # print(files.date.df$str[unqNearestTime])
    # print(grep(paste0("*", files.date.df$str[unqNearestTime], "*"), files.l, value=TRUE))
    # print(covpath)
    covfullnext <- raster(covpathnext)
    covnext <- crop(covfullnext, extent(c(lonmin, lonmax, latmin, latmax)))
    if(length(timeIndex) == 1){
      locsT <- rbind(locs[timeIndex, ], locs[timeIndex, ])
      covValsNext <- extract(covnext, SpatialPoints(locsT[,2:1]))[1]
      # outCovVals[timeIndex[[j]]] <- covVals
    }else{
      covValsNext <- extract(covnext, SpatialPoints(locsT[,2:1]))
      # outCovVals[timeIndex[[j]]] <- covVals
    }
    
    
    for(l in 1:length(timeIndex)){
      k <- timeIndex[l]
      # print(k)
      # print(outCovVals[k])
      # print(covVals * nearestTime[k, 2] + covValsNext * (1 - nearestTime[k, 2]))
      outCovVals[k] <- covVals[l] * nearestTime[k, 2] + covValsNext[l] * (1 - nearestTime[k, 2])
      # print(outCovVals)
    }
    
  }
  #print(i)
  if(namesSupplied) outVals[[covnames[[i]]]] <- outCovVals else outVals[[i]] <- outCovVals
  
}