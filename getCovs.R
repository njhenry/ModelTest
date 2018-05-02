###functions to automatically get covariate values from mastergrids
library(raster)
library(ggplot2)
library(lubridate)
runtests <- FALSE

##function to get static covariates
getStaticCovs <- function(locs, covpaths, covnames=NULL){
  outVals <- list()
  if(is.null(covnames)) namesSupplied <- FALSE else namesSupplied <- TRUE
  
  #load each covariate
  i <- 1
  for(covpath in covpaths){
    if(namesSupplied) print(covnames[[i]]) else print(paste0("cov ", i))
    cov <- raster(covpath)
    covVals <- extract(cov, SpatialPoints(locs[,2:1]))
    
    #put result into output
    if(namesSupplied) outVals[[covnames[[i]]]] <- covVals else outVals[[i]] <- covVals
    
    i <- i + 1
  }
  return(outVals)
}

getTempCovsMonth <- function(locs, year, month, startpaths, endpaths, covnames=NULL,
                             timelags=rep(0, length(startpaths))){
  outVals <- list()
  if(is.null(covnames)) namesSupplied <- FALSE else namesSupplied <- TRUE
  
  ##get ranges of loc dealt with

  
  ##sort data into blocks from same month and year
  times <- cbind(year, month)
  unqTimes <- unique(times, MARGIN=1)
  
  #list to keep track of which indices are at what time
  timeIndex <- list()
  
  for(i in 1:dim(unqTimes)[1]){
    unqTime <- unqTimes[i, ]
    #make function that tests equality with unqTime
    eqUnqTime <- function(x){
      return(all(x == unqTime))
    }
    timeIndex[[i]] <- which(apply(times, 1, eqUnqTime))
  }
  
  ##open each covariate for each time and get values  
  for(i in 1:length(startpaths)){
    if(namesSupplied) print(covnames[[i]]) else print(paste0("cov ", i))
    startpath <- startpaths[[i]]
    endpath <- endpaths[[i]]
    outCovVals <- rep(NA, length(year))
    #open cov at each time and extract relevant values
    for(j in 1:dim(unqTimes)[1]){
      locsT <- locs[timeIndex[[j]], ]
      latmin <- floor(min(locsT[,1]))
      latmax <- ceiling(max(locsT[,1]))
      lonmin <- floor(min(locsT[,2]))
      lonmax <- ceiling(max(locsT[,2]))
      unqTime <- unqTimes[j, ] + c(0, timelags[i])

      while(unqTime[2] > 12){
        unqTime[2] <- unqTime[2] - 12
        unqTime[1] <- unqTime[1] + 1
      }

      if(as.numeric(unqTime[2]) < 10){
        covpath <- paste0(startpath, unqTime[1], ".0", as.numeric(unqTime[2]), endpath)
      }else{
        covpath <- paste0(startpath, unqTime[1], ".", unqTime[2], endpath)
      }
      covfull <- tryCatch(raster(covpath),
                          error = function(e){
                            
                            print(unqTime[2] < 10)
                            print(unqTime[2])
                            print(covpath)
                            print(e)
                          })
      cov <- crop(covfull, extent(c(lonmin, lonmax, latmin, latmax)))
      if(length(timeIndex[[j]]) == 1){
        locsT <- rbind(locs[timeIndex[[j]], ], locs[timeIndex[[j]], ])
        covVals <- extract(cov, SpatialPoints(locsT[,2:1]))[1]
        outCovVals[timeIndex[[j]]] <- covVals
      }else{
        
        covVals <- extract(cov, SpatialPoints(locsT[,2:1]))
        outCovVals[timeIndex[[j]]] <- covVals
      }
      

    }
    if(namesSupplied) outVals[[covnames[[i]]]] <- outCovVals else outVals[[i]] <- outCovVals
  }
  return(outVals)
}


getTempCovsDate <- function(locs, obsdates, startpaths, endpaths, covnames=NULL, 
                            byMonth=rep(TRUE, length(startpaths))){
  outVals <- list()
  if(is.null(covnames)) namesSupplied <- FALSE else namesSupplied <- TRUE
  
  obsMonth <- c()
  ##get times for month covariates
  ##keep track of year, start month and what fraction through the month it is
  for(i in 1:length(obsdates)){
    obsdate <- obsdates[i]
    # print(obsdate)
    obsMonth <- rbind(obsMonth, c(year(obsdate), month(obsdate), day(obsdate)/(days_in_month(obsdate)+1)))
  }
  # print(obsMonth)
  unqTimesMonth <- unique(obsMonth[, 1:2], MARGIN=1)
  
  #list to keep track of which indices are at what time
  timeIndexMonth <- list()
  
  for(i in 1:dim(unqTimesMonth)[1]){
    unqTimeMonth <- unqTimesMonth[i, ]
    #make function that tests equality with unqTime
    eqUnqTimeMonth <- function(x){
      return(all(x == unqTimeMonth))
    }
    timeIndexMonth[[i]] <- which(apply(obsMonth[, 1:2], 1, eqUnqTimeMonth))
  }
  
  
  ##open each covariate for each time and get values  
  for(i in 1:length(startpaths)){
    if(namesSupplied) print(covnames[[i]]) else print(paste0("cov ", i))
    startpath <- startpaths[[i]]
    endpath <- endpaths[[i]]
    outCovVals <- rep(NA, dim(locs)[1])
    
    ##if MONTH
    # print(byMonth[i])
    if(byMonth[i]){
      #open cov at each time and extract relevant values
      for(j in 1:dim(unqTimesMonth)[1]){
        print(paste0(j, " of ", dim(unqTimesMonth)[1]))
        # print("j")
        # print(j)
        locsT <- locs[timeIndexMonth[[j]], ]
        latmin <- floor(min(locsT[,1]))
        latmax <- ceiling(max(locsT[,1]))
        lonmin <- floor(min(locsT[,2]))
        lonmax <- ceiling(max(locsT[,2]))
        
        unqTimeMonth <- unqTimesMonth[j, ]
        # unqTime <- unqTimes[j, ] + c(0, timelags[i])
        # while(unqTime[2] > 12){
        #   unqTime[2] <- unqTime[2] - 12
        #   unqTime[1] <- unqTime[1] + 1
        # }
        
        #open current month and extract values
        if(as.numeric(unqTimeMonth[2]) < 10){
          covpath <- paste0(startpath, unqTimeMonth[1], ".0", as.numeric(unqTimeMonth[2]), endpath)
        }else{
          covpath <- paste0(startpath, unqTimeMonth[1], ".", unqTimeMonth[2], endpath)
        }
        covfull <- tryCatch(raster(covpath),
                            error = function(e){
                              print("test2")
                              print(unqTimeMonth[2] < 10)
                              print(unqTimeMonth)
                              print(covpath)
                              print(e)
                            })
        cov <- crop(covfull, extent(c(lonmin, lonmax, latmin, latmax)))
        if(length(timeIndexMonth[[j]]) == 1){
          locsT <- rbind(locs[timeIndexMonth[[j]], ], locs[timeIndexMonth[[j]], ])
          covVals <- extract(cov, SpatialPoints(locsT[,2:1]))[1]
          # outCovVals[timeIndex[[j]]] <- covVals
        }else{
          covVals <- extract(cov, SpatialPoints(locsT[,2:1]))
          # outCovVals[timeIndex[[j]]] <- covVals
        }
        
        #get next month
        unqTimeMonthNext <- unqTimeMonth + c(0, 1)
        if(as.numeric(unqTimeMonthNext[2]) > 12){
             unqTimeMonthNext[2] <- unqTimeMonthNext[2] - 12
             unqTimeMonthNext[1] <- unqTimeMonthNext[1] + 1
        }
        #open current month and extract values
        if(as.numeric(unqTimeMonthNext[2]) < 10){
          covpath <- paste0(startpath, unqTimeMonthNext[1], ".0", as.numeric(unqTimeMonthNext[2]), endpath)
        }else{
          covpath <- paste0(startpath, unqTimeMonthNext[1], ".", unqTimeMonthNext[2], endpath)
        }
        covfull <- tryCatch(raster(covpath),
                            error = function(e){
                              print("test3")
                              covfull
                            })
        cov <- crop(covfull, extent(c(lonmin, lonmax, latmin, latmax)))
        if(length(timeIndexMonth[[j]]) == 1){
          locsT <- rbind(locs[timeIndexMonth[[j]], ], locs[timeIndexMonth[[j]], ])
          covValsNext <- extract(cov, SpatialPoints(locsT[,2:1]))[1]
          # outCovVals[timeIndex[[j]]] <- covVals
        }else{
          covValsNext <- extract(cov, SpatialPoints(locsT[,2:1]))
          # outCovVals[timeIndex[[j]]] <- covVals
        }
        
        for(l in 1:length(timeIndexMonth[[j]])){
          k <- timeIndexMonth[[j]][l]
           # print(k)
          # print(outCovVals[k])
          # print(obsMonth[k, 3])
          # print(covVals * obsMonth[k, 3] + covValsNext * (1 - obsMonth[k, 3]))
          outCovVals[k] <- covVals[l] * obsMonth[k, 3] + covValsNext[l] * (1 - obsMonth[k, 3])
          # print(outCovVals)
        }
        
      }
      if(namesSupplied) outVals[[covnames[[i]]]] <- outCovVals else outVals[[i]] <- outCovVals
    }
    ##if not by month
    else{
      
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
        print(paste0(j, " of ", length(timesIndex)))
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
        covfullnext <- tryCatch(raster(covpathnext),
                                error = function(e){
                                  print("test4")
                                  raster(covpath)
                                })
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
  }
  # print(timeIndexMonth)
  return(outVals)
}

##test
if(runtests){
  # covpaths <- list("Z:/CHAI/Modelling/Mozambique/AI.tif")
  # M <- 50000
  # locs <- cbind(runif(M, -15, -10), runif(M, 30, 37))
  # AI <- getStaticCovs(locs, list("Z:/CHAI/Modelling/Mozambique/AI.tif"))
  # 
  # df <- data.frame(x=locs[,2], y=locs[,1], AI=AI[[1]])
  # p <- ggplot(df, aes(x=x, y=y, color=AI)) + geom_point()
  # print(p)
  # plot(raster(covpaths[[1]]))

  # M <- 10000
  # locs <- cbind(c(runif(M, -15, -12.5), runif(M, -12.5, -10)), c(runif(M, 30, 32), runif(M, 32, 37)))
  # startpaths <- list("Z:/mastergrids/MODIS_Global/MOD11A2_LST/LST_Day/1km/Monthly/LST_Day.")
  # endpaths <- list(".Mean.1km.Data.tif")
  # year <- c(rep(2005, M), rep(2010, M))
  # month <- c(rep(1, M), rep(5, M))
  # v <- getTempCovsMonth(locs, year, month, startpaths, endpaths)
  # 
  # fullpath <- "Z:/mastergrids/MODIS_Global/MOD11A2_LST/LST_Day/1km/Monthly/LST_Day.2005.01.Mean.1km.Data.tif"
  # 
  # df <- data.frame(x=locs[,2], y=locs[,1], AI=v[[1]])
  # p <- ggplot(df, aes(x=x, y=y, color=AI)) + geom_point()
  # print(p)
  
  dts <- c(rep(as.Date("02/12/2014", format="%d/%m/%Y"), 2),
           rep(as.Date("12/03/2015", format="%d/%m/%Y"), 2)
  )
  locs <- cbind(runif(4, -15, -12.5), runif(4, 30, 32))
  startpaths <- list("Z:/mastergrids/MODIS_Global/MOD11A2_LST/LST_Day/1km/Monthly/LST_Day.",
                     "Z:/mastergrids/MODIS_Global/MOD11A2_LST/LST_Day/1km/Monthly/LST_Day.",
                     "Z:/mastergrids/Other_Global_Covariates/Rainfall/CHIRPS/1k/Monthly/chirps-v2-0.",
                     "Z:/mastergrids/MODIS_Global/MCD43B4_BRDF_Reflectance/TCB/1km/8-Daily/",
                     "Z:/mastergrids/MODIS_Global/MCD43B4_BRDF_Reflectance/EVI/1km/8-Daily/")
  endpaths <- list(".Mean.1km.Data.tif",
                   ".Mean.1km.Data.tif",
                   ".1km.tif",
                   "_TCB_Filled_Data.tif",
                   "_EVI_Filled_Data.tif")
  covnames <- list("Temp1", "Temp2", "Rain", "TCB", "EVI")
  byMonth <- c(TRUE, TRUE, TRUE, FALSE, FALSE)
  v <- getTempCovsDate(locs, dts, startpaths, endpaths, covnames, byMonth)
}
