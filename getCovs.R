###functions to automatically get covariate values from mastergrids
library(raster)
library(ggplot2)
runtests <- FALSE

##function to get static covariates
getStaticCovs <- function(locs, covpaths, covnames=NULL){
  outVals <- list()
  if(is.null(covnames)) namesSupplied <- FALSE else namesSupplied <- TRUE
  
  #load each covariate
  i <- 1
  for(covpath in covpaths){
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

  M <- 10000
  locs <- cbind(c(runif(M, -15, -12.5), runif(M, -12.5, -10)), c(runif(M, 30, 32), runif(M, 32, 37)))
  startpaths <- list("Z:/mastergrids/MODIS_Global/MOD11A2_LST/LST_Day/1km/Monthly/LST_Day.")
  endpaths <- list(".Mean.1km.Data.tif")
  year <- c(rep(2005, M), rep(2010, M))
  month <- c(rep(1, M), rep(5, M))
  v <- getTempCovsMonth(locs, year, month, startpaths, endpaths)

  fullpath <- "Z:/mastergrids/MODIS_Global/MOD11A2_LST/LST_Day/1km/Monthly/LST_Day.2005.01.Mean.1km.Data.tif"

  df <- data.frame(x=locs[,2], y=locs[,1], AI=v[[1]])
  p <- ggplot(df, aes(x=x, y=y, color=AI)) + geom_point()
  print(p)
}
