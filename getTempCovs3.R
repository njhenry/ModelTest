# obsdates <- c(rep(as.Date("02/12/2014", format="%d/%m/%Y"), 2),
#               rep(as.Date("12/03/2015", format="%d/%m/%Y"), 2)
# )

obsdates <- as.Date(c("02/12/2014", "09/12/2014", "16/12/2014", "23/12/2014"),
                    format="%d/%m/%Y")
locs <- cbind(rep(runif(2, -15, -12.5), 4), rep(runif(2, 30, 32), 4))
folderpaths <- list("Z:/mastergrids/MODIS_Global/MOD11A2_LST/LST_Day/1km/Monthly/",
                    "Z:/mastergrids/MODIS_Global/MCD43B4_BRDF_Reflectance/TCB/1km/8-Daily/",
                   "Z:/mastergrids/MODIS_Global/MCD43B4_BRDF_Reflectance/EVI/1km/8-Daily/")
endpaths <- list(".Mean.1km.Data.tif",
                 "_TCB_Filled_Data.tif",
                 "_EVI_Filled_Data.tif")
covnames <- list("LST", "TCB", "EVI")
byMonth <- c(TRUE, FALSE, FALSE)
timelags <- c(1,1,1)

library(raster)
library(ggplot2)
library(lubridate)

getTempCovsDate3 <- function(locs, obsdates, folderpaths, endpaths, covnames=NULL, byMonth=rep(TRUE, length(folderpaths)),
                             timelags=rep(list(NULL), length(folderpaths))){
  outVals <- list()
  if(is.null(covnames)) namesSupplied <- FALSE else namesSupplied <- TRUE
  
  ##create full list of names with time lags
  if(namesSupplied){
    covnamesfull <- list()
    count <- 1
    for(i in 1:length(covnames)){
      covnamesfull[[count]] <- covnames[[i]]
      count <- count + 1
      if(!is.null(timelags[[i]])){
        for(timelag in timelags[[i]]){
          covnamesfull[[count]] <- paste(covnames[[i]], timelag, sep="-")
          count <- count + 1
        }
      }
    }
  }
  # print(covnamesfull)
  # print(covnames)
  # print("test")
  
  ##for each covariate]
  for(i in 1:length(folderpaths)){
    if(namesSupplied) print(covnames[[i]]) else print(paste0("cov ", i))
    
    
    ###deal with timelags
    if(!is.null(timelags[[i]])){
      obsdates.w.lag <- obsdates
      for(timelag in timelags[[i]]){
        obsdates.lag <- obsdates
        week(obsdates.lag) <- week(obsdates) - timelag
        obsdates.w.lag <- c(obsdates.w.lag, obsdates.lag)
      }
      obsdates.use <- obsdates.w.lag
      locs.use <- matrix(rep(as.vector(t(locs)), length(timelags[[i]]) + 1), ncol=2, byrow=TRUE)
    }else{
      obsdates.use <- obsdates
      locs.use <- locs
    }
    print(obsdates.use)
    print(locs.use)
    outCovVals <- rep(0, length(obsdates.use))
    folderpath <- folderpaths[[i]]
    endpath <- endpaths[[i]]
    
    #get list of dates available
    file.names <- list.files(folderpath, pattern=paste0(endpath, "$") )
    if(byMonth[i]){
      file.year.month <- regmatches(file.names, regexpr("[0-9]{4}\\.[0-9]{2}", file.names))
      file.dates <- as.Date(paste(file.year.month, ".15", sep=""), format="%Y.%m.%d")
      # print(nearestFile)
      
    }else{
      file.dates <- as.Date(regmatches(file.names, regexpr("[0-9]{7}", file.names)), format="%Y%j")
      
      #sort file names by date in case not already in the right order
      file.names <- file.names[order(file.dates)]
      file.dates <- sort(file.dates)
    }
    
    
    
    ##for each obsdate associate two files, before and after, and the weighting of each
    nearestFile2 <- c()
    weights <- c()
    nonneg <- function(vec){
      vec[vec < 0] <- NA
      return(vec)
    }
    
    for(k in 1:length(obsdates.use)){
      obsdate <- obsdates.use[k]
      #file before date
      nearestFile2[2*k-1] <- which.min(nonneg(obsdate - file.dates)) 
      before.date <- file.dates[nearestFile2[2*k-1]]
      #file after date
      nearestFile2[2*k] <- which.min(nonneg(file.dates - obsdate))
      after.date <- file.dates[nearestFile2[2*k]]
      
      #weights
      a <- as.numeric(after.date - obsdate) / as.numeric(after.date - before.date)
      b <- as.numeric(obsdate - before.date) / as.numeric(after.date - before.date)
      
      weights <- rbind(weights, c(a, b))
    }
    
    nearestIndex2 <- sapply(1:length(file.names), function(k) ceiling(which(nearestFile2 == k)/2))
    nearestIndex2After <- sapply(1:length(file.names), function(k) which(nearestFile2==k) %% 2)
    nearestIndex2BeforeAfter <- sapply(nearestIndex2After, function(k) if(length(k) > 0) rbind(k, !k))
    
    #for each obsdate find which of the files are the nearest
    nearestFile <- c()
    for(k in 1:length(obsdates.use)){
      nearestFile[k] <- which.min(abs(file.dates - obsdates.use[k]))
    }
    nearestIndex <- sapply(1:length(file.names), function(k) which(nearestFile==k))
    
    #find number of files that need to be opened
    nonempty <- c()
    for(j in 1:length(nearestIndex)){
      if(length(nearestIndex2[[j]] > 0)) nonempty <- c(nonempty, j)
    }
    #for each file that has a non-zero number of nearest obsdates, open that covariate and
    #extract data at locations
    for(l in 1:length(nonempty)){
      j <- nonempty[l]
      print(paste0(l, " of ", length(nonempty)))
      nearIndex <- nearestIndex2[[j]]
      if(length(nearIndex) > 0){
        covpath <- paste0(folderpath, file.names[j])
        # print(covpath)
        cov <- raster(covpath)
        #extract values
        cov.vals <- extract(cov, SpatialPoints(locs.use[nearIndex, 2:1, drop=FALSE]))
        outCovVals[nearIndex] <- outCovVals[nearIndex] + 
          rowSums(t(nearestIndex2BeforeAfter[[j]]) * weights[nearIndex, ]) * cov.vals
        rm(cov)
      }
    }
    print(outCovVals)
    
    # outCovVals2 <- rep(NA, length(obsdates.use))
    # ##check old code
    # for(l in 1:length(nonempty)){
    #   j <- nonempty[l]
    #   print(paste0(l, " of ", length(nonempty)))
    #   nearIndex <- nearestIndex[[j]]
    #   if(length(nearIndex) > 0){
    #     covpath <- paste0(folderpath, file.names[j])
    #     # print(covpath)
    #     cov <- raster(covpath)
    #     #extract values
    #     cov.vals <- extract(cov, SpatialPoints(locs.use[nearIndex, 2:1, drop=FALSE]))
    #     outCovVals2[nearIndex] <- cov.vals
    #   }
    # }
    # print(outCovVals2)
    # 
    if(!is.null(timelags[[i]])){
      print(outCovVals)
      ##separate out timelags
      if(namesSupplied){
        outVals[[covnames[[i]]]] <- outCovVals[1:length(obsdates)]
        for(j in 1:length(timelags[[i]])){
          timelag <- timelags[[i]][j]
          outVals[[paste(covnames[[i]], timelag, sep="-")]] <- outCovVals[(1:length(obsdates)) + j*length(obsdates)]
        }
      }else{
        outVals[[length(outVals)+1]] <- outCovVals[1:length(obsdates)]
        for(j in 1:length(timelags[[i]])){
          timelag <- timelags[[i]][j]
          outVals[[length(outVals)+1]] <- outCovVals[(1:length(obsdates)) + j*length(obsdates)]
        }
      }
      print(outVals)
    }else{
      if(namesSupplied) outVals[[covnames[[i]]]] <- outCovVals else outVals[[length(outVals) + 1]] <- outCovVals
    }

  }
  return(outVals)
}





# 
# timelags <- list(c(1,2), c(1,2,3), NULL)
# v <- getTempCovsDate3(locs, obsdates, folderpaths, endpaths, byMonth=byMonth, timelags=timelags, covnames=covnames)