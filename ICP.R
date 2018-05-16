library(Rcpp)
sourceCpp("C:/Users/scro3122/Documents/Mozambique/prewhitenfunctions.cpp")
source("C:/Users/scro3122/Documents/Mozambique/pcalgparallelf.R")


icp <- function(obsDat, locs, times , alpha, target.index, environ.index, G_0 = NULL, supprMessages = FALSE, nSample = length(obsDat[[1]]),
                  plotgraph = TRUE, sampleList=NULL){
  if(!is.null(sampleList)){
    s <- sampleList
  }else{
    s <- sort(sample.int(length(obsDat[[1]]), nSample))
  }
  #print(s)
  ##prewhiten
  #print(times)
  obsDatpw <- list()
  for(i in 1:length(obsDat)){
    obsDatpw[[i]] <- prewhiten(obsDat[[i]][s], locs[s, 1], locs[s, 2], times[s])
  }
  
  ##normalise
  obsDat.n <- list()
  for(i in 1:length(obsDatpw)){
    v <- obsDatpw[[i]]
    obsDat.n[[i]] <- (v - mean(v, na.rm = T)) / sd(v, na.rm = T)
  }
  # print("o")
  # print(obsDatpw)
  
  ##get all subsets
  ss <- getAllSubsets((1:length(obs))[-c(target.index, environ.index)])
  print("total number of subsets")
  print(length(ss))
  
  ##test indep(E, Target) given each subset - save the ones that are
  ss.working <- list()
  count <- 1
  for(s in ss){
    print(paste0(count, " of ", length(ss)))
    count <- count + 1
    condData <- do.call(cbind, obsDat.n[s])
    if(!condtest(obsDat.n[[target.index]], obsDat.n[[environ.index]], condData, alpha=alpha)){
      ss.working[[length(ss.working) + 1]] <- s
    }
  }
  
  return(list(ss.working, listIntersect(ss.working), definingSets(ss.working)))
}


getAllSubsets <- function(set){
  ss <- list()
  for(i in 1:length(set)){
    ss <- c(ss, getSubsets(set, i))
  }
  return(ss)
}

listIntersect <- function(lInt){
  i <-lInt[[1]]
  for(j in 2:length(lInt)){
    i <- intersect(i, lInt[[j]])
  }
  return(i)
}

definingSets <- function(setList){
  ##get all elements in any list
  e <- setList[[1]]
  for(j in 2:length(setList)){
    e <- union(e, setList[[j]])
  }
  
  #get all subsets of e
  e.ss <- getAllSubsets(e)
  
  #check if each subset has this property
  e.ss.p <- list()
  for(e.s in e.ss){
    p <- TRUE
    for(setL in setList){
      if(!any(e.s %in% setL)){
        p <- FALSE
        break
      }
    }
    if(p) e.ss.p <- c(e.ss.p, list(e.s))
  }
  
  #for subset with property check if it is minimal
  e.ss.p.m <- list()
  for(i in 1:length(e.ss.p)){
    settest <- e.ss.p[[i]]
    minimal <- TRUE
    for(j in (1:length(e.ss.p))[-i]){
      settest2 <- e.ss.p[[j]]
      if(all(settest2 %in% settest)){
        minimal <- FALSE
        break
      }
    }
    if(minimal) e.ss.p.m <- c(e.ss.p.m, list(settest))
  }
  return(e.ss.p.m)
}