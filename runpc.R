library(Rcpp)
sourceCpp("C:/Users/scro3122/Documents/Mozambique/prewhitenfunctions.cpp")
source("C:/Users/scro3122/Documents/Mozambique/pcalgparallelf.R")


runpc <- function(obsDat, locs, times , alpha, G_0 = NULL, supprMessages = FALSE, nSample = length(obsDat[[1]]),
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
    obsDatpw[[i]] <- prewhiten(obsDat[[i]][s], locs[s, 1], locs[s, 2], times[s], alpha=0.9)
  }
  
  ##normalise
  obsDat.n <- list()
  for(i in 1:length(obsDatpw)){
    v <- obsDatpw[[i]]
    obsDat.n[[i]] <- (v - mean(v, na.rm = T)) / sd(v, na.rm = T)
  }
  # print("o")
  # print(obsDatpw)
  pc <- pcalg(obsDat.n, alpha, G_0 = G_0, supprMessages = supprMessages)
  if(plotgraph){
    library(graph)
    pcmat <- pc
    rownames(pcmat[[1]])<- names(obsDat)
    colnames(pcmat[[1]])<- names(obsDat)
    rownames(pcmat[[2]])<- names(obsDat)
    colnames(pcmat[[2]])<- names(obsDat)
    am.graph <-new("graphAM", adjMat=pcmat[[1]], edgemode="directed")
    am.graph2 <-new("graphAM", adjMat=pcmat[[2]], edgemode="directed")
    
    
    plot(am.graph, attrs = list(node = list(fillcolor = "white"),
                                edge = list(arrowsize=0.5)), main = paste0("Full_", length(s), "_pts_",alpha))


    plot(am.graph2, attrs = list(node = list(fillcolor = "white"),
                                 edge = list(arrowsize=0.5)), main = "Skeleton_parallel")
  }
  
  names(obsDatpw) <- names(obsDat)
  return(c(pc, list(obsDatpw)))
}

plot.minimal <- function(adjMat, targetIndex, obsNames=NULL){
  nvar <- dim(adjMat)[1]
    
  #find parents of target variable
  ps <- c()
  for(i in (1:nvar)[-targetIndex]){
    if(adjMat[i, targetIndex] || adjMat[targetIndex, i]){
      ps <- c(ps, i)
    }
  }
  
  #delete edges not between parent and target or parent and parent
  if(length(ps) < nvar){
    for(i in (1:nvar)[-c(ps, targetIndex)]){
      for(j in (1:nvar)[-c(ps, targetIndex)]){
        adjMat[i, j] <- 0
        adjMat[j, i] <- 0
      }
    }
  }
  adjMat2 <- adjMat[c(ps, targetIndex), c(ps, targetIndex)]
  if(!is.null(obsNames)){
    rownames(adjMat2) <- obsNames[c(ps, targetIndex)]
    colnames(adjMat2) <- obsNames[c(ps, targetIndex)]
  }
  
  library(graph)
  for(i in (1:nvar)[-targetIndex]){
    for(j in (1:nvar)[-targetIndex]){
      adjMat[i, j] <- 0
      adjMat[j, i] <- 0
    }
  }
  
  am.graph <-new("graphAM", adjMat=adjMat2, edgemode="directed")
  
  plot(am.graph, attrs = list(node = list(fillcolor = "white"),
                              edge = list(arrowsize=0.5)), main = paste0("Minimal"))
  
}
