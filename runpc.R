library(Rcpp)
sourceCpp("C:/Users/scro3122/Documents/Mozambique/prewhitenfunctions.cpp")
source("C:/Users/scro3122/Documents/Mozambique/pcalgparallelf.R")


runpc <- function(obsDat, locs, times , alpha, G_0 = NULL, supprMessages = FALSE, nSample = length(obsDat[[1]]),
                  plotgraph = TRUE, sampleList=NULL, pw=TRUE, returnSeed=FALSE, seed=NULL){
  ptm <- proc.time()
  #if seed supplied use it
  if(!is.null(seed)){
    set.seed(seed)
  }else if(returnSeed){
    seed <- runif(1)
    set.seed(seed)
  }
  #if returnSeed, generate and save a seed and return it at the end


  
  ##if pw (?prewhiten) is false we dont need to use locations and times at all
  ##or test them 
  
  ##check inputs
  #length of vectors in obsDat should be the same
  for(i in 2:length(obsDat)){
    if(length(obsDat[[1]]) != length(obsDat[[i]])) stop(paste0("observations 1 and ", i,  "are not of the same length"))
  }
  
  #check other lengths
  if(pw) if(length(obsDat[[1]]) != dim(locs)[1]) stop("different number of observations and locations")
  if(pw) if(length(obsDat[[1]]) != length(times)) stop("different number of observations and times")
  if(!is.null(G_0)) if(dim(G_0)[1] != dim(G_0)[2]) stop("G_0 is not square")
  if(!is.null(G_0)) if(length(obsDat) != dim(G_0)[1]) stop("G_0 is of the wrong dimension")
  if(!is.null(G_0)) if(sum(diag(G_0)) > 0) stop("invalid G_0, nodes cannot point to themselves")
  
  
  if(!is.null(sampleList)){
    s <- sampleList
  }else{
    s <- sort(sample.int(length(obsDat[[1]]), nSample))
  }
    #print(s)
  ##prewhiten
  #print(times)
  obsDatpw <- list()
  if(pw){
    for(i in 1:length(obsDat)){
      obsDatpw[[i]] <- prewhiten(obsDat[[i]][s], locs[s, 1], locs[s, 2], times[s], alpha=0.9)
    }
  }else{
    for(i in 1:length(obsDat)){
      obsDatpw[[i]] <- obsDat[[i]][s]
    }
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
    
    plot(am.graph2, attrs = list(node = list(fillcolor = "white"),
                                 edge = list(arrowsize=0.5)), main = "Skeleton_parallel")
    
    plot(am.graph, attrs = list(node = list(fillcolor = "white"),
                                edge = list(arrowsize=0.5)), main = paste0("Full_", length(s), "_pts_",alpha))


  }
  
  names(obsDatpw) <- names(obsDat)
  print(proc.time() - ptm)
  if(returnSeed) return(c(pc, list(obsDatpw), seed))
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


plot.minimal.n <- function(adjMat, targetIndex, obsNames=NULL, nDeg = 1,
                           fontsize=20, height=2, color.list=NULL,
                           fontcolor="black"){
  
  
  edges.keep <- c(targetIndex)
  edges.deg <- list()
  for(i in 1:nDeg){
    #print(adjMat[, edges.keep, drop=FALSE])
    edges.keep.old <- edges.keep
    #print(edges.keep.old)
    edges.keep <- unique(c(edges.keep, which(rowSums(adjMat[, edges.keep, drop=FALSE]) > 0))) 
    #print(setdiff(edges.keep, edges.keep.old))
    edges.deg[[i]] <- setdiff(edges.keep, edges.keep.old)
  }
  
  
  ##define some colors
  if(is.null(color.list)) color.list <- c("lightblue", "honeydew", "white", "blue", "gold", "lime")
  n.colors <- length(color.list)
  
  #create new adjacency matrix
  adjMat2 <- adjMat[edges.keep, edges.keep, drop=FALSE]
  if(!is.null(obsNames)){
    rownames(adjMat2) <- obsNames[edges.keep]
    colnames(adjMat2) <- obsNames[edges.keep]
  }
  
  
  #if the adjacency matrix has weights remove them
  adjMat2_use <- (adjMat2 != 0) * 1
  am.graph <-new("graphAM", adjMat=adjMat2_use, edgemode="directed")
  
  modn <- function(m, n){
    if(m==n){
      return(n)
    }else{
      return(m %% n)
    }
  }
  
  ##do colors for different degrees
  l <- c()
  names.l <- c()
  for(i in 1:nDeg){
    l <- c(l, rep(color.list[modn(i, n.colors)], length(edges.deg[[i]])))
    names.l <- c(names.l, obsNames[edges.deg[[i]]])
  }
  
  
  names(l) <- names.l
  l <- c(l, "firebrick")
  names(l) <- c(names.l, obsNames[targetIndex])
  #print(edges.keep)
  nAttrs <- list()
  nAttrs$fillcolor <- l
  nAttrs$fontcolor <- c(cases="white")
  
  
  #nAtrrs$fillcolor <- c(cases="blue")
  plot(am.graph, nodeAttrs=nAttrs, attrs=list(node=list(fillcolor="lightblue",
                                                        fontsize=fontsize,
                                                        height=height,
                                                        fontcolor=fontcolor),
                                              edge=list(arrowsize="0.5")))
  return(list(am.graph, edges.deg, adjMat2, edges.keep))
}

