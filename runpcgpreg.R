library(Rcpp)
sourceCpp("C:/Users/scro3122/Documents/Mozambique/prewhitenfunctions.cpp")
source("C:/Users/scro3122/Documents/Mozambique/pcalgparallelf.R")


runpc2 <- function(obsDat, locs, times , alpha, G_0 = NULL, supprMessages = FALSE, nSample = length(obsDat[[1]]),
                  plotgraph = TRUE, sampleList=NULL){
  if(!is.null(sampleList)){
    s <- sampleList
  }else{
    s <- sample.int(length(obsDat[[1]]), nSample)
  }
    #print(s)
  ##prewhiten
  #print(times)
  obsDatpw <- list()
  for(i in 1:length(obsDat)){
    obsDatpw[[i]] <- gpReg(cbind(locs, times)[s], obsDat[[i]][s])
  }
  # print("o")
  # print(obsDatpw)
  pc <- pcalg(obsDatpw, alpha, G_0 = G_0, supprMessages = supprMessages)
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
  
  return(pc)
}

