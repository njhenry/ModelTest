source("C:/Users/scro3122/Documents/Mozambique/pcalgparallelf.R")

runpc <- function(obsDat, locs, times , alpha, G_0 = NULL, supprMessages = FALSE, nSample = length(obsDat[[1]])){
    s <- sample.int(length(obsDat[[1]]), nSample)
    print(s)
  ##prewhiten
  obsDatpw <- list()
  for(i in 1:length(obsDat)){
    obsDatpw[[i]] <- gpReg(c(locs[s,], times[s]), obsDat[[i]][s])
  }
  return(pcalg(obsDatpw, alpha, G_0 = G_0, supprMessages = supprMessages))
}