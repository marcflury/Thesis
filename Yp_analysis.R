rm(list=ls(all=TRUE))
res_path <- "D:/R/Results/"
load("D:/R/Results/results44T_1_1.Robj")
load("D:/R/Results/results44T_22_1.Robj")
load("D:/R/Results/updates44T_1.Robj")
load("D:/R/Results/updates44T_2.Robj")
perUp <- 5

hist(updates44T_1$Beta - updates44T_2$Beta)

myLoadYp <- function(Index, Ps, res_name){
  dat_name <- paste(res_name, Ps, Index,  sep="_")
  load(paste(res_path, dat_name, ".Robj", sep=""))
  return(get(dat_name)$Yp)
  
}

YpL <- lapply(1:100, myLoadYp, Ps=22, res_name = "results44T")
Yp <- array(0, dim = c(dim(YpL[[1]])[1:2], length(YpL)*perUp ))

for(i in 1:length(YpL)){
  Yp[,, (1:perUp) + (i-1)*perUp] <- YpL[[1]]
}
rm("YpL")
gc()

YpMeans <- apply(Yp, 2, rowMeans)
