rm(list=ls(all=TRUE))
load("Data.Rdata")
source("C:/Users/Marc/SkyDrive/R/Thesis/Thesis/Prediction_functions.R")
source("C:/Users/Marc/SkyDrive/R/Thesis/Thesis/t2maturity.r")

colnames(Data) = as.numeric(colnames(Data))  + 1
Y <- cbind(maturity2tenor(Data[,1:522], rowDates = as.Date(rownames(Data)),
                          EndDatesBrent, bizdaysList)[,1:24],
           maturity2tenor(Data[,523:1046], rowDates = as.Date(rownames(Data)),
                          EndDatesNap, bizdaysList)[,1:24])

res_path <- "D:/R/Results/44T/"

load(paste(res_path, "results44T_1_1.Robj", sep=""))
load(paste(res_path, "results44T_1_20.Robj", sep=""))
load(paste(res_path, "results44T_1_50.Robj", sep=""))
load(paste(res_path, "updates44T_1.Robj", sep=""))
load(paste(res_path, "updates44T_50.Robj", sep=""))
load(paste(res_path, "updates44T_20.Robj", sep=""))

perUp <- 10

myLoadYp <- function(Index, Ps, res_name){
  dat_name <- paste(res_name, Ps, Index,  sep="_")
  load(paste(res_path, dat_name, ".Robj", sep=""))
  return(get(dat_name)$Yp)
}

YpL <- lapply(c(1:20,34:53, 68:87), myLoadYp, Ps=22, res_name = "results44T")
Yp <- array(0, dim = c(dim(YpL[[1]])[1:2], length(YpL)*perUp ))

for(i in 1:length(YpL)){
  Yp[,, (1:perUp) + (i-1)*perUp] <- YpL[[i]]
}

rm("YpL")
gc()
plot(density(Yp[1,1,]))
summary(Yp[1,1,])
YpMeans <- apply(Yp, 2, rowMeans)

diffs <- Y[-(1:23),]-YpMeans

diffs2 <- Y[-(1:23),c(1:23,25:47)]-Y[-((nrow(Y)-22):nrow(Y)),c(2:24,26:48)]

Yp2 <- createYp( array(0, c(1259,8,10)), postBetaAll[5000,,],
                 Phit = Phit, d= postDAll[5000,,], C= 2, K=4,
                 Et = postEtAll[5000,])


Beta <- t(apply(postBetaAll[2000:7000,,],2, colMeans))
bdiffs <- Beta - updates44T_1[[1]][1:1260,]
