rm(list=ls(all=TRUE))
library(ggplot2)
load("Data.Rdata")
source("C:/Users/Marc/SkyDrive/R/Thesis/Thesis/Prediction_functions.R")
source("C:/Users/Marc/SkyDrive/R/Thesis/Thesis/t2maturity.r")
Path <- "D:/R/Results/"
Data["2015-09-04", 870:1046] <- NA # 25% percent change in one day -> outlier

colnames(Data) = as.numeric(colnames(Data))  + 1

Y <- cbind(maturity2tenor(Data[, 1:522], rowDates = as.Date(rownames(Data)),
                          EndDatesBrent, bizdaysList)[,1:24],
           maturity2tenor(Data[, 523:1046], rowDates = as.Date(rownames(Data)),
                          EndDatesNap, bizdaysList)[,1:24])


save(list = "Y", file = paste(Path, "Y.RObj", sep=""))

for( i in c(1, 5, 22)){
    diffs2Name <- paste("diffs2", "_", i, sep="")
  assign(diffs2Name, diffs2fn(i, Y))
  save(list = diffs2Name,
       file = paste(Path, diffs2Name, ".RObj", sep=""))
}
perUp <- 10

combinations <- paste(rep(c("44T", "33T"), c(3,3)),
                      rep(c(1,5,22), 2),
                      sep="_")

lapply(combinations, evalRun, Path = "D:/R/Results/", perUp = perUp, Y = Y)

