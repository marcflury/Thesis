rm(list = ls(all=TRUE))
# Required packages:
library(fda)
library(KFAS)
library(MCMCpack)
library(stochvol)
library(truncnorm)
library(parallel)
library(timeDate)
library(dplyr)
library(tidyr)
library(lubridate)
library(xts)

setwd("C:/Users/Marc/SkyDrive/R/Thesis/Thesis/")
set.seed(42)
load("EstResultsFull3HMM.RObj")
load("EstResults3HMM.RObj")
load("Data.RData")
source("MF_MFDLM_ARIMA.R")
source("Prediction_functions.R")
source("t2maturity.r")
K <- 3
C <- 2
postBetaAll <- EstResultsFull3HMM$postBetaAll
postDAll <- EstResultsFull3HMM$postDAll
postLambdaAll <- EstResultsFull3HMM$postLambdaAll
postEtAll <- EstResultsFull3HMM$postEtAll
postHtAll <- EstResultsFull3HMM$postHtAll
postGammaSlopesAll <- EstResultsFull3HMM$postGammaSlopesAll
postSAll <- EstResultsFull3HMM$postSAll

postsvMu <- EstResultsFull3HMM$postsvMu
postsvPhi <- EstResultsFull3HMM$postsvPhi
postsvSigma <- EstResultsFull3HMM$postsvSigma
postPsi <- EstResultsFull3HMM$postPsi
postq10 <- EstResultsFull3HMM$postq10
postq01 <- EstResultsFull3HMM$postq01

rm(list = c("EstResultsFull3HMM"))
gc()


Data["2015-09-04", 870:1046] <- NA # 25% percent change in one day -> outlier
Y <- Data
rowDatesY <- rownames(Data)
T = nrow(Y)							# Total number of time points
tau0 = as.numeric(colnames(Y))			# Frequencies within each time bin, concatenated across c = 1,...,C
tau = (tau0 - min(tau0))/(diff(range(tau0)))	# Shift to [0,1] (for numerical stability)
allTaus0 = sort(unique(tau0))				# All unique frequencies (observation points)
allTaus = sort(unique(tau)) 				# All unique frequencies (observation points), restricted to [0,1]
m = length(allTaus)					# Total number of unique frequencies

splineInfo <- getSplineInfo(tau)
# Some useful indices:
outcomeInds = c(1, which(tau[-1] < tau[-length(tau)])+1, length(tau)+1)
# For outcome-specific observations, use Y[, yc.inds]
# where yc.inds = outcomeInds[c]:(outcomeInds[c+1]-1)
betaInds = seq(1, C*(K+1), by=K)								# For outcome-specific Beta, use  Beta[, bc.inds]
# where bc.inds = betaInds[c]:(betaInds[c+1]-1)

####################################
Phit = splineInfo$basisPhi(seq(splineInfo$a, splineInfo$b, length.out=523))
# Interpolate the data with using cubic splines
Y.full <- do.call("cbind",lapply(1:C,splineinterpol))
T = nrow(Y.full)

Y <- do.call("cbind",lapply(1:C,splineinterpol))

# MCMC parameters:
nsims =  7000		# Total number of simulations
burnin = 2000		# Burn-in

updateNum <- 100
postIs <- burnin + floor(runif(updateNum)*(nsims-burnin))

# Calculate the number of cores
num_cores <- detectCores() - 1


print(Sys.time())
if(TRUE){
  
  cl <- makeCluster(num_cores)
  
  clusterExport(cl, varlist=ls())
  t0 <- Sys.time()
  parLapply(cl, 67:100, fullPrediction)
  t1 <- Sys.time()
  t1 - t0
  stopCluster(cl)
  
}
save("postIs", file="postIs33T.Robj")
