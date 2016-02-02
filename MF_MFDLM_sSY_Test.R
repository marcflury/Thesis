rm(list=ls(all=TRUE)) 

set.seed(14850) # to reproduce (most of) the results from the paper

#########################################################################################################################

# MCMC parameters:
nsims =  1e4 		# Total number of simulations
burnin = 0.1*nsims		# Burn-in

# FDLM parameters:
K = 4			# Number of factors
K.hmm.sv = 4		# Number of factors for the common trend model
useHMM = FALSE		# Hidden Markov model (HMM), or common trend (CT) model? (CT in the paper)
# Note: HMM needs additional adjustments to store the relevant parameters
#########################################################################################################################

# Required packages:
library(fda)
library(KFAS)
library(MCMCpack)
library(stochvol)
library(truncnorm) 
#library(dlm)

# File containing functions
source("MF_MFDLM.R")

# This may improve efficiency in some cases:
#library(compiler);  enableJIT(3)

#########################################################################################################################

# Load the data and store key variables:
load("Data.RData"); Y = diff(Data)			# Difference the data for easier analysis, stationarity

T = nrow(Y)							# Total number of time points
tau0 = as.numeric(colnames(Y))			# Frequencies within each time bin, concatenated across c = 1,...,C
tau = (tau0 - min(tau0))/(diff(range(tau0)))	# Shift to [0,1] (for numerical stability)
allTaus0 = sort(unique(tau0))				# All unique frequencies (observation points)
allTaus = sort(unique(tau)) 				# All unique frequencies (observation points), restricted to [0,1]
m = length(allTaus)					# Total number of unique frequencies

C = 2							# Number of outcomes
cnames = c("Brent", "Naphtha")		# Names of the outcomes (central banks)
dates = as.Date(rownames(Y))				# Dates of the observations

#########################################################################################################################

# Some useful indices:
outcomeInds = c(1, which(tau[-1] < tau[-length(tau)])+1, length(tau)+1)		
# For outcome-specific observations, use Y[, yc.inds] 
# where yc.inds = outcomeInds[c]:(outcomeInds[c+1]-1)
betaInds = seq(1, C*(K+1), by=K)								# For outcome-specific Beta, use  Beta[, bc.inds]
# where bc.inds = betaInds[c]:(betaInds[c+1]-1)

#########################################################################################################################

# Initialize the main parameters:
inits = initParams(Y, tau, K=K, useAllTaus=TRUE)
Beta = inits$Beta0
d = inits$d0
Et = inits$Et0
lambda = inits$lambda0
splineInfo = inits$splineInfo
K = inits$K
betaInds = seq(1, C*(K+1), by=K)

# Initialize the AR(1) and HMM parameters:
initsHMMpar = initHMM(Beta, K.hmm.sv, useHMM)
psi = initsHMMpar$psi
q01 = initsHMMpar$q01
q10 = initsHMMpar$q10
S = initsHMMpar$S
gammaSlopes = initsHMMpar$gammaSlopes

# Initialize the SV parameters (just using independent AR(1) models for Beta_ck):
initsSVpar = initSV(apply(Beta, 2, function(x){arima(x, c(1,0,0), include.mean=FALSE)$resid[-1]}), C)
ht = initsSVpar$ht
svMu = initsSVpar$svMu
svPhi = initsSVpar$svPhi
svSigma = initsSVpar$svSigma
Wt = initsSVpar$Wt

#########################################################################################################################

# Set up the state space model
Gt = array(0, c(C*K, C*K, T)); Gt[1:(C*K), 1:(C*K),] = diag(psi); Wt = initsSVpar$Wt

Model = SSModel(Y~-1+SSMcustom(Z = array(0, c(ncol(Y), nrow(Gt))),
                               T = Gt, Q = Wt, P1 = diag(10^4, nrow(Gt))))

#########################################################################################################################

# Parameters to save:
postBetaAll = array(0, c(nsims, dim(Beta))); postDAll = array(0, c(nsims, dim(d)))
postLambdaAll = array(0, c(nsims, length(lambda))); postEtAll = array(0, c(nsims, C))  
postHtAll = array(0, c(nsims, dim(ht))); postGammaSlopesAll = array(0, c(nsims, C*K))
devAll = numeric(nsims)  # Deviance
postCountAll = array(0, c(nsims, T-1, C*K)); postCountAggAll = array(0, c(nsims, T-1, C-1)); # For outlier plot
#########################################################################################################################

# Now, run the MCMC
timer0 = proc.time()[3]			# for timing the sampler
for(nsi in 1:nsims){			# nsi is the simulation index
  
  # Sample Beta, d, and lambda:
  samples = mfdlm(Y, tau, Beta, Et, Gt, Wt, Model, d, splineInfo, lambda)
  Beta = samples$Beta; d = samples$d; lambda=samples$lambda
  
  # Cycle through the outcomes:
  for(c in 1:C) {
    bc.inds = betaInds[c]:(betaInds[c+1]-1)		# Subsets: \beta_{1,t}^{(c)}, ..., \beta_{K,t}^{(c)}
    yc.inds = outcomeInds[c]:(outcomeInds[c+1]-1) 	# Subsets: Y_t^{(c)}(\tau_1), ..., Y_t^{(c)}(\tau_m)
    
    # Conditional means of Y_t^{(c)}(\tau) (T x m)
    mu.c = tcrossprod(Beta[,bc.inds], splineInfo$Phi[match(tau[yc.inds], allTaus),]%*%d)
    
    # Observation error variance
    Et[c] = sampleEt(Y[, yc.inds], mu.c)
    
    # Sample the AR(1) and slope parameters (and other HMM parameters, if desired)
    shmm = sampleHMMpar(Beta, K.hmm.sv, S, gammaSlopes, psi, ht, c, useHMM)
    gammaSlopes=shmm$gammaSlopes
    psi = shmm$psi
    
    # Sample the stochastic volatility parameters:
    if(c > 1){
      resBeta.c = Beta[, bc.inds] - S[,bc.inds]*Beta[,1:K]*matrix(rep(gammaSlopes[bc.inds], T), nrow=T, byrow=TRUE)
      resBeta.c = resBeta.c[-1,] -  matrix(rep(psi[bc.inds], T-1), nrow=T-1, byrow=TRUE)*resBeta.c[-nrow(resBeta.c),]
    } else resBeta.c = Beta[-1, bc.inds] - matrix(rep(psi[bc.inds], T-1), nrow=T-1, byrow=TRUE)*Beta[-T, bc.inds]
    
    samples = sampleSV(resBeta.c,  ht[, bc.inds], svMu[bc.inds], svPhi[bc.inds], svSigma[bc.inds], svOffset=10^-6)
    ht[, bc.inds] <- samples$ht.c
    svMu[bc.inds] <- samples$svMu.c
    svPhi[bc.inds] <- samples$svPhi.c
    svSigma[bc.inds] <- samples$svSigma.c
    Wt[bc.inds, bc.inds,] <- samples$Wt.c
    
    if(c > 1){
      # Standardized residual diagnostic/outlier plot:
      sResid2 = resBeta.c^2/exp(ht[-1, bc.inds])
      postCountAll[nsi,,bc.inds] = sResid2 > qchisq(0.95, 1)
      postCountAggAll[nsi,, c-1] = rowSums(sResid2[,1:K.hmm.sv]) > qchisq(0.95, K.hmm.sv)
    }
    
    # Compute the (observation-level) deviance, summing over c=1,...,C
    devAll[nsi] = devAll[nsi] + 1/Et[c]*sum((Y[, yc.inds]-mu.c)^2,na.rm=TRUE) + sum(!is.na(Y[, yc.inds]))*log(2*pi*Et[c])
  }
  
  # Compute Gt and Wt matrices for HMM; common trend model is a submodel (Note: this is not efficient)
  GtWt <- computeGtWtHMM(Gt, Wt, S, gammaSlopes, psi, exp(ht))
  Gt <- GtWt$Gt
  Wt <- GtWt$Wt
  
  # Store the samples, respecting the burn-in and thinning
  postBetaAll[nsi,,] <- Beta
  postDAll[nsi,,] <- d
  postLambdaAll[nsi,] <- lambda
  postEtAll[nsi,] <- Et
  postHtAll[nsi,,] <- ht
  postGammaSlopesAll[nsi,] <- gammaSlopes 
  
  # Check the time remaining:
  computeTimeRemaining(nsi, timer0, nsims)
}
print(paste('Total time: ', round((proc.time()[3] - timer0)/60), 'minutes'))

#########################################################################################################################
# Discard burn-in in duplicate arrays (for some variables)
postBeta <- postBetaAll[-(1:burnin),,]; postD = postDAll[-(1:burnin),,]
postEt = postEtAll[-(1:burnin),]
postHt = postHtAll[-(1:burnin),,]
postGammaSlopes = postGammaSlopesAll[-(1:burnin),]
dev = devAll[-(1:burnin)];
# postCount = postCountAll[-(1:burnin),,]
# postCountAgg = postCountAggAll[-(1:burnin),,]

# Fix important parameters at their posterior means:
Beta = colMeans(postBeta); d = colMeans(postD)
Et = colMeans(postEt); ht = colMeans(postHt)
gammaSlopes = colMeans(postGammaSlopes)


# Test Data

# Plot the Joint loading curves:
taugrid = seq(allTaus0[1], allTaus0[m], length.out=24)
Phit = splineInfo$basisPhi(seq(splineInfo$a, splineInfo$b, length.out=length(taugrid)))
par(mfrow=c(1,1), mai=c(1.1,1,1,.1))
plot(taugrid, Phit%*%d[,1], type='l', lwd=10, ylim=1.1*range(Phit%*%d[,1:K.hmm.sv]), xlab=expression(tau), ylab=expression(f[k](tau)), main='Common Factor Loading Curves', cex.axis=1.5, cex.lab=2, cex.main=2)
abline(h=0, col='gray', lwd=4)
#for(k in 1:K.hmm.sv){
#dci = HPDinterval(as.mcmc(postD[,,k]))
#polygon(c(taugrid, rev(taugrid)), c(Phit%*%dci[,2], rev(Phit%*%dci[,1])), col='grey', border=NA)
#lines(taugrid, Phit%*%dci[,2], col=k, lty=k) # border for HPD intervals
#lines(taugrid, Phit%*%dci[,1], col=k, lty=k)
#}
for(k in 1:K.hmm.sv) lines(taugrid, Phit%*%d[,k], type='l', lwd=10, col=k, lty=k)
legend('topright', paste('k =', 1:K.hmm.sv), col=1:K.hmm.sv, lty=1:K.hmm.sv, lwd=10, cex=2)
