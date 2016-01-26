rm(list=ls(all=TRUE)) 

set.seed(14850) # to reproduce (most of) the results from the paper

#########################################################################################################################

# MCMC parameters:
nsims =  2e4 		# Total number of simulations
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
library(dlm)

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

Model = SSModel(Y~-1+SSMcustom(Z = array(0, c(ncol(Y), nrow(Gt))), T = Gt, Q = Wt, P1 = diag(10^4, nrow(Gt))))

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
      resBeta.c = Beta[, bc.inds] - S[,bc.inds]*Beta[,1:K]*matrix(rep(gammaSlopes[bc.inds], T), nrow=T, byrow=TRUE);
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
postCount = postCountAll[-(1:burnin),,]
postCountAgg = postCountAggAll[-(1:burnin),,]

# Fix important parameters at their posterior means:
Beta = colMeans(postBeta); d = colMeans(postD)
Et = colMeans(postEt); ht = colMeans(postHt)
gammaSlopes = colMeans(postGammaSlopes)

#########################################################################################################################

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

#########################################################################################################################

# Diagnostics for the FLCs:  trace plots and effective size 
# Examine quantiles of allTaus for fewer plots
# FLC trace plot for paper
quants = c(.02,.25, .5, .75); tau.inds = floor(quantile(1:m, quants))	# equally spaced
dev.new(); par(mfrow=c(2,2))
for(k in 1:K.hmm.sv) { 
  pd = as.mcmc(postDAll[,,k]%*%t(splineInfo$Phi[tau.inds,]))
  j=1; traceplot(pd[,j], cex.lab=2, cex.axis=1.5, cex.main=2, main = paste('f_k(tau): k = ',k, ', ', quants[j]*100, 'nd quantile of tau',  sep=''))
  abline(v=burnin, col='gray', lwd=2)
  for(j in 2:length(tau.inds)) {traceplot(pd[,j],cex.lab=2, cex.axis=1.5, cex.main=2, main = paste('f_k(tau): k = ',k, ', ', quants[j]*100, 'th quantile of tau',  sep='')); 	abline(v=burnin, col='gray', lwd=4)}
}

# FLC table (exclude burn-in):
effSize = matrix(0, nrow=K, ncol=length(tau.inds)); rownames(effSize) = paste('k =', 1:K); colnames(effSize) = round(allTaus0[tau.inds],2)
for(k in 1:K) { 
  pd = as.mcmc(postD[,,k]%*%t(splineInfo$Phi[tau.inds,])); colnames(pd) = paste('f_k(tau): k = ',k, ', ', quants*100, 'th quantile of tau',  sep='')
  colnames(pd)[1] =paste('f_k(tau): k = ',k, ', 2nd quantile of tau',  sep='')
  effSize[k,] = effectiveSize(pd) 
}
print(effSize/(nsims-burnin))	# column names give the maturities

# Gamma trace plot for paper
effSize = matrix(0, nrow=K.hmm.sv, ncol=C-1)
dev.new();par(mfrow=c(3,4))
for(c in 2:C){
  for(k in 1:K.hmm.sv){
    traceplot(as.mcmc(postGammaSlopesAll[,(c-1)*K + k]), cex.lab=2, cex.axis=1.5, cex.main=2, main = paste(cnames[c], 'slope, k =', k))
    abline(v=burnin, col='gray', lwd=4)
    effSize[k,c-1] = effectiveSize(as.mcmc(postGammaSlopes[,(c-1)*K + k]))
  }
}
# Gamma table (exclude burnin):
print(effSize/(nsims-burnin))

#########################################################################################################################

# DIC computations:
Dbar = mean(dev) # posterior mean of dev 
Dhat = 0
for(c in 1:C) {
  bc.inds = betaInds[c]:(betaInds[c+1]-1);	yc.inds = outcomeInds[c]:(outcomeInds[c+1]-1)
  
  # As before, but now using the posterior means	
  mu.c = tcrossprod(Beta[,bc.inds], splineInfo$Phi[match(tau[yc.inds], allTaus),]%*%d)
  Dhat = Dhat + 1/Et[c]*sum((Y[, yc.inds]-mu.c)^2,na.rm=TRUE) + sum(!is.na(Y[, yc.inds]))*log(2*pi*Et[c]) 
}
pD = Dbar - Dhat 		# effective number of parameters
pV = 1/2*var(dev)		# alternative estimate of the effective number of parameters
DIC = Dbar + pD	 	# DIC

#########################################################################################################################

# Diagnostics for the Betas: trace plots effective size
# Examine quantiles of 1:T for fewer plots
quants = c(.15, 0.3, .45, .6, .75, .9)
t.inds = floor(quantile(1:T, quants))	# equally spaced

# trace plot 
dev.new();par(mfrow=c(4,4))
for(c in 1:C){
  for(k in 1:K.hmm.sv) { 
    t.date = sample(1:T,1)
    pb = as.mcmc(postBetaAll[,t.date,(c-1)*K + k])
    traceplot(pb, cex.lab=2, cex.axis=1.5, cex.main=2, main = paste(cnames[c], ': Beta, k = ', k, ', ', dates[t.date], sep=''))
    abline(v=burnin, col='gray', lwd=4)
  }
}

# table
effSize = matrix(0, nrow=C*K, ncol=length(t.inds)); colnames(effSize) = paste(dates[t.inds])
tempnames = NULL
for(c in 1:C){
  for(k in 1:K) { 
    t.inds = sample(1:T,3)
    pb = as.mcmc(postBeta[,t.inds,(c-1)*K + k])
    colnames(pb) = paste('beta_{k,t}^{(c)}: k = ',k, ', ', ' c = ',c, ', t = ', dates[t.inds],  sep='')
    effSize[(c-1)*K + k, ] = effectiveSize(pb)
    tempnames = c(tempnames, paste('k = ',k, ', ', ' c = ',c, sep=''))
  }
}
rownames(effSize) = tempnames; print(effSize/(nsims-burnin))	

#########################################################################################################################
probS = colMeans(postCount); probSagg = colMeans(postCountAgg);  probS = rbind(0, probS); probSagg = rbind(0, probSagg)
recDates = c(as.Date("2010-08-06"), as.Date("2015-08-25")) #http://www.nber.org/cycles.htm
inds = (dates>=recDates[1])*(dates<=recDates[2])==1
dev.new(); par(mfrow=c(2,3), mai=c(1 ,1, .5,.5))
for(k in 1:K.hmm.sv){
  c=2; 	ck.inds = (c-1)*K + k
  
  plot(dates[inds], probS[inds, ck.inds], ylim=c(0,1), xaxt='n', lwd=6, col=2, xlab='Dates', ylab='Probability of Outlier', main=paste('k =', k),cex.lab=3, cex.axis=2, cex.main=4, tck = -.01)
  axis.Date(side = 1, dates[inds], at=seq(dates[inds][1], dates[inds][length(dates[inds])], by='months'),format = "%m/%y", cex.axis=2)
  abline(h=0.5, lwd=3, col='gray')
}
plot(dates[inds], probSagg[inds, 1], ylim=c(0,1), xaxt='n', lwd=6, col=2, xlab='Dates', ylab='Probability of Outlier', main='Aggregate',cex.lab=3, cex.axis=2, cex.main=4, tck = -.01)
axis.Date(side = 1, dates[inds], at=seq(dates[inds][1], dates[inds][length(dates[inds])], by='months'),format = "%m/%y", cex.axis=2)
abline(h=0.5, lwd=3, col='gray')
lines(dates[inds],probSagg[inds, 1], type='h', lwd=10, col=2)
lines(dates[inds],probSagg[inds, 2], type='h', lwd=10, col=3, lty=2)
lines(dates[inds],probSagg[inds, 3], type='h', lwd=10, col=4, lty=4)
legend('topleft', cnames[-1], lwd=10, col=2:C, cex=2 , lty=c(1,2,4))


#########################################################################################################################

# Plot the volatilities:
dev.new(); par(mfrow=c(C,K.hmm.sv), mai=c(.5,.4,.4,.1))
for(c in 1:C){
  for(k in 1:K.hmm.sv){
    color = 1 # or k
    plot(dates, exp(ht[, (c-1)*K + k]), tck=.04,  cex.lab=2, cex.axis=2.5, cex.main=4, type='l', main=paste(cnames[c], ': k = ',k,sep=''), xlab='', ylab='', ylim=quantile(exp(postHt[,,(c-1)*K + k]), probs=c(.001, .9995)))
    Bci = HPDinterval(as.mcmc(exp(postHt[,,(c-1)*K + k])))
    polygon(c(dates, rev(dates)), c(Bci[,2], rev(Bci[,1])), col='grey', border=NA)
    lines(dates, Bci[,1], col=color); lines(dates, Bci[,2], col=color) # border for HPD intervals
    lines(dates,  exp(ht[, (c-1)*K + k]), type='l', col=color, lwd=6)
  }
}

#########################################################################################################################

# Also plot the volatilities together for comparison:
dev.new(); par(mfrow=c(2,2))
for(k in 1:K.hmm.sv){
  plot(dates, exp(ht[, (c-1)*K + k]), xaxt='n',cex.lab=3, cex.axis=2, cex.main=5, type='l', main=paste('k = ',k,sep=''), xlab='Dates', ylab='', ylim=quantile(exp(postHt[,,((1:C)-1)*K + k]), probs=c(.0001, .9999)))
  axis.Date(side = 1, dates, at=seq(dates[1], dates[length(dates)], by='quarter'),format = "%m/%y", cex.axis=2)
  
  for(c in 1:C){
    color = c 
    Bci = HPDinterval(as.mcmc(exp(postHt[,,(c-1)*K + k])))
    polygon(c(dates, rev(dates)), c(Bci[,2], rev(Bci[,1])), col='lightgrey', border=NA)
    lines(dates, Bci[,1], col=color); lines(dates, Bci[,2], col=color, lwd=4, lty=5) # border for HPD intervals
    #lines(dates,  exp(ht[, (c-1)*K + k]), type='l', col=color, lwd=4)
  }
  for(c in C:1)lines(dates,  exp(ht[, (c-1)*K + k]), type='l', col=c, lwd=8)
}
legend('topleft', cnames, lwd=8, col=1:C, cex=1.5)

#########################################################################################################################

# Plot the Beta's:
dev.new(); par(mfrow=c(C,K.hmm.sv))
for(c in 1:C){
  for(k in 1:K.hmm.sv){
    plot(dates, Beta[, (c-1)*K + k], type='l', main=paste('Beta: ', cnames[c], ', k = ',k,sep=''), xlab='Dates', ylab='Beta')
    
    #Bci = HPDinterval(as.mcmc(postBeta[,,(c-1)*K + k]))
    #polygon(c(dates, rev(dates)), c(Bci[,2], rev(Bci[,1])), col='grey', border=NA)
    
    lines(dates,  Beta[, (c-1)*K + k], type='l', col=(k), lwd=3)
  }
}
#########################################################################################################################

