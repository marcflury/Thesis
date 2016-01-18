#####################################################################################################
# mfdlm() samples (once) the factors Beta AND the factor loading curves f, with the following inputs:
# Y: the data observation matrix
# T x TAU, NAs allowed; T is number of time points, TAU is the number of the union of observation points
# Assume t = 1,...,T are equally spaced (user should input NAs for completely missing observations at time t)
# tau: all of the observation points for all outcomes
# TAU-dimensional vector  
# Sorted within each outcome and concatenated from outcome to outcome
# Beta: the matrix of factors (from the previous iteration):
# T x p matrix
# where p is the number of elements in \beta_t (usually p = C*K)
# Et: the observation error variance matrix; since we assume independence, 
# C-dimensional vector  
# Gt: evolution matrix for the Beta's, either
# (1) p x p x 1 array, or
# (2) p x p x T array
# Wt: the evolution error variance matrix, either
# (1) p x p x 1 array, or
# (2) p x p x T array
# Model: the SSModel() object from the previous iteration
# Note: we update corresponding params w/in mfdlmBeta()
# d: the matrix of (common) basis coefficients for the FLCs (from the previous iteration):
# (M+4) x K matrix, where M is number of interior knots
# splineInfo: (named) list--obtained from getSplineInfo(tau)--consisting of the following:
# (1) a = min(tau), the lower endpoint for which f is evaluated
# (2) b = max(tau), the upper endpoint for which f is evaluated
# (3) intKnots, the location of interior knots
# (4) Phi, the (# ! TAU) x (M+4) matrix of spline basis functions evaluated at the unique tau's
# (5) Jkl, the (M+4) x (M+4) matrix of integrals for the inner product between dk and dl, i.e.
# Jkl = \int \phi_k(\tau) \phi_l(\tau)' d\tau =  \int \phi(\tau) \phi(\tau)' d\tau
# (6) basisPhi, a function to evaluate the basis for any tau in [a,b]
# NOTE: d_k is transformed so that the first 2 elements are the linear components and the remaining (M+4)-2 are nonlinear, with prior
# d_k ~ N(0, D), and D = diag(10^8, 10^8, 1/lambda_k, ..., 1/lambda_k)
# Phi and Jkl are transformed to match this parametrization
# Example: for VAR(1) model on Beta with K factors and SV model on Wt,
# p = K*C
# Gt = G is (K*C x K*C x 1) array
# Wt is (K*C x K*C x T) array
# for univariate SV, each Wt[,,t] is diagonal, t = 1,...,T
#####################################################################################################
mfdlm = function(Y, tau, Beta, Et, Gt, Wt, Model, d, splineInfo, lambda){
  # Sample d_k and lambda_k, k = 1,...,K, common to all outcomes:
  Fall = mfdlmF(Beta, Y, Et, tau, d, splineInfo, lambda)
  d = Fall$d # for updating Beta
  
  # Sample \beta_t, t = 1,...,T, which contains all C*K factors
  Beta = mfdlmBeta(Y, Et, Gt, Wt, Model, tau, d, splineInfo)
  
  # return Beta, d, lambda
  list(Beta=Beta, d=d, lambda = Fall$lambda)
}
#####################################################################################################
#####################################################################################################


#####################################################################################################
# mfdlmBeta() is a helper function for mfdlm() and samples the Beta's
# Simple (non-hierarchical) state space/DLM structure for Beta (T x C*K)
# One could modify this function to incorporate:
# Covariates or linear transformations of Beta:
# Sample the untransformed "Beta" (adjust the observation equation accordingly)
# Transform the resulting sample
# More sophisticated Et (correlation wrt c)
# Antithetic sampling of Beta
# Multiple lags, e.g. VAR(2)
#####################################################################################################
mfdlmBeta = function(Y, Et, Gt, Wt, Model, tau, d, splineInfo){
  # All observation points across time and outcome:
  allTaus = sort(unique(tau))
  
  # For outcome-specific subset, use yc.inds = outcomeInds[c]:(outcomeInds[c+1]-1)
  outcomeInds = c(1, which(tau[-1] < tau[-length(tau)])+1, length(tau)+1)
  
  # Define these locally	
  C = length(outcomeInds) - 1
  K = ncol(d)
  
  # F evaluates the FLCs at the outcome-specific observation points
  F = array(0, c(length(tau), nrow(Gt))); inds = 1:K
  
  # Cycle through the outcomes c = 1,...,C:
  for(c in 1:C){
    # Indices of allTaus that correspond to c-specific tau's
    tauCinds = match(tau[outcomeInds[c]:(outcomeInds[c+1]-1)], allTaus) 
    
    # Subsets: Y_t^{(c)}(\tau_1), ..., Y_t^{(c)}(\tau_m)
    yc.inds = outcomeInds[c]:(outcomeInds[c+1]-1) 	
    
    F[yc.inds, inds] = splineInfo$Phi[tauCinds,]%*%d
    
    # Update the column indexing for appropriate block structure
    if(K==1){inds = inds + 1} else {inds[1] = inds[length(inds)] + 1; inds[length(inds)] = inds[length(inds)] + K; inds = inds[1]:inds[length(inds)]}
  }
  
  # Repeat Et[c] for each outcome-specific number of tau's
  Model$H[,,1] = diag(rep(Et, diff(outcomeInds)))
  
  Model$Z[,,1] = F  # F is a matrix of correct dimensions
  Model$T = Gt 	# Gt is an array of correct dimensions
  Model$Q = Wt	# same
  
  # Check for errors
  if(!is.SSModel(Model)) stop("Error: Model has incorrect dimensions")
  
  # Run the sampler
  simulateSSM(Model, "states", nsim = 1, antithetics=FALSE, filtered=FALSE)[,,1]
  
  # Could instead use a "streamlined" version of the KFAS function, with redundant computations removed (Note: no checks for errors!)
  # Must first define modelInfo = getModelInfo(Model), which only needs to be computed (and stored) once, before running MCMC
  #simulateSSM2(Model, modelInfo)[,,1]
}
#####################################################################################################
#####################################################################################################


#####################################################################################################
# mfdlmF() is a helper function for mfdlm() and samples common FLCs and the smoothing parameters lambda
# With minor adjustments, we could incorporate:
# Semiparametric FLCs
# Independent, non-common FLCs
# More sophisticated Et
#####################################################################################################
mfdlmF = function(Beta, Y, Et, tau, d, splineInfo, lambda, orderLambdas = TRUE){
  outcomeInds = c(1, which(tau[-1] < tau[-length(tau)])+1, length(tau)+1)
  allTaus = sort(unique(tau))
  
  Phi = splineInfo$Phi; Jkl = splineInfo$Jkl
  
  # Define these locally	
  C = length(outcomeInds) - 1 
  K = ncol(Beta)/C 
  
  # For outcome-specific Beta indices, subset with bc.inds = betaInds[c]:(betaInds[c+1]-1)
  betaInds = seq(1, C*(K+1), by=K)
  
  for(k in sample(1:K)){
    # Unconstrained posterior is N(solve(Bksum)%*%bksum, solve(Bksum))
    # Adjust bksum for orthogonality constraints
    # Then normalize the posterior samples
    bksum = 0 # for summing over c = 1,...,C
    Bksum = 0 # same
    
    for(c in 1:C){
      # Indices of the c-specific tau's:
      tauCinds = match(tau[outcomeInds[c]:(outcomeInds[c+1]-1)], allTaus) 
      
      # Subsets: Y_t^{(c)}(\tau_1), ..., Y_t^{(c)}(\tau_m)
      yc.inds = outcomeInds[c]:(outcomeInds[c+1]-1) 	
      
      # Subsets: \beta_{1,t}^{(c)}, ..., \beta_{K,t}^{(c)}
      bc.inds = betaInds[c]:(betaInds[c+1]-1)		
      
      # The NAs in Y, combined w/ na.rm=TRUE, will only sum over non-missing data for bksum
      bksum = bksum + crossprod(Phi[tauCinds,], colSums((1/Et[c]*(Y[,yc.inds] - tcrossprod(Beta[,bc.inds[-k]], Phi[tauCinds,]%*%d[,-k])))*Beta[,bc.inds[k]], na.rm=TRUE))
      
      # Obtain outcome-specific complete cases:
      cc = complete.cases(Y[,yc.inds]) 
      
      # Sum over complete cases, but only look at outcome-specific tau's 
      Bksum  = Bksum + crossprod(Phi[tauCinds,])*sum(Beta[cc, bc.inds[k]]^2)/Et[c]
      
      # Only loop over incomplete cases:
      for(i in which(!cc)) Bksum  = Bksum + Beta[i,bc.inds[k]]^2/Et[c]*crossprod(Phi[tauCinds,][!is.na(Y[i,yc.inds]),])
    } 
    
    # Sample \lambda_k (precisions) using uniform priors on standard deviations, sd_k = lambda_k^-1/2, for identifiability:
    # 0 < sd_1 < sd_2 < ... < sd_K < 10^4 
    # or equivalently, Inf > \lambda_1 > \lambda_2 > ... > \lambda_K > 10^-8
    
    # Shape and rate parameters needed to enforce ordering:
    #shape0 = (ncol(Phi) - 2)/2 + 0.001; rate0 = crossprod(d[-(1:2),k])/2 + 0.001 # for Gamma(0.001, 0.001) prior
    shape0 = (ncol(Phi) + 1)/2; rate0 = crossprod(d[-(1:2),k])/2 			# for uniform prior on sd_k
    
    # Lower and upper bounds, w/ ordering constraints (if specified):
    lam.l = 10^-8; lam.u = Inf; if(orderLambdas){if(k != K) lam.l = lambda[k+1]; if(k != 1) lam.u = lambda[k-1]}
    u = runif(1, min = pgamma(lam.l, shape=shape0, rate=rate0), max = pgamma(lam.u, shape=shape0, rate=rate0))
    lambda[k] = qgamma(u, shape=shape0, rate=rate0)
    
    # Prior precision matrix implied by these \lambda_k's:
    priorPrec = diag(c(rep(10^-8, 2), rep(lambda[k], (ncol(Phi) - 2))))
    
    # Cholesky decomposition is useful for efficient sampling below
    cholFac = chol(Bksum + priorPrec)
    
    # Joint orthogonality constraints:
    Lcon = Jkl%*%d[,-k]  	
    
    # Sample the unconstrained vector:
    dk = backsolve(cholFac,forwardsolve(t(cholFac),bksum) + rnorm(length(bksum)))		
    
    # Compute the shift:
    BLcon = backsolve(cholFac,forwardsolve(t(cholFac),Lcon))
    
    # Incorporate the constraint:
    dk = dk - BLcon%*%chol2inv(chol(crossprod(Lcon, BLcon)))%*%crossprod(Lcon,dk)
    
    # Compute the norm for normalization
    Dnorm = sqrt(as.numeric(crossprod(dk,Jkl)%*%dk))
    d[,k] = dk/Dnorm
    
    # Also rescale Beta 
    Beta[, seq(from = k, to = K*C, by=K)] = Beta[, seq(from = k, to = K*C, by=K)]*Dnorm
  }
  list(d=d, lambda=lambda, Beta=Beta) 
}
#####################################################################################################
#####################################################################################################


#####################################################################################################
# initParams() initializes the main parameters in the model:
# Beta0 via an SVD of the (completed) data matrices
# Et0 as the conditional MLE, using the SVD
# d0 and lambda0 via mfdlmFinit()
# splineInfo via getSplineInfo()
# Also, specify whether to initialize using all observation points (slow) or a subset
# Note: we also run 10 initial simulations of Beta0 and d0 without the ordering constraints on lambda 
#####################################################################################################
initParams = function(Y, tau, K=NULL, tolCPV = 0.99,  useAllTaus = FALSE){
  allTaus = sort(unique(tau))
  
  # For outcome-specific subset, use outcomeInds[c]:(outcomeInds[c+1]-1)
  outcomeInds = c(1, which(tau[-1] < tau[-length(tau)])+1, length(tau)+1)
  
  # also can infer the number of outcomes
  C = length(outcomeInds) - 1	
  
  # and the number of times:
  T = nrow(Y)
  
  # relevant spline info:
  splineInfo = getSplineInfo(tau)
  
  # Specificy the tau values to interpolate for the Beta and F initialization (SVD)
  if(useAllTaus){
    useTaus = allTaus
  } else useTaus = c(splineInfo$a, splineInfo$intKnots, splineInfo$b)
  
  # For initialization:
  Y0 = array(0, c(C*T, length(useTaus)))
  
  for(c in 1:C){
    # just look at outcome c:
    Y0c = Y[, outcomeInds[c]:(outcomeInds[c+1]-1)]
    
    # if Y_t is partially missing, smooth across tau (for each such t):
    notCC = which(rowSums(!is.na(Y0c)) != 0 )
    Y0c[notCC,] = t(apply(Y0c[notCC,], 1, function(x) splinefun(tau[outcomeInds[c]:(outcomeInds[c+1]-1)], x, method='natural')(tau[outcomeInds[c]:(outcomeInds[c+1]-1)])))
    
    # if Y_t is completely missing, smooth across t for each observation point
    Y0c = apply(Y0c, 2, function(x){splinefun(1:T, x, method='natural')(1:T)})
    
    # Y0 is a stacked data matrix, with C*T rows and length(allTaus) columns
    Y0[(1:T) + (c-1)*T, ] = t(apply(Y0c, 1, function(x) splinefun(tau[outcomeInds[c]:(outcomeInds[c+1]-1)], x, method='natural')(useTaus)))
  }
  
  # Compute SVD of the (completed) data matrix:
  singVal = svd(Y0)
  
  # Cumulative sums of the s^2 proportions
  # More reasonable when Y has been centered
  cpv = cumsum(singVal$d^2/sum(singVal$d^2))
  print(t(matrix(cpv,dimnames=list(paste('k =',1:length(singVal$d))))))	
  
  # If K is unspecified, select based on cpv 
  if(is.null(K)) K = max(2, which(cpv >= tolCPV)[1])
  
  # For outcome-specific Beta indices, subset with betaInds[c]:(betaInds[c+1]-1)
  betaInds = seq(1, C*(K+1), by=K)
  F0 = vector("list", C); Beta0 = array(0,c(T, C*K)); Et0 = numeric(C); 
  
  for(c in 1:C){
    # If we use all tau's, then just subset for each outcome; otherwise, smooth the SVD estimates and interpolate the outcome-specific tau's
    if(useAllTaus){
      F0[[c]] = as.matrix(singVal$v[match(tau[outcomeInds[c]:(outcomeInds[c+1]-1)], allTaus),1:K])
    } else F0[[c]] = apply(as.matrix(singVal$v[,1:K]), 2, function(x){splinefun(useTaus, x, method='natural')(tau[outcomeInds[c]:(outcomeInds[c+1]-1)])})
    
    Beta0[,betaInds[c]:(betaInds[c+1]-1)] = (singVal$u%*%diag(singVal$d))[(1:T) + (c-1)*T , 1:K]
    
    # Estimate the observation-level error variance:
    Et0[c] = 1/sum(!is.na(Y[,outcomeInds[c]:(outcomeInds[c+1]-1)])) * sum((Y[, outcomeInds[c]:(outcomeInds[c+1]-1)] - tcrossprod(Beta0[,betaInds[c]:(betaInds[c+1]-1)], F0[[c]]))^2, na.rm=TRUE)
  }
  # Initialize the common basis coefficients
  Fall = mfdlmFinit(Beta0, Y, Et0, tau, F0, splineInfo)
  d0 = Fall$d; lambda0 = Fall$lambda; Beta0 = Fall$Beta 
  
  # To obtain ordering, run a few simple simulations:
  dArray = array(diag(C*K), c(C*K, C*K, 1)) # diagonal 
  Model0 = SSModel(Y~-1+SSMcustom(Z = array(0, c(ncol(Y), C*K)), T = dArray, Q = dArray, P1 = diag(10^4, C*K)))
  for(nsi in 1:10){
    samples = mfdlmF(Beta0, Y, Et0, tau, d0, splineInfo, lambda0, orderLambdas = FALSE)
    d0 = samples$d; lambda0 = samples$lambda;
    Beta0 = mfdlmBeta(Y, Et0, dArray, dArray, Model0, tau, d0, splineInfo)
    print(paste('Running preliminary simulations: ', nsi,'/10', sep=''))
  }
  
  # Order the k components based on the smoothing parameters, lambda:
  adjOrder = order(lambda0, decreasing = TRUE)
  lambda0 = lambda0[adjOrder]; d0 = d0[,adjOrder]; 
  for(c in 1:C){ bc.inds = betaInds[c]:(betaInds[c+1]-1)
  Beta0[,bc.inds] = Beta0[,bc.inds][, adjOrder]}
  
  # Make sure the k = 1 curve is positive (usually an overall level or intercept)
  # Might as well initialize all curves to have positive sums
  for(k in 1:K){ if(sum(splineInfo$Phi%*%d0[,k]) < 0){ d0[,k] = - d0[,k]
  Beta0[,seq(from = k, to = K*C, by=K)] = - Beta0[,seq(from = k, to = K*C, by=K)]}}
  
  list(Beta0=Beta0, d0=d0, Et0=Et0, lambda0=lambda0, splineInfo=splineInfo, K=K)
}
#####################################################################################################
#####################################################################################################


#####################################################################################################
# mfdlmFinit() initializes the common FLCs and the smoothing parameters 
# Based on the SVD of the data matrix 
# Nearly identical to mfldlmF()
# Here, we use sequential orthogonality constraints
#####################################################################################################
mfdlmFinit = function(Beta, Y, Et, tau, F0, splineInfo){
  outcomeInds = c(1, which(tau[-1] < tau[-length(tau)])+1, length(tau)+1)
  allTaus = sort(unique(tau))
  
  Phi = splineInfo$Phi; Jkl = splineInfo$Jkl
  
  C = length(outcomeInds) - 1 
  K = ncol(Beta)/C 
  
  # For outcome-specific Beta indices, subset with betaInds[c]:(betaInds[c+1]-1)
  betaInds = seq(1, C*(K+1), by=K)
  
  # Start w/ a small value for smoothness
  lambda = rep(10^-8, K) 
  
  d = array(0, c(ncol(Phi), K))
  
  for(k in 1:K){
    bksum = 0; Bksum = 0 
    for(c in 1:C){
      tauCinds = match(tau[outcomeInds[c]:(outcomeInds[c+1]-1)], allTaus)
      yc.inds = outcomeInds[c]:(outcomeInds[c+1]-1)
      bc.inds = betaInds[c]:(betaInds[c+1]-1)		
      bksum = bksum + crossprod(Phi[tauCinds,], colSums((1/Et[c]*(Y[,yc.inds] - tcrossprod(Beta[,bc.inds[-k]], F0[[c]][,-k])))*Beta[,bc.inds[k]], na.rm=TRUE))
      cc = complete.cases(Y[,yc.inds]) 
      Bksum  = Bksum + crossprod(Phi[tauCinds,])*sum(Beta[cc,(betaInds[c]:(betaInds[c+1]-1))[k]]^2)/Et[c]
      for(i in which(!cc)) Bksum  = Bksum + Beta[i,bc.inds[k]]^2/Et[c]*crossprod(Phi[tauCinds,][!is.na(Y[i,yc.inds]),])
    } 
    priorPrec = diag(c(rep(10^-8, 2), rep(lambda[k], (ncol(Phi) - 2))))
    cholFac = chol(Bksum + priorPrec)
    dk = backsolve(cholFac,forwardsolve(t(cholFac),bksum))	
    # Use sequential orthogonality for initialization:	
    if(k > 1){
      Lcon = Jkl%*%d[,1:(k-1)]  	
      BLcon = backsolve(cholFac,forwardsolve(t(cholFac),Lcon))
      dk = dk - BLcon%*%chol2inv(chol(crossprod(Lcon, BLcon)))%*%crossprod(Lcon,dk)
    }
    Dnorm = sqrt(as.numeric(crossprod(dk,Jkl)%*%dk))
    d[,k] = dk/Dnorm
    Beta[, seq(from = k, to = K*C, by=K)] = Beta[, seq(from = k, to = K*C, by=K)]*Dnorm
    
    # Conditional MLE of lambda:
    lambda[k] = (ncol(Phi) - 2)/crossprod(d[-(1:2),k])
    
    for(c in 1:C){ tauCinds = match(tau[outcomeInds[c]:(outcomeInds[c+1]-1)], allTaus)
    F0[[c]][,k] = Phi[tauCinds,]%*%d[,k]}
  }
  list(d=d, lambda=lambda, Beta=Beta) 
}
#####################################################################################################
#####################################################################################################


#####################################################################################################
# initSV() initializes the stochastic volatility model (Yield curve application)
# resBeta: (T-1) x (C*K) matrix of residuals from the Beta-level equation
# Should be approximately mean zero
# C: number of outcomes
#####################################################################################################
initSV = function(resBeta, C){
  # Define locally
  K = ncol(resBeta)/C
  
  # The evolution error variance to be passed to KFAS (state space sampler):
  Wt = array(0,c(C*K,C*K,T))
  
  # Log-volatilities, plus an offset
  # First row is t=0 (update later)
  ht = rbind(1, log(resBeta^2 + 10^-8))
  
  # SV parameters:
  svMu = svPhi = svSigma =  numeric(C*K) 
  
  for(c in 1:C){ 
    for(k in 1:K){
      # Useful index to access outcome c, factor k, within the vector(s)
      ck.inds = (c-1)*K + k
      
      # SV model for k = 1,...,K
      
      # AR(1) model for initialization:
      modh = arima(ht[-1,ck.inds], c(1,0,0))
      
      # The AR(1) coefficient  
      svPhi[ck.inds] = coef(modh)[1]
      
      # The unconditional mean 
      svMu[ck.inds] = coef(modh)[2]
      
      # The error SD
      svSigma[ck.inds] = sqrt(modh$sigma2)
      
      # Also initial h_1 at the unconditional mean:
      ht[1,ck.inds] = svMu[ck.inds]
      
      # For the KFAS Model:
      Wt[ck.inds, ck.inds, ] = exp(ht[,ck.inds])
    }			
  }
  list(ht = ht, svMu = svMu, svPhi = svPhi, svSigma = svSigma, Wt = Wt)
}
#####################################################################################################
#####################################################################################################


#####################################################################################################
# sampleSV() samples the stochastic volatility model for fixed c, across k (Yield curve application)
# resBeta.c: (T-1) x K matrix of residuals from the Beta-level for outcome c
# ht.c: the T x K matrix of log-volatilities for outcome c
# svMu.c: the K-dimensional vector of the unconditional means of the log-volatilities for outcome c
# svPhi.c: the K-dimensional vector of the AR(1) coefficients of the log-volatilities for outcome c
# svSigma.c: the K-dimensional vector of the error variances of the log-volatilities for outcome c
# svOffset: an offset we can add to the squared residuals (before taking logs) for computational stability 
# then adjust the implied variances after sampling (Note: other SV parameters are not rescaled, but do not appear elsewhere in the model)
#####################################################################################################
sampleSV = function(resBeta.c,  ht.c, svMu.c, svPhi.c, svSigma.c, svOffset = 0){
  # Define locally
  K = ncol(resBeta.c)
  
  # To update the Beta-level error variance for KFAS
  Wt.c = array(0, c(K, K, nrow(ht.c)))
  
  for(k in 1:K){
    resBeta.ck = resBeta.c[,k]
    
    # Sample the SV model
    
    # Update the "start parameters" to most recent MCMC sample for outcome c, factor k
    startpara = list(mu = svMu.c[k], phi = svPhi.c[k], sigma = svSigma.c[k])
    
    # Include an offset for computational stability (then rescale after sampling):
    resBeta.ck = sqrt(1+svOffset/resBeta.ck^2)*resBeta.ck
    
    sv = svsample2(resBeta.ck, draws = 1, priormu=c(0,100), priorphi = c(5,1.5), priorsigma=1, startpara=startpara , startlatent = ht.c[-1, k])
    
    # Store the results
    svMu.c[k] = para(sv)[1]; svPhi.c[k] = para(sv)[2]
    svSigma.c[k] = para(sv)[3]; ht.c[,k] = c(sv$latent0, sv$latent) 
    Wt.c[k, k, ] = exp(ht.c[,k])
    Wt.c[k, k, -1] = Wt.c[k, k, -1]/(1+svOffset/resBeta.ck^2)
  }
  list(ht.c = ht.c, svMu.c = svMu.c, svPhi.c = svPhi.c, svSigma.c = svSigma.c, Wt.c = Wt.c)
}
#####################################################################################################
#####################################################################################################



#####################################################################################################
# sampleEt() samples the observation-level error variance 
# Assumes diagonality, so it returns a scalar (for outcome c)
# Prior on the precision is Gamma(gamma1, gamma2)
# Inputs:
# Y.c: the observation matrix Y restricted to outcome c
# mu.c: the MFDLM estimate of Y.c (conditional mean)
# gamma1, gamma2: the priors
#####################################################################################################
sampleEt = function(Y.c, mu.c, gamma1 = 0.001, gamma2 = 0.001){
  1/rgamma(1, shape = (gamma1 + sum(!is.na(Y.c))/2), rate = (gamma2 +sum((Y.c - mu.c)^2, na.rm=TRUE)/2))
}
#####################################################################################################
#####################################################################################################



#####################################################################################################
# sampleWk() samples the evolution-level error variance (LFP application)
# resBeta: matrix of residuals from the Beta-level equation
# Two options:
# useDiagonal = TRUE: assumes Wt is digaonal
# Components are independent w/ priors (on precisions) Gamma(0.001, 0.001)
# useDiagonal = FALSE: assumes Wt is a (permutated) block matrix, i.e.
# For each k, the error covariance matrix Wk is C x C (must then be permuted to match our definition of Beta)
# Prior on each Wk is inverse Wishart w/ identity prior precision and scale parameter C
#####################################################################################################
sampleWk = function(resBeta, Rinv = diag(1, C), rh0 = C, useDiagonal=FALSE){
  
  if(useDiagonal){
    # Assumes independent Gamma(0.001, 0.001) priors for each component
    diag(apply(resBeta, 2, function(x) 1/rgamma(n=1, shape = 0.001 + (length(x)-1)/2, rate = sum(x^2)/2 + 0.001)))
  } else {
    Wtemp = diag(C*K) # storage
    for(k in 1:K){
      # Need the correct indices for fixed k, across outcomes c=1,...,C
      c.inds = seq(from = k, to = K*C, by = K)
      
      # Wishart full conditional posterior:
      # The inversions are for K x K matrices, which should be small
      Wtemp[c.inds, c.inds] = chol2inv(chol(rwish(rh0 + nrow(resBeta), chol2inv(chol(Rinv*rh0 + crossprod(resBeta[, c.inds]))))))
    }
    Wtemp
  }
}
#####################################################################################################
#####################################################################################################


#####################################################################################################
# getSplineInfo() initializes (and transforms) the spline basis
# Uses quantile-based placement of knots for a cubic spline basis
# Enfoces a penalty on the integrated squared second derivative 
# Computes the matrix of integrals for the orthonormality constraints
# Transforms the basis, penalty, and matrix of integrals so that:
# d_k is decomposed into linear (first 2 elements) and nonlinear components
# the resulting prior for d_k is diagonal, which helps with posterior mixing, and proper
# Follows Wand and Ormerod (2008)
# Inputs:
# tau: all of the observation points for all outcomes
#####################################################################################################
getSplineInfo = function(tau){
  allTaus = sort(unique(tau)) 	# all observation points
  a = min(allTaus)        	# lower endpoint
  b = max(allTaus)        	# upper endpoint
  
  numIntKnots = 20													# number of interior knots (= M in paper)
  intKnots = quantile(allTaus,seq(0,1,length= (numIntKnots+2))[-c(1,(numIntKnots+2))]) 	# interior knots
  
  basis = create.bspline.basis(c(a,b),breaks=c(a,intKnots,b))						# spline basis
  blin = create.monomial.basis(c(a,b), nbasis=2) 								# linear basis
  
  Phi = bs(allTaus,knots=intKnots,degree=3,Boundary.knots=c(a,b),intercept=TRUE) 		# basis matrix
  Omega = eval.penalty(basis, Lfdobj=2, rng=c(a,b)) 					 		# derivative matrix
  Jkl = eval.penalty(basis, Lfdobj=0, rng=c(a,b))						 		# integral matrix
  
  # Now, transform these matrices to align with the linear/nonlinear decomposition of d_k:
  # The first two entries of d_k are the linear components, the remaining are nonlinear (see Wand and Ormerod, 2008)
  eigOmega = eigen(Omega)
  indsZ = 1:(numIntKnots+2)
  UZ = eigOmega$vectors[, indsZ] 			# Linear part
  LZ = t(t(UZ)/sqrt(eigOmega$values[indsZ])) 	# Nonlinear part
  
  # Basis matrices
  PhiL = cbind(1, allTaus)				# Linear part
  PhiN = Phi%*%LZ						# Nonlinear part
  Phi = cbind(PhiL, PhiN)
  
  # Basis functions:	
  basisPhiL = function(x){cbind(1,x)}
  basisPhiN = function(x){bs(x,knots=intKnots,degree=3,Boundary.knots=c(a,b),intercept=TRUE)%*%LZ}
  basisPhi = function(x){cbind(1,x, bs(x,knots=intKnots,degree=3,Boundary.knots=c(a,b),intercept=TRUE)%*%LZ)}
  
  # Integrals:
  #JLL = inprod(blin, blin)
  JLL =  matrix(c(b-a, (b^2 - a^2)/2, (b^2 - a^2)/2, (b^3 - a^3)/3), nrow=2) 	# Linear x Linear
  JLN = inprod(blin, basis)%*%LZ								# Linear x Nonlinear
  JNN = crossprod(LZ, Jkl)%*%LZ									# Nonlinear x Nonlinear
  
  # Combine the integrals into one matrix:
  Jkl = cbind(rbind(JLL, t(JLN)), rbind(JLN, JNN))
  
  list(a=a, b=b, intKnots=intKnots, Phi=Phi, Jkl= Jkl, basisPhi = basisPhi)
}
#####################################################################################################
#####################################################################################################



#####################################################################################################
# initHMM() initializes the HMM/CT model
# Beta: factors (already initialized)
# K.hmm.sv: how many of the K factors to model with the HMM/CT model
# the other factors will still have an AR(1)-SV(1) model
# useHMM: TRUE for HMM, FALSE for CT (CT model in paper)
# if TRUE, u00, u10, u01, u11 are prior "sample sizes" that count state transitions from 0 to 0, 1 to 0, etc.
#####################################################################################################
initHMM = function(Beta, K.hmm.sv, useHMM = FALSE, u01 = 1, u00 = 1, u10 = 1, u11 = 1){
  
  # AR(1) parameters:
  psi = apply(Beta, 2, function(x){arima(x, c(1,0,0), include.mean=FALSE)$coef}) 
  
  # Matrix of states in {0,1}
  S = matrix(0, nrow=T, ncol=C*K)
  
  if(useHMM){
    # For each outcome and factor, define
    # q01 = P(s_t = 1|s_{t-1} = 0) 
    # q10 = P(s_t = 0|s_{t-1} = 1)
    
    # Initialize with the prior means:
    q01 = rep(u01/(u01+u00), C*K) 
    q10 = rep(u10/(u10+u11), C*K) 
    q01[1:K] = q10[1:K] = NA  	# c=1 is a placeholder for simpler indexing
  }
  
  # Slope terms:
  gammaSlopes = rep(0, C*K)
  
  for(c in 2:C){ 
    for(k in 1:K){
      # Useful indices for factor k, outcome c:
      ck.inds = (c-1)*K + k
      
      if(useHMM && k<= K.hmm.sv){
        # Initialize as a random sample, using prior values:
        S[, ck.inds] = sample(c(0,1), T, replace=TRUE, prob=c(u01+u00,u10+u11))
        
        # Now actually sample the states:
        S[, ck.inds] = sampleHMMstates(Beta[, ck.inds], Beta[,k], psi[ck.inds], S[,ck.inds], q01[ck.inds], q10[ck.inds], rep(1, T))
        
        # Sample q01 and q10
        qtemp = sampleTransitProbs(S[, ck.inds])
        q01[ck.inds] = qtemp$q01.ck; q10[ck.inds] = qtemp$q10.ck
      } else S[,ck.inds] = 1
      
      # Sample the slopes
      if(k <= K.hmm.sv) gammaSlopes[ck.inds] = sampleSlopes(Beta[, ck.inds], Beta[,k], psi[ck.inds], S[,ck.inds], rep(1, T))		
    }	
  }
  
  if(useHMM){
    list(psi=psi, q01=q01, q10=q10, S=S, gammaSlopes=gammaSlopes)
  } else list(psi=psi, S=S, gammaSlopes=gammaSlopes) # still return S for simpler computations elsewhere
}
#####################################################################################################
#####################################################################################################


#####################################################################################################
# samplePsi() samples the AR(1) coefficient in the HMM/CT model
# resBeta.ck: T-dimensional vector of residuals from Beta-level for outcome c, factor k
# sigma.t2: T-dimensional vector of error variances from the Beta-level for outcome c, factor k
# Assumes prior N(0, 10^8) truncated to (-1,1) for stationarity
#####################################################################################################
samplePsi = function(resBeta.ck, sigma.t2){
  
  n = length(resBeta.ck)
  
  Vi = 1/(sum(resBeta.ck[-n]^2/sigma.t2[-1]) + 10^-8)
  u = sum(resBeta.ck[-n]*resBeta.ck[-1]/sigma.t2[-1])
  
  rtruncnorm(1,a=-1, b=1, mean= Vi*u, sd = sqrt(Vi))
}
#####################################################################################################
#####################################################################################################



#####################################################################################################
# sampleSlopes() samples the slopes in the HMM/CT model
# B.ck: Beta for outcome c (> 1), factor k
# B.1k: Beta for outcome 1, factor k
# psi.ck: AR(1) coefficient for outcome c, factor k
# St: T-dimensional vector of states for outcome c, factor k; St = 1 for CT model
# sigma.t2: T-dimensional vector of error variances from the Beta-level for outcome c, factor k
# Assumes prior N(0, 10^8)
# We could instead use a hierarchical prior to incorporate some shrinkage
#####################################################################################################
sampleSlopes = function(B.ck, B.1k, psi.ck, St, sigma.t2){
  
  n = length(B.ck)
  
  # Simple regression framework:
  Ytemp = (B.ck[-1] - psi.ck*B.ck[-n])/sqrt(sigma.t2[-1])
  Xtemp = (St[-1]*B.1k[-1] - psi.ck*St[-n]*B.1k[-n])/sqrt(sigma.t2[-1])
  
  Vi = 1/(sum(Xtemp^2) + 10^-8)
  u = sum(Xtemp*Ytemp) 
  
  rnorm(1, mean= Vi*u, sd = sqrt(Vi))
}
#####################################################################################################
#####################################################################################################


#####################################################################################################
# sampleHMMpar() samples the HMM/CT parameters
# Beta: T x (C*K) matrix of factors
# K.hmm.sv: how many factors for the HMM/CT model (other factors will still have AR(1)-SV(1))
# S: T x (C*K) matrix of states
# gammaSlopes: (C*K)-dimensional vector of slopes (w/ appropriate placement of zeros)
# psi: (C*K)-dimesional vector AR(1) coefficients
# ht: T x (C*K) matrix of log-volatilities from the Beta-level
# c: which outcome; defaults to 1
# useHMM: TRUE for HMM, FALSE for CT (CT model in paper)
# q01: transition probability from 0 to 1 in HMM model
# q10: transition probability from 1 to 0 in HMM model
#####################################################################################################
sampleHMMpar = function(Beta, K.hmm.sv, S, gammaSlopes, psi, ht, c=1, useHMM=FALSE, q01 = NULL, q10 = NULL){
  
  for(k in 1:K) {
    # Useful indices for factor k, outcome c:
    ck.inds = (c-1)*K + k
    
    if(c > 1){
      if(useHMM){
        # Sample the states
        S[, ck.inds] = sampleHMMstates(Beta[, ck.inds], gammaSlopes[ck.inds]*Beta[,k], psi[ck.inds], S[,ck.inds], q01[ck.inds], q10[ck.inds], exp(ht[, ck.inds]))
        
        # Sample q01 and q10
        qtemp = sampleTransitProbs(S[, ck.inds])
        q01[ck.inds] = qtemp$q01.ck
        q10[ck.inds] = qtemp$q10.ck
      } 
      
      # Sample the slopes
      if(k <= K.hmm.sv) gammaSlopes[ck.inds] = sampleSlopes(Beta[, ck.inds], Beta[,k], psi[ck.inds], S[,ck.inds], exp(ht[, ck.inds]))            
      
      # Sample the AR(1) coefficients:
      resBeta.ck = Beta[,ck.inds] - S[, ck.inds]*gammaSlopes[ck.inds]*Beta[,k]
      psi[ck.inds] = samplePsi(resBeta.ck, exp(ht[,ck.inds]))
      
    } else  psi[ck.inds] = samplePsi(Beta[, ck.inds], exp(ht[,ck.inds]))  
  }	
  if(useHMM){
    list(S = S, q01 = q01, q10 = q10, gammaSlopes=gammaSlopes, psi = psi)
  } else list(gammaSlopes=gammaSlopes, psi = psi)
  
}
#####################################################################################################
#####################################################################################################


#####################################################################################################
# computeGtWtHMM() computes the Gt and Wt matrices that appear in the KFAS sampler for the HMM/CT model
# Gt: evolution coefficients array, (C*K) x (C*K) x T
# Wt: evoution error variance array, (C*K) x (C*K) x T
# S: T x (C*K) matrix of states
# gammaSlopes: (C*K)-dimensional vector of slopes (w/ appropriate placement of zeros)
# psi: (C*K)-dimesional vector AR(1) coefficients
# sigma.t2: T-dimensional vector of error variances from the Beta-level for outcome c, factor k
# Note: this function is NOT optimized for the CT model, and is rather slow
#####################################################################################################
computeGtWtHMM = function(Gt, Wt, S, gammaSlopes, psi, sigma.t2){
  
  # Compute S_t*gammaSlopes*psi, and for psi with c=1 only
  S.gamma.psi_t = S*matrix(rep(gammaSlopes*psi, T), nrow=T, byrow=TRUE)
  S.gamma.psi1_t = S*matrix(rep(gammaSlopes*rep(psi[1:K], C), T), nrow=T, byrow=TRUE)
  
  # For outcome-specific Beta indices, subset with bc.inds = betaInds[c]:(betaInds[c+1]-1)
  betaInds = seq(1, C*(K+1), by=K)
  
  # i = 1:
  Wt[,,1] = diag(sigma.t2[1,])
  Gt[,, 1] = diag(psi)# diag(0,C*K)
  
  for(i in 2:T){
    #Set the usual diagonal SV components:
    Wt[1:(C*K),1:(C*K),i] = diag(sigma.t2[i,])
    
    # And the usual diagonal AR(1) model:
    Gt0 = diag(psi)
    
    for(c in 2:C){
      bc.inds = betaInds[c]:(betaInds[c+1]-1)
      
      # Update the common component based on outcome c
      Gt0[bc.inds,1:K] = diag(S.gamma.psi1_t[i, bc.inds] - S.gamma.psi_t[i-1, bc.inds])
      
      # Now update the variance for columns 1:K (and rows 1:K by symmetry)
      Wt[bc.inds, 1:K,i] = Wt[1:K, bc.inds, i] = diag(S[i,bc.inds]*gammaSlopes[bc.inds]*sigma.t2[i, 1:K])
      
      # Nested loop without the outcomes
      for(c2 in 2:C){
        bc2.inds = betaInds[c2]:(betaInds[c2+1]-1)
        # Diagonal and off-diagonal (symmetric!)
        if(all(bc.inds == bc2.inds)){
          Wt[bc.inds, bc2.inds,i] = Wt[bc.inds, bc2.inds,i] + diag(S[i,bc.inds]*S[i,bc2.inds]*gammaSlopes[bc.inds]*gammaSlopes[bc2.inds]*sigma.t2[i, 1:K])
        } else Wt[bc.inds, bc2.inds,i] = Wt[bc2.inds, bc.inds,i] = diag(S[i,bc.inds]*S[i,bc2.inds]*gammaSlopes[bc.inds]*gammaSlopes[bc2.inds]*sigma.t2[i, 1:K])
      }
    } 
    Gt[,, i] = Gt0
  }
  list(Gt = Gt, Wt=Wt)
}
#####################################################################################################
#####################################################################################################



#####################################################################################################
# computeTransitProbs() computes transition probabilities between states for the HMM model
# Inputs:
# stm1 = s_{t-1}
# st = s_t
# q01 = P(st = 1|stm1 = 0)
# q10 = P(st = 0|stm1 = 1)
#####################################################################################################
computeTransitProbs = function(stm1,st, q01, q10){
  (1-q01)^((1-st)*(1-stm1))*q10^((1-st)*stm1)*q01^(st*(1-stm1))*(1-q10)^(st*stm1)
}
#####################################################################################################
#####################################################################################################



#####################################################################################################
# sampleTransitProbs() samples transition probabilities q01 and q10
# Inputs:
# S.ck: states for all times t=1,...,T, for outcome c, factor k
# u01+1, u00+1, u10+1, u11+1 are the hyperaparameters in the Beta priors
#####################################################################################################
sampleTransitProbs = function(S.ck, u01 = 1, u00 = 10, u10 = 1, u11 = 10){
  
  T = length(S.ck)
  
  # Count the state transitions
  n01 = sum((S.ck[1:(T-1)]==0)*(S.ck[2:T]==1))
  n00 = sum((S.ck[1:(T-1)]==0)*(S.ck[2:T]==0))
  n10 = sum((S.ck[1:(T-1)]==1)*(S.ck[2:T]==0))
  n11 = sum((S.ck[1:(T-1)]==1)*(S.ck[2:T]==1))
  
  q01.ck = rbeta(1, n01 + u01+1, n00 + u00+1)
  q10.ck = rbeta(1, n10 + u10+1, n11 + u11+1)
  
  list(q01.ck = q01.ck, q10.ck = q10.ck)
}
#####################################################################################################
#####################################################################################################



#####################################################################################################
# sampleHMMstates() samples the HMM states
# B.ck: Beta for outcome c (> 1), factor k
# B.1k: Beta for outcome 1, factor k
# psi.ck: AR(1) coefficient for outcome c, factor k
# St: T-dimensional vector of states for outcome c, factor k 
# q01.ck: transition probability from 0 to 1 in HMM model for outcome c, factor k
# q10.ck: transition probability from 1 to 0 in HMM model for outcome c, factor k
# sigma.t2: T-dimensional vector of error variances from the Beta-level for outcome c, factor k
# See Albert and Chib, 1993, for more details
#####################################################################################################
sampleHMMstates = function(B.ck, B.1k, psi.ck, St, q01.ck, q10.ck, sigma.t2){
  
  n = length(B.1k)
  
  # Residuals from AR(1) model:
  resBc = B.ck[2:n] - psi.ck* B.ck[1:(n-1)] 
  resBc = c(0,resBc) # for indexing 
  
  Stemp = St; counter = 0; notDiverse=TRUE # all 0s or all 1s
  
  while(notDiverse && counter < 50){
    counter = counter + 1
    
    St = Stemp
    
    # Special case: i=n 
    i = n; 
    
    # First, consider state = 1
    St[i] = 1
    # Residual contribution from Beta_1k in Beta_ck equation
    resB1 =  B.1k[i:(i-1)]*St[i:(i-1)]; resB1 = resB1[1] - psi.ck*resB1[2] 
    ldist1 =  log(computeTransitProbs(St[i-1], 1, q01.ck, q10.ck))-.5*(resBc[i]-resB1)^2/sigma.t2[i]
    
    # Next, consider state = 0
    St[i] = 0
    resB1 =  B.1k[i:(i-1)]*St[i:(i-1)]; resB1 = resB1[1] - psi.ck*resB1[2] 
    ldist0 = log(computeTransitProbs(St[i-1], 0, q01.ck, q10.ck))-.5*(resBc[i]-resB1)^2/sigma.t2[i]
    
    # Compare 
    if((ldist0 - ldist1) < log(1/runif(1)-1)) St[i] = 1
    
    for(i in (n-1):2){
      St[i] = 1
      resB1 = B.1k[(i+1):(i-1)]*St[(i+1):(i-1)]
      resB1 = c(resB1%*%c(1, -psi.ck,0), resB1%*%c(0,1, -psi.ck))
      ldist1 = log(computeTransitProbs(1, St[i+1], q01.ck, q10.ck))+log(computeTransitProbs(St[i-1], 1, q01.ck, q10.ck))-.5*sum((resBc[(i+1):(i)] - resB1)^2/sigma.t2[(i+1):i])
      
      
      St[i] = 0
      resB1 = B.1k[(i+1):(i-1)]*St[(i+1):(i-1)]
      resB1 = c(resB1%*%c(1, -psi.ck,0), resB1%*%c(0,1, -psi.ck))
      ldist0 = log(computeTransitProbs(0, St[i+1], q01.ck, q10.ck)) + log(computeTransitProbs(St[i-1], 0, q01.ck, q10.ck))-.5*sum((resBc[(i+1):(i)] - resB1)^2/sigma.t2[(i+1):i])
      if((ldist0 - ldist1) < log(1/runif(1)-1)) St[i] = 1			
    }
    
    # Special case: i=1 (ignore stationary component for simplicity)
    i = 1
    pi0 = q10.ck/(q01.ck+q10.ck) # = P(S_t = 0) 
    
    St[i] = 1
    resB1 = B.1k[(i+1):i]*St[(i+1):i]; resB1 = resB1[1] - psi.ck*resB1[2] 
    ldist1 =  log(1-pi0) + log(computeTransitProbs(1, St[i+1], q01.ck, q10.ck))-.5*(resBc[i+1]-resB1)^2/sigma.t2[i+1]
    
    St[i] = 0
    resB1 = B.1k[(i+1):i]*St[(i+1):i]; resB1 = resB1[1] - psi.ck*resB1[2] 
    ldist0 =  log(pi0) + log(computeTransitProbs(0, St[i+1], q01.ck, q10.ck))-.5*(resBc[i+1]-resB1)^2/sigma.t2[i+1]
    if((ldist0 - ldist1) < log(1/runif(1)-1)) St[i] = 1
    
    if(!all(St==0) && !all(St==1)) notDiverse=FALSE
  } 
  if(counter==50) print('Note: did not sample S_t')
  
  # Return the states:
  St
}
#####################################################################################################
#####################################################################################################



#####################################################################################################
# simulateSSM2(): a reduced version of simulateSSM() from the KFAS package, 
# which eliminates redundant computations within the MCMC 
# object: KFAS model
# modelInfo: computed from getModelInfo(); stores the variables that only need to be computed once
#####################################################################################################
simulateSSM2 = function (object, modelInfo){
  sim.what = "states"
  nsim = 1
  p <- attr(object, "p") # sum m_t
  m <- attr(object, "m") # K*C
  k <- attr(object, "k") # K*C
  n <- attr(object, "n") # T
  tv <- attr(object, "tv") # which are time-varying?
  
  simtmp = simHelper2(object, modelInfo)
  
  sim.what <- which(c("epsilon", "eta", "disturbances", "states", "signals", "observations") == sim.what)
  simdim <- as.integer(switch(sim.what, p, k, p + k, m, p, p))
  
  
  out <- .Fortran(KFAS:::fsimgaussian, NAOK = TRUE, modelInfo$ymiss, tv, 
                  object$y, object$Z, object$H, object$T, object$R, 
                  object$Q, object$a1, object$P1, object$P1inf, modelInfo$nNonzeroP1, 
                  as.integer(nsim), simtmp$epsplus, simtmp$etaplus, 
                  simtmp$aplus1, p, n, m, k, info = as.integer(0), 
                  modelInfo$nNonzeroP1inf, object$tol, modelInfo$zeroP1inf, 
                  length(modelInfo$zeroP1inf), sim = array(0, c(simdim, n, nsim)), 
                  modelInfo$c2, sim.what, simdim, 0)
  
  if(out$info != 0) {
    stop(switch(as.character(out$info), `-3` = "Couldn't compute LDL decomposition of P1.", 
                `-2` = "Couldn't compute LDL decomposition of Q."))
  }
  rownames(out$sim) <- switch(sim.what, rep("eps", p), rep("eta", k), c(rep("eps", p), rep("eta", k)), rownames(object$a1), colnames(object$y), colnames(object$y))
  aperm(out$sim, c(2, 1, 3))
}
#####################################################################################################
#####################################################################################################


#####################################################################################################
# simHelper2(): a reduced version of simHelper() from the KFAS package
# model: KFAS model
# modelInfo: computed from getModelInfo(); stores the variables that only need to be computed once
#####################################################################################################
simHelper2 = function (model, modelInfo){
  nsim = 1
  epsplus <- array(0, c(attr(model, "p"), attr(model, "n"), nsim))
  etaplus <- array(0, c(attr(model, "k"), attr(model, "n"), nsim))
  aplus1 <- array(0, dim = c(attr(model, "m"), nsim))
  
  u <- rnorm(modelInfo$dfu, mean = 0, sd = 1)
  if (modelInfo$dfeps > 0) epsplus[modelInfo$x] <- u[1:(modelInfo$dfeps)]
  if (modelInfo$dfeta > 0) etaplus[modelInfo$x2] <- u[(modelInfo$dfeps + 1):(modelInfo$dfeps + modelInfo$dfeta)]
  if (modelInfo$nNonzeroP1 > 0) aplus1[modelInfo$nonzeroP1, ] <- u[(modelInfo$dfeps + modelInfo$dfeta + 1):(modelInfo$dfu)]
  
  list(epsplus = epsplus, etaplus = etaplus, aplus1 = aplus1, nonzeroP1 = as.integer(modelInfo$nonzeroP1))
}
#####################################################################################################
#####################################################################################################


#####################################################################################################
# getModelInfo(): computes many of the parameters required for the simulateSSM() state space sampler,
# specifically those which only need to be computed once
# model: the KFAS model
#####################################################################################################
getModelInfo = function(model){
  nsim = 1
  ymiss <- is.na(model$y)
  storage.mode(ymiss) <- "integer"
  
  x <- array(abs(apply(model$H, 3, diag)) > model$tol, c(attr(model, "p"), attr(model, "n"))) & (!t(ymiss))
  x <- array(x, c(attr(model, "p"), attr(model, "n"), nsim))
  
  dfeps <- sum(x)/nsim
  
  x2 <- array(abs(apply(model$Q, 3, diag)) > model$tol, c(attr(model, "k"), (attr(model, "n") - 1) * attr(model, "tv")[5] + 1))
  x2 <- array(x2, c(attr(model, "k"), attr(model, "n"), nsim))
  
  dfeta <- sum(x2)/nsim
  
  nonzeroP1 <- which(diag(model$P1) > model$tol)
  nNonzeroP1 <- length(nonzeroP1)
  dfu <- dfeps + dfeta + nNonzeroP1
  
  nNonzeroP1inf = as.integer(sum(model$P1inf))
  c2 = numeric(1)
  
  zeroP1inf = which(diag(model$P1inf) > 0)
  
  list(ymiss=ymiss,  x=x, dfeps=dfeps, x2=x2, dfeta = dfeta, nonzeroP1 = nonzeroP1, nNonzeroP1 = nNonzeroP1, dfu = dfu, c2 = c2, nNonzeroP1inf = nNonzeroP1inf, zeroP1inf = zeroP1inf)
}
#####################################################################################################
#####################################################################################################



#####################################################################################################
# computeTimeRemaining() estimates the remaining time in the MCMC based on previous samples
#####################################################################################################
computeTimeRemaining = function(nsi, timer0, nsims){
  
  # Only print occasionally:
  if(nsi%%100 == 0 || nsi==20) {
    # Current time:
    timer = proc.time()[3]
    
    # Simulations per second:
    simsPerSec = nsi/(timer - timer0)
    
    # Seconds remaining, based on extrapolation:
    secRemaining = (nsims - nsi -1)/simsPerSec
    
    # Print the results:
    if(secRemaining > 3600) {
      print(paste(round(secRemaining/3600, 1), "hours remaining"))
    } else {
      if(secRemaining > 60) print(paste(round(secRemaining/60, 2), "minutes remaining"))
    }
  }
  
}	
#####################################################################################################
#####################################################################################################
