
# Required packages:
library(fda)
library(KFAS)
library(MCMCpack)
library(stochvol)
library(truncnorm)
library(abind)



# nsims: Total number of simulations
# burnin: Burn-in
# K.hmm.sv: Number of factors for the common trend model
# useHMM: Hidden Markov model (HMM), or common trend (CT) model? (CT in the paper)


runModelfun <- function(nsims, burnin, K, K.hmm.sv, useHMM){
  set.seed(1990)

    # File containing functions
  source("MF_MFDLM_ARIMA.R")

  #########################################################################################################################

  # Load the data and store key variables:
  load("Data.RData")

  Data["2015-09-04",525:1047] <- NA # 25% percent change in one day -> outlier
  Y <- Data


  C = 2							# Number of outcomes
  cnames = c("Brent", "Naphtha")		# Names of the outcomes (central banks)
  dates = as.Date(rownames(Y))				# Dates of the observations

  #########################################################################################################################

  T = nrow(Y)							# Total number of time points
  tau0 = as.numeric(colnames(Y))			# Frequencies within each time bin, concatenated across c = 1,...,C
  tau = (tau0 - min(tau0))/(diff(range(tau0)))	# Shift to [0,1] (for numerical stability)
  allTaus0 = sort(unique(tau0))				# All unique frequencies (observation points)
  allTaus = sort(unique(tau)) 				# All unique frequencies (observation points), restricted to [0,1]
  m = length(allTaus)					# Total number of unique frequencies


  # Some useful indices:
  outcomeInds = c(1, which(tau[-1] < tau[-length(tau)])+1, length(tau)+1)
  # For outcome-specific observations, use Y[, yc.inds]
  # where yc.inds = outcomeInds[c]:(outcomeInds[c+1]-1)
  betaInds = seq(1, C*(K+1), by=K)								# For outcome-specific Beta, use  Beta[, bc.inds]
  # where bc.inds = betaInds[c]:(betaInds[c+1]-1)

  ####################################

  # Interpolate the data with using cubic splines
  Y <- do.call("cbind",lapply(1:C,splineinterpol))
  T = nrow(Y)
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
    dBeta  <- rbind(0, Beta[2:T, ]- Beta[1:(T-1), ])
  # Initialize the AR(1) and HMM parameters:
  initsHMMpar = initHMM(dBeta, K.hmm.sv, useHMM)
  psi = initsHMMpar$psi
  q01 = initsHMMpar$q01
  q10 = initsHMMpar$q10
  S = initsHMMpar$S
  gammaSlopes = initsHMMpar$gammaSlopes


  # Initialize the SV parameters (just using independent ARIMA(1,0,0) models for Beta_ck):
  initsSVpar = initSV(apply(Beta, 2, function(x){arima(x, c(1,0,0), include.mean=FALSE)$resid[-1]}), C)
  ht = initsSVpar$ht
  svMu = initsSVpar$svMu
  svPhi = initsSVpar$svPhi
  svSigma = initsSVpar$svSigma
  Wt = initsSVpar$Wt

  #########################################################################################################################

  # Set up the state space model
  Gt = array(0, c(2*C*K, 2*C*K, T))

  Gt[,,1:T] <- rbind(cbind(diag(1, C*K, C*K), - diag(C*K)),
                     cbind(diag(C*K), array(0, dim = c(C*K, C*K))))
  Wt = initsSVpar$Wt

  # Compute Gt and Wt matrices for HMM; common trend model is a submodel (Note: this is not efficient)
  GtWt <- computeGtWtHMM(s, gammaSlopes, psi, exp(ht))

  Wt <- GtWt$Wt
  Gt <- GtWt$Gt

  dWt <- ifelse(length(dim(Wt))==3, dim(Wt)[3], 1)
  Q <- array(0, c(2*C*K, 2*C*K, dWt))
  Q[1:(C*K), 1:(C*K), ] <- Wt
  Q[(C*K+1):(2*C*K), (C*K+1):(2*C*K),] <- diag(C*K)

  Model = SSModel(Y~-1+SSMcustom(Z = array(0, c(ncol(Y), nrow(Gt))), T = Gt,
                                 R = diag(c(rep(1,C*K),rep(0,C*K))),
                                 Q = Q, P1 = diag(10^4, nrow(Gt))))

  #########################################################################################################################

  # Parameters to save:
  postBetaAll = array(0, c(nsims, dim(Beta)))
  postDAll = array(0, c(nsims, dim(d)))
  postLambdaAll = array(0, c(nsims, length(lambda)))
  postEtAll = array(0, c(nsims, C))
  postHtAll = array(0, c(nsims, dim(ht)))
  postGammaSlopesAll = array(0, c(nsims, C*K))
  devAll = numeric(nsims)  # Deviance
  postCountAll = array(0, c(nsims, T-1, C*K))
  postCountAggAll = array(0, c(nsims, T-1, C-1)) # For outlier plot
  postSAll <-  array(0, c(nsims, dim(S)))

  #########################################################################################################################

  # Now, run the MCMC
  timer0 = proc.time()[3]			# for timing the sampler
  for(nsi in 1:nsims){
    # nsi is the simulation index
    # Sample Beta, d, and lambda:
    samples = mfdlm(Y, tau, Beta, Et, Gt, Wt, Model, d, splineInfo, lambda)
    Beta = samples$Beta
    d = samples$d
    lambda=samples$lambda
    dBeta  <- rbind(0, Beta[2:T, ]- Beta[1:(T-1), ])
    Theta <- samples$Theta

    # Cycle through the outcomes:
    for(c in 1:C) {
      bc.inds = betaInds[c]:(betaInds[c+1]-1)		# Subsets: \beta_{1,t}^{(c)}, ..., \beta_{K,t}^{(c)}
      yc.inds = outcomeInds[c]:(outcomeInds[c+1]-1) 	# Subsets: Y_t^{(c)}(\tau_1), ..., Y_t^{(c)}(\tau_m)

      # Conditional means of Y_t^{(c)}(\tau) (T x m)
      mu.c = tcrossprod(Beta[,bc.inds], splineInfo$Phi[match(tau[yc.inds], allTaus),]%*%d)

      # Observation error variance
      Et[c] = sampleEt(Y[, yc.inds], mu.c)

      # Sample the AR(1) and slope parameters (and other HMM parameters, if desired)
      shmm = sampleHMMpar(dBeta, K.hmm.sv, S, gammaSlopes, psi, ht, c, useHMM, q01, q10)
      gammaSlopes <- shmm$gammaSlopes
      psi         <- shmm$psi

      if(useHMM){
        S           <- shmm$S
        q01         <- shmm$q01
        q10         <- shmm$q10
      }

      # Sample the stochastic volatility parameters:
      if(c > 1){
        wt <- dBeta[, bc.inds] - S[,bc.inds] * dBeta[,1:K] *
          matrix(rep(gammaSlopes[bc.inds], T), nrow=T, byrow=TRUE)

        resBeta.c = wt[-1,] -
          matrix(rep(psi[bc.inds], T-1), nrow=T-1, byrow=TRUE) *
          wt[-T,]
      } else{
        resBeta.c = dBeta[-1, bc.inds] -
          matrix(rep(psi[bc.inds], T-1), nrow=T-1, byrow=TRUE) *
          dBeta[-T, bc.inds]
      }

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
      devAll[nsi] = devAll[nsi] + 1/Et[c] * sum((Y[, yc.inds]-mu.c)^2,na.rm=TRUE) +
        sum(!is.na(Y[, yc.inds]))* log(2 * pi * Et[c])
    }

    # Compute Gt and Wt matrices for HMM; common trend model is a submodel (Note: this is not efficient)
    GtWt <- computeGtWtHMM(Gt, Wt, S, gammaSlopes, psi, exp(ht))

    Wt <- GtWt$Wt
    Gt <- GtWt$Gt


    # Store the samples, respecting the burn-in and thinning
    postBetaAll[nsi,,] <- Beta
    postDAll[nsi,,] <- d
    postLambdaAll[nsi,] <- lambda
    postEtAll[nsi,] <- Et
    postHtAll[nsi,,] <- ht
    postGammaSlopesAll[nsi,] <- gammaSlopes
    postSAll[nsi,, ] <- S

    # Check the time remaining:
    computeTimeRemaining(nsi, timer0, nsims)
  }
  print(paste('Total time: ', round((proc.time()[3] - timer0)/60), 'minutes'))
resNames <- paste("Results", K,K.hmm.sv, useHMM, sep="_")
  assign(resNames,
  list(postBetaAll = postBetaAll,
    postDAll = postDAll,
    postLambdaAll = postLambdaAll,
    postEtAll = postEtAll,
    postHtAll = postHtAll,
    postGammaSlopesAll = postGammaSlopesAll,
    postSAll = postSAll, 
    d = d,
    splineInfo = splineInfo))
save(resNames,
file= paste("D:/R/Thesis_Data/",
            resNames,
 ".robj", sep="" ))  
}
runModelfun(2,1,4,4,TRUE)