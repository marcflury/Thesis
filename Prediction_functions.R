# k is the factor of the beta to be predicted
# hT is the last ht calculated
# Capital T stands for the index at from which to predict,
# psi is the coefficient of the AR(1) process for w (\omega)
# gammaSlopes denotes the magnitude with which a shock to beta_{k}^(1) influnces beta_{k}^(2)
# i.e. the last observed value
# The parameters of the stochastic volatility model start with sv
# Ps is the number of ahead forcast steps
# q01 and q10 are the transition probabilities for the HMM
# q01 is the porbability to transition from state 0 to state 1
# The last 4 arguments are predicted values
predictBetas <- function(k, C, wT, hT, psi, gammaSlopes,
                         svMu, svPhi, svSigma, Ps, sT, q10, q01,
                         Beta.p, ht.p, wt.p, st.p){
  for( i in 1:Ps){
    if(i ==1){
      # predict ht as a draw from a normal distribution
      ht.p[1, k] <- rnorm(1, mean = svMu[k] + svPhi[k] *
                            (hT[k] - svMu[k]), sd = svSigma[k])
      # predict w_{T+1} using the predicted log standard deviation from before
      wt.p[1, k] <- psi[k]*wT[k] + exp(ht.p[i, k]/2)*rnorm(1)
      Beta.p[1, k] <- wt.p[1, k]
    } else{
      ht.p[i, k] <- rnorm(1, mean = svMu[k] + svPhi[k] *
                            (ht.p[i-1, k] - svMu[k]), sd = svSigma[k])
      
      wt.p[i,k] <- psi[k]*wt.p[i-1,k] + exp(ht.p[i, k]/2)*rnorm(1)
      
      Beta.p[i,k] <- wt.p[i,k]
    }
    # For c != 1 we use the values from above
    for(c in 2:C) {
      ck.inds <- betaInds[c]-1+k
      if(i == 1){
        ht.p[i, ck.inds] <- rnorm(1, mean = svMu[ck.inds] + svPhi[ck.inds] *
                                    (hT[ck.inds] - svMu[ck.inds]),
                                  sd = svSigma[ck.inds])
        
        wt.p[i, ck.inds] <- psi[ck.inds] * wT[ck.inds] +
          exp(hT[ck.inds]/2) * rnorm(1)
        Beta.p[i, ck.inds] <-   gammaSlopes[ ck.inds] * st.p[i, ck.inds] * wt.p[i,k] +
          wt.p[i, ck.inds]
      } else {
        ht.p[i, ck.inds] <- rnorm(1, mean = svMu[ck.inds] + svPhi[ck.inds] *
                                    (ht.p[i-1, ck.inds] - svMu[ck.inds]), sd = svSigma[ck.inds])
        
        wt.p[i, ck.inds] <- psi[ck.inds] * wt.p[i-1, ck.inds] +
          exp(ht.p[i, ck.inds]/2) * rnorm(1)
        
        Beta.p[i, ck.inds] <- gammaSlopes[ck.inds] * st.p[i, ck.inds] * wt.p[i, k] +
          wt.p[i, ck.inds]
      }
    }
  }
  ret <- list(Beta.p = Beta.p[, c(betaInds[1:C]-1+k)],
              wt.p = wt.p[, c(betaInds[1:C]-1+k)],
              ht.p = ht.p[, c(betaInds[1:C]-1+k)],
              st.p = st.p[, c(betaInds[1:C]-1+k)])
  return(ret)
}

# Predicts the next state given the transition probabilities and the current state
predictState <- function(st, q10, q01){
  as.numeric(runif(length(st))< (1-st)*q01 + st*(1-q10))
}

# Predict several states in advance
predictStatefn <- function(ind.ck, Ps, q10, q01, sT){
  st.p.ck <- rep(0,Ps)
  for(i in 1:Ps){
    if(i == 1){
      st.p.ck[i] <-  predictState(sT[ind.ck], q10[ind.ck], q01[ind.ck])
    } else{
      st.p.ck[i] <-  predictState(st.p.ck[i-1], q10[ind.ck], q01[ind.ck])
    }
  }
  return(st.p.ck)
}

# Predicts the states and betas for one date t
# Index is the row index (date) from which to predict forward
predictAll <- function(Index, Beta, ht, K, C, psi, gammaSlopes,
                       svMu, svPhi, svSigma, Ps, S, q10, q01){
  
  wT <- diff(Beta)[Index-1, ]
  
  for(c in 2:C){
    ind.c <- (1:K) + (c-1)*K
    wT[ind.c] <- wT[ind.c ] -
      gammaSlopes[ind.c] * S[Index, ind.c] * wT[1:K]
  }
  
  Beta.p <- array(0, c(Ps, ncol(Beta)))
  ht.p = array(0, c(Ps, ncol(ht)))
  wt.p = array(0, c(Ps, ncol(Beta)))
  st.p = array(0, c(Ps, ncol(S)))
  
  st.p[, (K+1):(C*K)] <- do.call("cbind", lapply((K+1):(C*K), predictStatefn,
                                                 Ps = Ps, q10 = q10, q01 = q01,
                                                 sT = S[Index, ]))
  
  for(k in 1:K){
    predList  <- predictBetas(k, C, wT, hT= ht[Index, ],
                              psi, gammaSlopes, svMu,
                              svPhi, svSigma, Ps,
                              sT= S[Index, ], q10, q01,
                              Beta.p =  Beta.p,
                              ht.p = ht.p,
                              wt.p = wt.p,
                              st.p = st.p
    )
    
    Beta.p[, c(betaInds[1:C]-1+k)] <- predList$Beta.p
    wt.p[, c(betaInds[1:C]-1+k)] <- predList$wt.p
    ht.p[, c(betaInds[1:C]-1+k)] <- predList$ht.p
    st.p[, c(betaInds[1:C]-1+k)] <- predList$st.p
  }
  
  return(colSums(Beta.p))
}


pred1 <- function(x,Beta=Beta, ht=ht, K=K, C=C, psi=psi,
                  gammaSlopes=gammaSlopes, svMu=svMu, svPhi=svPhi, svSigma=svSigma,
                  Ps=Ps, S=S, q10=q10, q01=q01){
  do.call("rbind", lapply(2:nrow(Beta), predictAll, Beta=Beta, ht=ht, K=K, C=C, psi=psi,
                          gammaSlopes=gammaSlopes, svMu=svMu, svPhi=svPhi, svSigma=svSigma,
                          Ps=Ps, S=S, q10=q10, q01=q01 )
  )
}

npred <- function(nsims, Beta=Beta, ht=ht, K=K, C=C, psi=psi,
                  gammaSlopes=gammaSlopes, svMu=svMu, svPhi=svPhi, svSigma=svSigma,
                  Ps=Ps, S=S, q10=q10, q01=q01){
  
  
  bigList <- lapply(1:nsims,pred1,Beta=Beta, ht=ht, K=K, C=C, psi=psi,
                    gammaSlopes=gammaSlopes, svMu=svMu, svPhi=svPhi,
                    svSigma=svSigma,
                    Ps=Ps, S=S, q10=q10, q01=q01)
  
  ret <- array(0, dim = c(dim(bigList[[1]]),nsims))
  
  for(i in 1:nsims){
    ret[,,i] <- bigList[[i]]
  }
  return(ret)
}

# Predict the log standard deviation using the SV model
predictht <- function(ck.inds, hT, svMu, svPhi, svSigma){
  ht <- rnorm(1, mean = svMu[ck.inds] + svPhi[ck.inds] *
                (hT[ck.inds] - svMu[ck.inds]), sd = svSigma[ck.inds])
  return(ht)
}



#####################################################################################################
# updateS1 updates the HMM states
# dB.ck: delta Beta for outcome c (> 1), factor k
# dB.1k: delta Beta for outcome 1, factor k
# psi.ck: AR(1) coefficient for outcome c, factor k
# St: T-dimensional vector of states for outcome c, factor k
# q01.ck: transition probability from 0 to 1 in HMM model for outcome c, factor k
# q10.ck: transition probability from 1 to 0 in HMM model for outcome c, factor k
# sigma.t2: T-dimensional vector of error variances from the Beta-level for outcome c, factor k
# See Albert and Chib, 1993, for more details
#####################################################################################################
updateS1 <- function(ind.ck, K, dBeta, psi, St, q01, q10, sigma.t2,
                     gammaSlopes){
  n = nrow(dBeta)
  # Residuals from AR(1) model:
  resBc = dBeta[2:n, ind.ck] - psi[ind.ck] * dBeta[1:(n-1), ind.ck]
  resBc = c(0,resBc) # for indexing
  
  # First, consider state = 1
  St[n, ind.ck] = 1
  # Residual contribution from Beta_1k in Beta_ck equation
  resB1 = dBeta[n:(n-1), (ind.ck-1) %% K + 1]*St[n:(n-1), ind.ck] * gammaSlopes[ind.ck]
  resB1 = resB1[1] - psi[ind.ck] * resB1[2]
  ldist1 =  log(computeTransitProbs(St[n-1, ind.ck], 1,
                                    q01[ind.ck], q10[ind.ck]))-
    .5*(resBc[n]-resB1)^2/sigma.t2[n, ind.ck]
  
  # Next, consider state = 0 
  St[n, ind.ck] = 0 
  resB1 =  dBeta[n:(n-1), (ind.ck-1) %% K + 1] * St[n:(n-1), ind.ck] * gammaSlopes[ind.ck]
  resB1 = resB1[1] - psi[ind.ck]*resB1[2]
  ldist0 = log(computeTransitProbs(St[n-1, ind.ck], 0,
                                   q01 [ind.ck],
                                   q10 [ind.ck]))-
    .5*(resBc[n]-resB1)^2/ sigma.t2[n, ind.ck]
  
  # Compare
  if((ldist0 - ldist1) < log(1/runif(1)-1)){
    St[n, ind.ck] = 1
  }
  return(St[n, ind.ck])
}


updateS <- function(Beta, psi, St, q01, q10, sigma.t2, gammaSlopes, C, K){
  sapply((K+1):(C*K), updateS1, K=K, dBeta = diff(Beta), psi = psi, St = St,
         q01 = q01, q10 = q10, sigma.t2 =sigma.t2, gammaSlopes = gammaSlopes)
}

#####################################################################################################
# updateht() samples the stochastic volatility model for fixed c, across k (Yield curve application)
# resBeta.c: (T-1) x K matrix of residuals from the Beta-level for outcome c
# ht.c: the T x K matrix of log-volatilities for outcome c
# svMu.c: the K-dimensional vector of the unconditional means of the log-volatilities for outcome c
# svPhi.c: the K-dimensional vector of the AR(1) coefficients of the log-volatilities for outcome c
# svSigma.c: the K-dimensional vector of the error variances of the log-volatilities for outcome c
# svOffset: an offset we can add to the squared residuals (before taking logs) for computational stability
# then adjust the implied variances after sampling (Note: other SV parameters are not rescaled, but do not appear elsewhere in the model)
#####################################################################################################
updateht = function(resBeta.c,  ht.c, svMu.c, svPhi.c, svSigma.c, svOffset = 0){
  # Define locally
  K = ncol(resBeta.c)
  
  # To update the Beta-level error variance for KFAS
  Wt.c = array(0, c(K, K, nrow(ht.c)))
  
  for(k in 1:K){
    resBeta.ck = resBeta.c[,k]
    
    # Sample the SV model
    
    # Update the "start parameters" to most recent MCMC sample for outcome c, factor k
    startpara = list(mu = svMu.c[k],
                     phi = svPhi.c[k],
                     sigma = svSigma.c[k])
    
    
    # Include an offset for computational stability (then rescale after sampling):
    resBeta.c[,k] = sqrt(1+svOffset/resBeta.c[,k]^2)*resBeta.c[,k]
    
    sv = svsample2(resBeta.c[,k], draws = 1, startpara=startpara ,
                   startlatent = ht.c[-1, k])
    
    # Store the results
    svMu.c[k] = para(sv)[1]
    svPhi.c[k] = para(sv)[2]
    svSigma.c[k] = para(sv)[3]
    ht.c[,k] = c(sv$latent0, sv$latent)
    
    Wt.c[k, k, ] = exp(ht.c[,k])
    Wt.c[k, k, -1] = Wt.c[k, k, -1]/(1+svOffset/resBeta.c[,k]^2)
    ht.c[-1,k] <- log(exp(ht.c[-1,k])/(1+svOffset/resBeta.c[,k]^2))
  }
  list(ht.c = ht.c, svMu.c = svMu.c, svPhi.c = svPhi.c, svSigma.c = svSigma.c, Wt.c = Wt.c)
}
#####################################################################################################
# inSampleLength: length of the training sample
# Beta: averaged Betas from the training stage
# ht: latent states from the SV part of the training step
# K: number of coefficients
# C: number of Commodities
# psi: coefficient from ARIMA(1,1,0) model
# gammaSlopes: interaction term from HMM
# svMu, svPhi, svSigma: parameters from the SV step
# Ps: number of forecast steps
# S: states matrix for HMM
# q10, q01: transition probabilites for the states in S
# Y.full: complete set of data in Y
#####################################################################################################
updateStep <- function(inSampleLength, Beta, ht, K, C, psi, gammaSlopes,
                       svMu, svPhi, svSigma, S, q10, q01, Y.full, Et, d){
  # Append matricies for ht, states and beta for updated values
  #
  ht_ <- rbind(ht, array(0, c(nrow(Y.full)-nrow(ht), ncol(ht))))
  S_ <- rbind(S, array(0, c(nrow(Y.full)-nrow(S), ncol(S))))
  Beta_ <- rbind(Beta, array(0, c(nrow(Y.full)-nrow(Beta), ncol(Beta))))
  
  for(i in (inSampleLength+1):nrow(Y.full)){
    ht <- rbind(ht,
                t(sapply(1:(C*K), predictht, hT=ht[i-1,],
                         svMu=svMu, svPhi=svPhi, svSigma=svSigma))
    )
    
    S <- rbind(S, predictState(S_[i-1,], q10, q01))
    
    GtWt <- computeGtWtHMM(S, gammaSlopes, psi, exp(ht))
    Wt <- GtWt$Wt
    Gt <- GtWt$Gt
    
    dWt <- ifelse(length(dim(Wt))==3, dim(Wt)[3], 1)
    Q <- array(0, c(2*C*K, 2*C*K, dWt))
    Q[1:(C*K), 1:(C*K), ] <- Wt
    Q[(C*K+1):(2*C*K), (C*K+1):(2*C*K),] <- diag(C*K)
    
    Model = SSModel(Y.full[1:i, ]~ -1 +
                      SSMcustom(Z = array(0, c(ncol(Y.full), nrow(Gt))),
                                T = Gt, R = diag(c(rep(1,C*K),rep(0,C*K))),
                                Q = Q, P1 = diag(10^4, nrow(Gt))))
    
    Theta <- mfdlmBeta(Y.full[rownames(Y.full)[1:i], ], Et, Gt, Wt, Model, tau, d, splineInfo)
    Beta <- Theta[,1:(C*K)]
    dBeta <- diff(Beta)
    ind.T <- nrow(dBeta)
    
    # Cycle through the outcomes:
    for(c in 1:C) {
      bc.inds = betaInds[c]:(betaInds[c+1]-1)		# Subsets: \beta_{1,t}^{(c)}, ..., \beta_{K,t}^{(c)}
      yc.inds = outcomeInds[c]:(outcomeInds[c+1]-1) 	# Subsets: Y_t^{(c)}(\tau_1), ..., Y_t^{(c)}(\tau_m)
      
      # Conditional means of Y_t^{(c)}(\tau) (T x m)
      mu.c = tcrossprod(Beta[,bc.inds], splineInfo$Phi[match(tau[yc.inds], allTaus),]%*%d)
      
      # Observation error variance
      Et[c] = sampleEt(Y.full[1:i, yc.inds], mu.c)
      
      # Sample the stochastic volatility parameters:
      if(c > 1){
        wt <- dBeta[, bc.inds] - S[-1,bc.inds] * dBeta[,1:K] *
          matrix(rep(gammaSlopes[bc.inds], ind.T), nrow=ind.T, byrow=TRUE)
        
        resBeta.c = wt[-1,] -
          matrix(rep(psi[bc.inds], ind.T -1), nrow=ind.T -1, byrow=TRUE) *
          wt[-ind.T ,]
      } else{
        resBeta.c = dBeta[-1, bc.inds] -
          matrix(rep(psi[bc.inds], ind.T -1), nrow=ind.T -1, byrow=TRUE) *
          dBeta[-ind.T , bc.inds]
      }
      
      samples = updateht(resBeta.c,  ht[-1, bc.inds], svMu[bc.inds], svPhi[bc.inds],
                         svSigma[bc.inds], svOffset=10^-6)
      
      ht[-1, bc.inds] <- samples$ht.c
    }
    Beta_[i, ] <- Beta[i, ]
    ht_[i, ] <- ht[i, ]
    S_[i, ] <- updateS(Beta, psi, S, q01, q10, sigma.t2 = exp(ht), gammaSlopes, C, K)
  }
  return(list(Beta_ = Beta_, ht_ = ht_, S_ = S_))
}

#################################################
predOoS1 <- function(nsim, inSampleLength, Beta, ht, K, C, psi, gammaSlopes, svMu,
                     svPhi, svSigma, Ps, S, q10, q01, Y.full, d){
  
  updates <- updateStep(inSampleLength = inSampleLength, Beta = Beta, ht = ht, K = K,
                        C = C, psi = psi, gammaSlopes = gammaSlopes, svMu = svMu, svPhi = svPhi,
                        svSigma = svSigma, S = S, q10 = q10, q01 = q01,
                        Y.full = Y.full, d=d)
  
  
  pred1(Beta=updates$Beta_, ht=updates$ht_, K=K, C=C, psi=psi,
        gammaSlopes=gammaSlopes, svMu=svMu, svPhi=svPhi,
        svSigma=svSigma, Ps=Ps, S=updates$S_, q10=q10, q01=q01)
}


npredooS <- function(nsims, inSampleLength, Beta, ht, K, C, psi, gammaSlopes, svMu,
                     svPhi, svSigma, Ps, S, q10, q01, Y.full, d){
  bigL <- lapply(1:nsims, predOoS1, inSampleLength = inSampleLength, Beta = Beta, ht = ht, K = K,
                 C = C, psi = psi, gammaSlopes = gammaSlopes, svMu = svMu, svPhi = svPhi,
                 svSigma = svSigma, Ps = Ps, S = S, q10 = q10, q01 = q01,
                 Y.full = Y.full, d=d)
  
  ret <- array(0, c(dim(bigL[[1]]), nsims))
  
  for(i in 1:nsims){
    ret[,,i] <- bigL[[i]]
  }
  return(ret)
}
### Run full preditcions in and out of sample
fullPrediction <- function(i){
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
  
  set.seed(i)
  
  postI <- postIs[i]
  d <- postDAll[postI,,]
  updates <- updateStep(inSampleLength = 1260,
                        Beta = postBetaAll[postI,,],
                        ht = postHtAll[postI,,], K = K, C = C,
                        psi = postPsi[postI,],
                        gammaSlopes = postGammaSlopesAll[postI,],
                        svMu = postsvMu[postI,], svPhi = postsvPhi[postI,],
                        svSigma = postsvSigma[postI,],
                        S = postSAll[postI,,], q10 = postq10[postI,],
                        q01 = postq01[postI,],
                        Y.full = Y, Et = postEtAll[postI, ], d=d)
  for(Ps in c(1,22)){
    print(Ps)
    numPreds <- 10
    preds <- npred(numPreds, Beta=updates$Beta_, ht=updates$ht_, K=K, C=C,
                   psi = postPsi[postI,],
                   gammaSlopes = postGammaSlopesAll[postI,],
                   svMu = postsvMu[postI,], svPhi = postsvPhi[postI,],
                   svSigma = postsvSigma[postI,],
                   Ps = Ps,
                   S = updates$S_, q10 = postq10[postI,],
                   q01 = postq01[postI,])
    Yp <- createYp(dBeta.P=preds, Beta = updates$Beta_,
                   Phit= Phit, d=d, C=C, K=K, Et = postEtAll[postI, ]
    )[1:(NROW(Y)-1-Ps), ,]
    rowDatesYp <- rowDatesY[(1+1+Ps):NROW(Y)]
    
    name <- paste("results44T", Ps, i, sep="_")
    tmpYp <- array(0,dim=c(NROW(Yp),48, numPreds))
    for(numPred in 1: numPreds){
      tmpYp[, , numPred] <- cbind(maturity2tenor(Yp[,1:523, numPred],
                                                 rowDatesYp,
                                                 EndDatesBrent, bizdaysList)[, 1:24],
                                  maturity2tenor(Yp[,524:1046, numPred],
                                                 rowDatesYp,
                                                 EndDatesNap, bizdaysList)[, 1:24])
    }
    res <- list(Yp =tmpYp ,
                preds = preds)
    rm(list = c("Yp", "tmpYp", "preds"))
    gc()
    assign(name, res)
    save(list = name, file = paste(name, ".Robj", sep=""))
  }
  upname <- paste("updates44T", i, sep="_")
  assign(upname, updates)
  save(list = upname, file = paste(paste("updates44T", i, sep="_"), ".Robj", sep=""))
}

# allPredictions is used for prediction if the updates are already available
# Written to use on multiple cores
# Run full preditcions in and out of sample
allPredictions <- function(i, resPath, numPreds= 10, PsVec = C(1, 5, 22)){
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
  
  set.seed(i)
  
  postI <- postIs[i]
  d <- postDAll[postI,,]
  updates <- load(paste(resPath, paste("updates44T", i, sep="_"), ".Robj", sep=""))
  assign("updates", get(paste("updates44T", i, sep="_")))
  
  for(Ps in PsVec){
    print(Ps)
    preds <- npred(numPreds, Beta=updates$Beta_, ht=updates$ht_, K=K, C=C,
                   psi = postPsi[postI,],
                   gammaSlopes = postGammaSlopesAll[postI,],
                   svMu = postsvMu[postI,], svPhi = postsvPhi[postI,],
                   svSigma = postsvSigma[postI,],
                   Ps = Ps,
                   S = updates$S_, q10 = postq10[postI,],
                   q01 = postq01[postI,])
    Yp <- createYp(dBeta.P=preds, Beta = updates$Beta_,
                   Phit= Phit, d=d, C=C, K=K, Et = postEtAll[postI, ]
    )[1:(NROW(Y)-1-Ps), ,]
    rowDatesYp <- rowDatesY[(1+1+Ps):NROW(Y)]
    
    name <- paste("results44T", Ps, i, sep="_")
    tmpYp <- array(0,dim=c(NROW(Yp),48, numPreds))
    for(numPred in 1: numPreds){
      tmpYp[, , numPred] <- cbind(maturity2tenor(Yp[,1:523, numPred],
                                                 rowDatesYp,
                                                 EndDatesBrent, bizdaysList)[, 1:24],
                                  maturity2tenor(Yp[,524:1046, numPred],
                                                 rowDatesYp,
                                                 EndDatesNap, bizdaysList)[, 1:24])
    }
    res <- list(Yp =tmpYp ,
                preds = preds)
    rm(list = c("Yp", "tmpYp", "preds"))
    gc()
    assign(name, res)
    save(list = name, file = paste(name, ".Robj", sep=""))
  }
  upname <- paste("updates44T", i, sep="_")
  assign(upname, updates)
  save(list = upname, file = paste(paste("updates44T", i, sep="_"), ".Robj", sep=""))
}


####################
# Function to create predicted Y, i.e. Y.p
createYp1 <- function(psim, dBeta.P, Beta, Phit, d, C, K, Et){
  Taumax <- nrow(Phit)
  Y.p1 <- array(0, c(nrow(dBeta.P), Taumax*C))
  for(c in 1:C){
    inds.c <- 1:K + K*(c-1)
    Y.p1[,1:Taumax + Taumax*(c-1)] <-
      (dBeta.P[, inds.c, psim] + Beta[2:(nrow(Beta)), inds.c]) %*%
      t(Phit%*%d) +
      matrix(rnorm(Taumax*nrow(Y.p1), mean = 0, sd = sqrt(Et[c])),
             ncol = Taumax)
    gc()
  }
  
  #Et  <- sampleEt()
  
  # Y.p1 <- Y.p1
  return(Y.p1)
}

createYp <- function(dBeta.P, Beta, Phit, d, C, K, Et){
  
  if(length(dim(dBeta.P))==3){
    Y.p <- array(0,   c(nrow(dBeta.P), nrow(Phit)*2, dim(dBeta.P)[3] ))
  } else{
    Y.p <- array(0,   c(nrow(dBeta.P), nrow(Phit)*2 ))
  }
  Y.pL <- lapply(1:dim(dBeta.P)[3], createYp1, dBeta.P=dBeta.P, Beta=Beta,
                 Phit=Phit, d=d, C=C, K=K, Et=Et)
  if(length(dim(dBeta.P))==3){
    for(i in 1:dim(dBeta.P)[3]){
      Y.p[,,i] <- Y.pL[[i]]
    }
    
  }else {
    Y.p <- Y.pL
  }
  return(Y.p)
}
#############################################################################
##################### Function to turn Yp into Tx24 Matrix ##################
#############################################################################
# Yp = Predictions of the Y (Matrix)
# Indices = Indecies to Keep from original data
# Ps = Prediction steps (Positive interger)
# Diffs = How many times the data is differencieted
Yp2Tenor <- function(Yp, Dates, Endates, Ps, Diffs){
  data.frame(Date = Dates[(1+Diffs+Ps):length(Dates)],
             Yp[1:(nrow(Yp-Ps)),]) %>%
    gather(t2m, value, -Date) %>%
    mutate(Contract = ifelse(Date - value %in% Endates$Dates,
                             Endates$FC[which(Date - value %in% Endates$Dates)],
                             NA)) %>%
    return()
  
}




#######################################################
############# Evaluation functions
precentile1 <- function(ind.tau,obs, preds){
  return(rowMeans(preds[,ind.tau, ] > obs[, ind.tau]))
}


precentileFun <- function(preds, obs){
  return(do.call("cbind",
                 lapply(1:ncol(preds), precentile1, preds=preds, obs=obs)
  )
  )
}