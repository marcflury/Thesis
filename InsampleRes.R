# in sample Results
rm(list = ls(all = TRUE))
library(dplyr)
library(tidyr)
library(lubridate)
library(ggplot2)
library(forecast)
library(fda)
library(ggthemes)
source("t2maturity.r")
source("MF_MFDLM_ARIMA.R")
load("Data.Rdata")
source("C:/Users/Marc/SkyDrive/R/Thesis/Thesis/Prediction_functions.R")
source("C:/Users/Marc/SkyDrive/R/Thesis/Thesis/t2maturity.r")

Data["2015-09-04", 870:1046] <- NA # 25% percent change in one day -> outlier
Y <- Data
Ytenor <- Y
colnames(Ytenor) = as.numeric(colnames(Ytenor))  + 1

Ytenor <- cbind(maturity2tenor(Ytenor[,1:522], rowDates = as.Date(rownames(Data)),
                               EndDatesBrent, bizdaysList)[,1:24],
                maturity2tenor(Ytenor[,523:1046], rowDates = as.Date(rownames(Data)),
                               EndDatesNap, bizdaysList)[,1:24])



T = nrow(Y[1:1260, ])							# Total number of time points
tau0 = as.numeric(colnames(Y))			# Frequencies within each time bin, concatenated across c = 1,...,C
tau = (tau0 - min(tau0))/(diff(range(tau0)))	# Shift to [0,1] (for numerical stability)
allTaus0 = sort(unique(tau0))				# All unique frequencies (observation points)
allTaus = sort(unique(tau)) 				# All unique frequencies (observation points), restricted to [0,1]
m = length(allTaus)					# Total number of unique frequencies

splineInfo <- getSplineInfo(tau)

####################################
Phit = splineInfo$basisPhi(seq(splineInfo$a, splineInfo$b, length.out=522))

load("EstResultsFull3HMM.RObj")
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


Yhat <- cbind(maturity2tenor(colMeans(postBetaAll[2001:7000, , ])[,1:3] %*%
                               t(Phit%*%colMeans(postDAll[2001:7000, ,]))[,1:522],
                             rowDates = as.Date(rownames(Data[1:1260,])),
                             EndDatesBrent, bizdaysList)[,1:24],
              maturity2tenor((colMeans(postBetaAll[2001:7000, , ])[,4:6] %*%
                                t(Phit%*%colMeans(postDAll[2001:7000, ,])))[,1:522],
                             rowDates = as.Date(rownames(Data[1:1260,])),
                             EndDatesNap, bizdaysList)[,1:24])

diffshat <- Ytenor[1:1260, ] - Yhat

Ydiffs <- data.frame(Date = rownames(diffshat), diffshat) %>%
  gather(Tenor, value, -Date) %>% 
  mutate(Tenor = gsub("X", "", Tenor)) %>% 
  mutate(Commodity = ifelse(grepl("\\.1$", Tenor), "Naphtha", "Brent")) %>%
  mutate(Tenor = factor(gsub("\\.1$", "", Tenor), levels = as.character(1:24)))

ggplot(Ydiffs, aes(x = Tenor, y = value, colour = Commodity))+
  geom_boxplot()+
  theme(legend.position = "bottom")+
  labs(x = "Continuous Contract", y = "USD/bbl")+
  scale_color_fivethirtyeight() 


meanD <- t(apply(postDAll[2001:7000, , ], 2, colMeans))

ftau <- Phit %*% meanD

data.frame(Tau = 1:522, ftau) %>% 
  gather(variable, value, -Tau) %>%
  mutate(variable = gsub("X", "Beta ", variable)) %>%
  ggplot(aes(x=Tau, y=value, colour= variable))+
  geom_line()+
  scale_y_continuous(breaks = seq(-2,10,1))+
  labs(x = "Days to maturity", y = "USD/bbl")+
  theme(legend.position = "bottom", legend.title = element_blank())+
  scale_color_fivethirtyeight() 





meanBeta = t(apply(postBetaAll[2001:7000, , ], 2, colMeans))


brentBeta <- data.frame(Date = as.Date(rownames(Y[1:1260, ])),
                        meanBeta[,1:3], 
                        Commodity = "Brent")%>%
  gather(betaNum, value, -Date, -Commodity)


naphthaBeta <- data.frame(Date = as.Date(rownames(Y[1:1260, ])),
                          meanBeta[,4:6], 
                          Commodity = "Naphtha") %>%
  gather(betaNum, value, -Date, -Commodity)

betaDf <- rbind(brentBeta, naphthaBeta) %>%
  mutate(Facet = gsub("X", "Beta ", betaNum))

ggplot(betaDf, aes(x =Date, y= value, colour =  Commodity))+
  geom_line()+
  labs(x = NULL, y = "USD/bbl")+
  facet_grid(Facet~., scale = "free_y")+
  theme(legend.position = "bottom", legend.title = element_blank())+
  scale_color_fivethirtyeight() 


meanS = t(apply(postSAll[2001:7000, , ], 2, colMeans))

naphthaS <- data.frame(Date = as.Date(rownames(Y[1:1260, ])),
                       meanS[,4:6], 
                       Commodity = "Naphtha") %>%
  gather(SNum, value, -Date, -Commodity)%>%
  mutate(Facet = gsub("X", "S ", SNum))



ggplot(naphthaS, aes(x =Date, y= value, col = Facet))+
  geom_line()+
  labs(x = NULL, y = "USD/bbl")+
  theme(legend.position = "bottom", legend.title = element_blank())+
  scale_color_fivethirtyeight() 


data.frame(postGammaSlopesAll[2001:7000, 4:6])%>%
  gather(GNum, value)%>%
  mutate(Facet = gsub("X", "Gamma ", GNum)) %>%
  ggplot(aes(x = value, fill = Facet))+
  geom_density(alpha = 0.5)+
  facet_wrap(~Facet, scale = "free_x")+
  labs(x = NULL)+
  theme(legend.position = "none", legend.title = element_blank())+
  scale_fill_fivethirtyeight() 


meanSigma = t(apply(exp(postHtAll[2001:7000, , ]/2), 2, colMeans))

brentSigma <- data.frame(Date = as.Date(rownames(Y[1:1260, ])),
                        meanSigma[,1:3], 
                        Commodity = "Brent")%>%
  gather(sigmaNum, value, -Date, -Commodity)


naphthaSigma <- data.frame(Date = as.Date(rownames(Y[1:1260, ])),
                          meanSigma[,4:6], 
                          Commodity = "Naphtha")[-(1:22),] %>%
  gather(sigmaNum, value, -Date, -Commodity)



sigmaDf <- rbind(brentSigma, naphthaSigma) %>%
  mutate(Facet = gsub("X", "sigma ", sigmaNum))

ggplot(sigmaDf, aes(x =Date, y= value, colour =  Commodity))+
  geom_line()+
  labs(x = NULL, y = "USD/bbl")+
  facet_grid(Facet~., scale = "free_y")+
  theme(legend.position = "bottom", legend.title = element_blank())+
  scale_color_fivethirtyeight() 



