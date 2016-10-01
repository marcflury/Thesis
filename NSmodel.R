rm(list = ls(all=TRUE))
library(dplyr)
library(ggplot2)
library(tidyr)
library(forecast)

load("C:/Users/Marc/SkyDrive/R/Thesis/Thesis/Data.RData")

NS <- function(beta, tau,t2m){
  beta[1] +
    beta[2]*((1-exp(-t2m/tau))/(t2m/tau)) +
    beta[3]*((1-exp(-t2m/tau))/(t2m/tau)-exp(-t2m/tau))
}

NSoptim <- function(beta, tau, t2m, price){
  sum((price - NS(beta, tau, t2m))^2)
}

myOptim <- function(IndexDate, df, tau){
  cutdf <- dplyr::filter(df, Date == IndexDate, t2m!=0)
  t(cbind(optim(c(1, 1, 1), NSoptim, tau = tau, price = cutdf$price, t2m = cutdf$t2m)$par))
}


bigOptim <- function(tau, dat){
  datNS <- do.call("rbind", lapply(unique(dat$Date), myOptim, df=dat, tau=tau))
  datNS <- data.frame(unique(dat$Date), datNS) 
  colnames(datNS) <- c("Date", "beta0", "beta1", "beta2")
  datFull <- datNS %>%
    merge(dat, all=TRUE) %>% 
    rowwise() %>%
    mutate(pricehat = NS(c(beta0, beta1, beta2), tau, t2m),
           Error2 = (price - pricehat)^2) 
  return(sum(datFull$Error2, na.rm=TRUE))
}


Brent <- data.frame(Date = rownames(Data), Data[,1:523]) %>%
  gather(t2m, price, -Date) %>%
  dplyr::arrange(Date) %>%
  na.omit() %>%
  mutate(t2m = as.numeric(gsub("[^0-9]", "", as.character(t2m))),
         Date = as.Date(as.character(Date)))

bigopt <- sapply(seq(1,1e5, 1), bigOptim, dat = Brent[1:1260, ])

# tau =  seq(1,1e5, 1)[which.min(bigopt)]

BrentNS <- do.call("rbind", lapply(unique(Brent$Date), myOptim, df=Brent,
                                   tau = tau))

BrentNS <- data.frame(unique(Brent$Date), BrentNS) 
colnames(BrentNS) <- c("Date", "beta0", "beta1", "beta2")

BrentFull <- BrentNS %>%
  merge(Brent, all=TRUE) %>% 
  rowwise() %>%
  mutate(pricehat = NS(c(beta0, beta1, beta2), tau =tau, t2m),
         Diff = price - pricehat)


mean(abs(BrentFull$Diff), na.rm=TRUE)
mean((BrentFull$Diff)^2, na.rm=TRUE)
sqrt(mean((BrentFull$Diff)^2, na.rm=TRUE))


ggplot(BrentFull, aes(x=Diff)) +
  geom_density()

plot(ecdf(BrentFull$Diff))

summary(regBrent0 <- Arima(BrentFull$beta0, c(1,0,0)))
summary(regBrent1 <- Arima(BrentFull$beta1, c(1,0,0)))
summary(regBrent2 <- Arima(BrentFull$beta2, c(1,0,0)))
summary(regBrent0 <- auto.arima(BrentFull$beta0))
summary(regBrent1 <- auto.arima(BrentFull$beta1))
summary(regBrent2 <- auto.arima(BrentFull$beta2))

## Random Walk

BrentNSPreds22 <- mutate(BrentNS,
                       Date = BrentNS$Date[which(unique(BrentNS$Date)==Date) + 22]) %>%
  dplyr::filter(!is.na(Date))%>%
  merge(Brent, all=TRUE)%>% 
  rowwise() %>%
  mutate(pricehat = NS(c(beta0, beta1, beta2), tau =tau, t2m),
         Diff = price - pricehat)

mean(abs(BrentNSPreds22$Diff), na.rm=TRUE)
mean((BrentNSPreds22$Diff)^2, na.rm=TRUE)
sqrt(mean((BrentNSPreds22$Diff)^2, na.rm=TRUE))


BrentNSPreds1 <- mutate(BrentNS,
                         Date = BrentNS$Date[which(unique(BrentNS$Date)==Date) + 1]) %>%
  dplyr::filter(!is.na(Date))%>%
  merge(Brent, all=TRUE)%>% 
  rowwise() %>%
  mutate(pricehat = NS(c(beta0, beta1, beta2), tau =tau, t2m),
         Diff = price - pricehat)

mean(abs(BrentNSPreds1$Diff), na.rm=TRUE)
mean((BrentNSPreds1$Diff)^2, na.rm=TRUE)
sqrt(mean((BrentNSPreds1$Diff)^2, na.rm=TRUE))

#########################################################################
#########################################################################
###########################    Naptha    ################################
#########################################################################
#########################################################################
#########################################################################


Naphtha <- data.frame(Date = rownames(Data), Data[,524:1046]) %>%
  gather(t2m, price, -Date) %>%
  dplyr::arrange(Date) %>%
  na.omit() %>%
  mutate(t2m = as.numeric(gsub("[^0-9]", "", as.character(t2m))),
         Date = as.Date(as.character(Date)))

bigopt <- sapply(seq(1,1e5, 1), bigOptim, dat = Naphtha[1:1260, ])

tau = seq(1,1e5, 1)[which.min(bigopt)]

NaphthaNS <- do.call("rbind", lapply(unique(Naphtha$Date), myOptim, df=Naphtha,
                                   tau = tau))

NaphthaNS <- data.frame(unique(Naphtha$Date), NaphthaNS) 
colnames(NaphthaNS) <- c("Date", "beta0", "beta1", "beta2")

NaphthaFull <- NaphthaNS %>%
  merge(Naphtha, all=TRUE) %>% 
  rowwise() %>%
  mutate(pricehat = NS(c(beta0, beta1, beta2), tau =tau, t2m),
         Diff = price - pricehat)


mean(abs(NaphthaFull$Diff), na.rm=TRUE)
mean((NaphthaFull$Diff)^2, na.rm=TRUE)
sqrt(mean((NaphthaFull$Diff)^2, na.rm=TRUE))


ggplot(NaphthaFull, aes(x=Diff)) +
  geom_density()

plot(ecdf(NaphthaFull$Diff))

summary(regNaphtha0 <- Arima(NaphthaFull$beta0, c(1,0,0)))
summary(regNaphtha1 <- Arima(NaphthaFull$beta1, c(1,0,0)))
summary(regNaphtha2 <- Arima(NaphthaFull$beta2, c(1,0,0)))
summary(regNaphtha0 <- auto.arima(NaphthaFull$beta0))
summary(regNaphtha1 <- auto.arima(NaphthaFull$beta1))
summary(regNaphtha2 <- auto.arima(NaphthaFull$beta2))

## Random Walk

NaphthaNSPreds22 <- mutate(NaphthaNS,
                         Date = NaphthaNS$Date[which(unique(NaphthaNS$Date)==Date) + 22]) %>%
  dplyr::filter(!is.na(Date))%>%
  merge(Naphtha, all=TRUE)%>% 
  rowwise() %>%
  mutate(pricehat = NS(c(beta0, beta1, beta2), tau =tau, t2m),
         Diff = price - pricehat)

mean(abs(NaphthaNSPreds22$Diff), na.rm=TRUE)
mean((NaphthaNSPreds22$Diff)^2, na.rm=TRUE)
sqrt(mean((NaphthaNSPreds22$Diff)^2, na.rm=TRUE))


NaphthaNSPreds1 <- mutate(NaphthaNS,
                        Date = NaphthaNS$Date[which(unique(NaphthaNS$Date)==Date) + 1]) %>%
  dplyr::filter(!is.na(Date))%>%
  merge(Naphtha, all=TRUE)%>% 
  rowwise() %>%
  mutate(pricehat = NS(c(beta0, beta1, beta2), tau =tau, t2m),
         Diff = price - pricehat)

mean(abs(NaphthaNSPreds1$Diff), na.rm=TRUE)
mean((NaphthaNSPreds1$Diff)^2, na.rm=TRUE)
sqrt(mean((NaphthaNSPreds1$Diff)^2, na.rm=TRUE))


#########################################################################
#########################################################################
###########################    Crack     ################################
#########################################################################
#########################################################################
#########################################################################


Crack <- data.frame(Date = rownames(Data), Data[, 524:1046] - Data[, 1:5243]) %>%
  gather(t2m, price, -Date) %>%
  dplyr::arrange(Date) %>%
  na.omit() %>%
  mutate(t2m = as.numeric(gsub("[^0-9]", "", as.character(t2m))),
         Date = as.Date(as.character(Date)))

 bigopt <- sapply(seq(1,1e5, 1), bigOptim, dat = Crack[1:1260, ])
 seq(1,1e5, 1)[which.min(bigopt)]
 tau = tau

CrackNS <- do.call("rbind", lapply(unique(Crack$Date), myOptim, df=Crack,
                                     tau = tau))

CrackNS <- data.frame(unique(Crack$Date), CrackNS) 
colnames(CrackNS) <- c("Date", "beta0", "beta1", "beta2")

CrackFull <- CrackNS %>%
  merge(Crack, all=TRUE) %>% 
  rowwise() %>%
  mutate(pricehat = NS(c(beta0, beta1, beta2), tau =tau, t2m),
         Diff = price - pricehat)


mean(abs(CrackFull$Diff), na.rm=TRUE)
mean((CrackFull$Diff)^2, na.rm=TRUE)
sqrt(mean((CrackFull$Diff)^2, na.rm=TRUE))


ggplot(CrackFull, aes(x=Diff)) +
  geom_density()

plot(ecdf(CrackFull$Diff))

summary(regCrack0 <- Arima(CrackFull$beta0, c(1,0,0)))
summary(regCrack1 <- Arima(CrackFull$beta1, c(1,0,0)))
summary(regCrack2 <- Arima(CrackFull$beta2, c(1,0,0)))
summary(regCrack0 <- auto.arima(CrackFull$beta0))
summary(regCrack1 <- auto.arima(CrackFull$beta1))
summary(regCrack2 <- auto.arima(CrackFull$beta2))

## Random Walk

CrackNSPreds22 <- mutate(CrackNS,
                           Date = CrackNS$Date[which(unique(CrackNS$Date)==Date) + 22]) %>%
  dplyr::filter(!is.na(Date))%>%
  merge(Crack, all=TRUE)%>% 
  rowwise() %>%
  mutate(pricehat = NS(c(beta0, beta1, beta2), tau =tau, t2m),
         Diff = price - pricehat)

mean(abs(CrackNSPreds22$Diff), na.rm=TRUE)
mean((CrackNSPreds22$Diff)^2, na.rm=TRUE)
sqrt(mean((CrackNSPreds22$Diff)^2, na.rm=TRUE))


CrackNSPreds1 <- mutate(CrackNS,
                          Date = CrackNS$Date[which(unique(CrackNS$Date)==Date) + 1]) %>%
  dplyr::filter(!is.na(Date))%>%
  merge(Crack, all=TRUE)%>% 
  rowwise() %>%
  mutate(pricehat = NS(c(beta0, beta1, beta2), tau =tau, t2m),
         Diff = price - pricehat)

mean(abs(CrackNSPreds1$Diff), na.rm=TRUE)
mean((CrackNSPreds1$Diff)^2, na.rm=TRUE)
sqrt(mean((CrackNSPreds1$Diff)^2, na.rm=TRUE))

