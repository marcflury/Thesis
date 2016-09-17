rm(list = ls(all=TRUE))
library(dplyr)
library(ggplot2)
library(tidyr)
library(pipeR)

load("C:/Users/Marc/SkyDrive/R/Thesis/Thesis/Data.RData")

NS <- function(betaTau,t2m){
  betaTau[1] +
    betaTau[2]*((1-exp(-t2m/betaTau[4]))/(t2m/betaTau[4])) +
    betaTau[3]*((1-exp(-t2m/betaTau[4]))/(t2m/betaTau[4])-exp(-t2m/betaTau[4]))
}

NSoptim <- function(betaTau, t2m, price){
  sum((price - NS(betaTau,t2m))^2)
}

myOptim <- function(IndexDate, df){
  cutdf <- dplyr::filter(df, Date == IndexDate, t2m!=0)
  t(optim(c(1, 1, 1, 20), NSoptim, price = cutdf$price, t2m = cutdf$t2m)$par)
}

Brent <- data.frame(Date = rownames(Data), Data[,1:523]) %>%
  gather(t2m, price, -Date) %>%
  dplyr::arrange(Date) %>%
  na.omit() %>%
  mutate(t2m = as.numeric(gsub("[^0-9]", "", as.character(t2m))),
         Date = as.Date(as.character(Date)))

cutdf <- dplyr::filter(Brent, Date == Brent[1, "Date"], t2m!=0)
myopts <- optim(c(1, 1, 1, 20), NSoptim, price = cutdf$price, t2m = cutdf$t2m)

NS(myopts$par, cutdf$t2m)


BrentNS <- do.call("rbind", lapply(unique(Brent$Date), myOptim, df=Brent))
BrentNS <- data.frame(unique(Brent$Date), BrentNS) 
colnames(BrentNS) <- c("Date", "beta0", "beta1", "beta2", "tau")

BrentFull <- BrentNS %>%
  merge(Brent, all=TRUE) %>% 
  rowwise() %>%
  mutate(pricehat = NS(c(beta0, beta1, beta2, tau), t2m),
         Diff = price - pricehat)



