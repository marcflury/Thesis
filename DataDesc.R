rm(list=ls(all=TRUE))
library(dplyr)
library(tidyr)
library(lubridate)
library(ggplot2)
library(forecast)
library(fda)
source("t2maturity.r")
source("MF_MFDLM_ARIMA.R")
load("Data.Rdata")
Data["2015-09-04", 870:1046] <- NA # 25% percent change in one day -> outlier
Y <- cbind(maturity2tenor(Data[, 1:522], rowDates = as.Date(rownames(Data)),
                          EndDatesBrent, bizdaysList, startT2m = 0)[, 1:24],
           maturity2tenor(Data[, 523:1046], rowDates = as.Date(rownames(Data)),
                          EndDatesNap, bizdaysList, startT2m = 0)[, 1:24])

Brent <- data.frame(Date = as.Date(rownames(Y)),
                    Y[, 1:24],
                    Commodity = "Brent")


Naphtha <- data.frame(Date = as.Date(rownames(Y)),
                      Y[, 25:48],
                      Commodity = "Naphtha")

Crack  <- data.frame(Date = as.Date(rownames(Y)),
                     Y[, 25:48] - Y[, 1:24],
                     Commodity = "Crack")


df <- rbind(Brent, Naphtha, Crack) %>%
  gather(Tenor, value, -Date, -Commodity) %>%
  mutate(Tenor = factor(gsub("X", "M", Tenor), levels = paste("M", 1:24, sep="")))


dplyr::filter(df, Tenor %in% c("M1", "M12", "M24")) %>%
ggplot(aes(x=Date, y= value, col=Tenor))+
  geom_line()+
  facet_grid(Commodity~., scale="free_y")+
  theme(legend.position = "bottom")

mutate(df, Facet = ifelse(Commodity == "Crack",
                          "Crack",
                          "Brent & Naphtha")) %>%
  dplyr::filter(Tenor == "M1") %>%
  ggplot(aes(x=Date, y = value, colour = Commodity))+
  geom_line()+
  facet_grid(Facet~., scale="free_y")+
  labs(x=NULL, y = "USD/bbl")+
  theme(legend.position = "bottom")+
  scale_x_date(date_breaks = "1 year", date_minor_breaks = "1 month")

na.omit(df) %>%
  dplyr::filter(Commodity != "Crack") %>%
  group_by(Commodity, Tenor) %>%
  summarise(Volatility = sd(diff(log(value)))*sqrt(260)*100) %>%
  ggplot(aes(x=Tenor, y= Volatility, color=Commodity))+
  geom_point()+geom_line()+
  labs(x="Continuous Futures", y = "Volatility in Percent")+
  theme(legend.position = "bottom")
