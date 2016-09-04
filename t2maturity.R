library(timeDate)
library(dplyr)
library(tidyr)
library(lubridate)
library(xts)

FCnum2FC <- function(x){
  return(paste(FCcodes[round((x - floor(x))*12) +1], floor(x), sep=""))
}

source("C:/Users/Marc/SkyDrive/R/Thesis/Thesis/helpingMatrices.R")

# Calculates the time to maturity for a given commodity
# Matrix of type dates x Tenors
# bizdaysList list with all the business days
t2maturity <- function(x, EndDatesdf, bizdaysList){
  df <- gather(x, Tenor, Value, -Date) %>%
    dplyr::mutate(Tenor = as.numeric(gsub("[^0123456789]", "", Tenor))) %>%
    merge(EndDatesdf[,c("FCnum","EndDate")],by.x="Date",
          by.y = "EndDate", all=TRUE) %>%
    transform(FCnum = na.locf(FCnum, fromLast=TRUE)) %>%
    dplyr::filter(!is.na(Tenor)) %>%
    transform(FC = FCnum2FC(FCnum+ (Tenor-1)/12)) %>%
    dplyr::select(Date, Value, FC) %>%
    merge(EndDatesdf[,c("FC","EndDate")],by.x=c("FC"),
          by.y = c("FC"), all.x=TRUE) %>%
    transform(Maturity = bizdaysList[as.character(EndDate), ] - 
                bizdaysList[as.character(Date), ]) %>%
    dplyr::select(Date, Value, FC, EndDate, Maturity) %>%
    arrange(Date) 
  
  return(tidyr::spread(df[,c("Date","Value", "Maturity")], Maturity, Value))
}

# Calculates the for a given time to maturity and end dates df
# Matrix of type dates x Tenors
# bizdaysList list with all the business days
maturity2tenor <- function(x, EndDatesdf, bizdaysList){
  
  df <- 
    data.frame( Date = as.Date(rownames(x)), x) %>%
    gather(t2m, Value, -Date) %>% 
    mutate(t2m = as.numeric(gsub("X","", as.character(t2m))))  %>%
    dplyr::mutate(expiryDate = as.Date(rownames(
      bizdaysList)[bizdaysList[format(Date), ] + t2m])) %>% 
    inner_join(EndDatesdf[,c("FCnum","EndDate")], by = c("expiryDate" = "EndDate")) %>%
    group_by(Date) %>%
    arrange(FCnum) %>%
    mutate(Tenor = paste("M", 1:n(), sep="")) 
  
  return(tidyr::spread(df[,c("Date","Value", "Tenor")], Tenor, Value)) 
} 

