library(timeDate)
library(dplyr)
library(tidyr)
library(lubridate)
library(xts)

FCnum2FC <- function(x){
  return(paste(FCcodes[(x - floor(x))*12 +1], floor(x)-2000, sep=""))
}

source("C:/Users/Marc/SkyDrive/R/Thesis/Thesis/helpingMatrices.R")

# Calculates the time to maturity for a given commodity
# Matrix of type dates x Tenors
# bizdaysList list with all the business days
t2maturity <- function(x, EndDatesdf, bizdaysL){
  EndDatesdf <- dplyr::arrange(EndDatesdf, EndDate)
  df <- gather(x, Tenor, Value, -Date) %>%
    dplyr::mutate(Tenor = as.numeric(gsub("[^0123456789]", "", Tenor))) %>%
    rowwise() %>%
    mutate(EndDate = base::subset(EndDatesdf,
                                  EndDatesdf$EndDate>=Date)[Tenor, "EndDate"])
  df <- mutate(df, 
               Maturity = bizdaysL[format(EndDate), ] -
                 bizdaysL[format(Date), ])
  return(tidyr::spread(df[,c("Date","Value", "Maturity")], Maturity, Value))
  #return(df$EndDate[which(!(format(df$EndDate) %in% rownames(bizdaysL))) ]) 
}

# Calculates the for a given time to maturity and end dates df
# Matrix of type dates x Tenors
# bizdaysList list with all the business days
maturity2tenor <- function(x, rowDates, EndDatesdf, bizdaysList){
  df <- 
    data.frame( Date = as.Date(rowDates), x) %>%
    gather(t2m, Value, -Date) %>% 
    mutate(t2m = as.numeric(gsub("X","", as.character(t2m)))-1)  %>%
    dplyr::mutate(expiryDate = as.Date(rownames(
      bizdaysList)[bizdaysList[format(Date), ] + t2m])) %>% 
    dplyr::filter(expiryDate %in% EndDatesdf$EndDate) %>%
    group_by(Date) %>%
    arrange(expiryDate) %>%
    mutate(Tenor = 1:n()) %>%
    ungroup() %>%
    dplyr::select(Date, Value, Tenor) %>%
    tidyr::spread(Tenor, Value)

  
  mat <- as.matrix(df[,-1])
  rownames(mat) <- format(df$Date)
  return(mat) 
} 

