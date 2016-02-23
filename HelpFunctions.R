library(timeDate)
library(dplyr)
library(tidyr)
library(lubridate)
library(xts)

FCnum2FC <- function(x){
  return(paste(FCcodes[round((x - floor(x))*12) +1], floor(x), sep=""))
}

source("helpingMatrices.R")

# Calculates the time to maturity for a given commodity
# Matrix of type dates x Tenors
# com character string with brent or naphtha
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
