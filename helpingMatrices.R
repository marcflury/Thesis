# Create matricies for calculating maturity dates
tS <- seq(timeDate("2010-01-01"),
          timeDate("2020-05-31"),
          by="day")

# lookup matrix
bizdays <- tS[isBizday(tS, holidays =holidayNYSE())]
bizdaysList <- matrix(1:length(bizdays),ncol=1)
rownames(bizdaysList) <- format(bizdays)

############ FC and EndDate matrix
FCcodes <- c("F", "G", "H", "J", "K", "M", "N",	"Q", "U",	"V", "X",	"Z")

bizdaysdf <- data.frame(as.Date(tS[isBizday(tS, holidayLONDON())]))
colnames(bizdaysdf) <-  "Date"

EndDatesNap <- transform(bizdaysdf,
                         Year = year(Date),
                         Month = month(Date),
                         Day = day(Date)) %>%
  group_by(Year, Month) %>%
  dplyr::summarise(EndDate=last(Date)) %>%
  dplyr::mutate(FCnum=Year + (Month - 1)/12,
                FC= FCnum2FC(FCnum)) %>%
  ungroup()


####### Load End Dates ############################
EndDatesBrent <- read.csv("ICE_Brent_Expiry_Calendar.csv", sep=",", 
                          stringsAsFactors = FALSE) %>%
  mutate(EndDate = as.Date(EndDate)) %>%
  dplyr::filter(year(EndDate)>=2010)





##### Delete if everything works
if(FALSE){
  
  
  ## Brent for contracts up to and including the February 2016 contract
  EndDatesBre1 <- transform(bizdaysdf,
                            Index = 1:nrow(bizdaysdf),
                            Year = year(Date),
                            Month = month(Date),
                            Day = day(Date)) %>%
    transform(temp1 = ifelse(Day == 15,Index-2,ifelse(Day < 15, Index -3, 0))) %>%
    dplyr::filter(Year + Month/12 <= 2016 + 1/12) %>%
    group_by(Year, Month) %>%
    dplyr::summarise( EndDate = bizdaysdf$Date[max(temp1)]) %>%
    dplyr::mutate(FCnum=Year + (Month + 1 - 1)/12,
                  FC1= FCnum2FC(FCnum)) %>%
    dplyr::mutate(FC=paste(rep(FCcodes,2)[Month+1],
                           Year + (unique(Month+1)-1) %/% 12,
                           sep="")) %>%
    ungroup()
  
  ## Brent for contracts starting with March 2016
  ## -unique(Month) %/% 12 subtracts 1 if the month is 12 New Year's rule
  EndDatesBre2 <- transform(bizdaysdf,
                            Index = 1:nrow(bizdaysdf),
                            Year = year(Date),
                            Month = month(Date),
                            Day = day(Date)) %>%
    dplyr::filter(Year >= 2016) %>%
    group_by(Year, Month) %>%
    dplyr::summarise( EndDate = bizdaysdf$Date[dplyr::last(Index) - 1
                                               - unique(Month) %/% 12]
    ) %>%
    dplyr::mutate(FCnum=Year + (Month + 2 - 1)/12,
                  FC1= FCnum2FC(FCnum)) %>%
    dplyr::mutate(FC=paste(rep(FCcodes,2)[Month+2],
                           Year + (unique(Month+2)-1) %/% 12,
                           sep="")) %>%
    ungroup()
  
  EndDatesBre <- rbind(EndDatesBre1[,c("FC","FCnum", "EndDate")],
                       EndDatesBre2[,c("FC", "FCnum", "EndDate")])
  
  
  
  
  
  library(dplyr)
  Sys.setlocale("LC_TIME", "C")
  EndDatesBrent <- read.csv("ICE_Brent_Expiry_Calendar.csv", sep=";", 
                            stringsAsFactors = FALSE) %>%
    mutate(EndDate = as.Date(EXPIRY.DATE, format="%B %d, %Y"))
  
  write.csv(EndDatesBrent[,c(1,3)], file = "ICE_Brent_Expiry_Calendar_hist.csv")
  
  GerMonths <- c( "Dez", "Nov", "Okt", "Sep", "Aug", "Jul", "Jun",
                  "Mai", "Apr", "Mrz", "Feb", "Jan")
  
  EngMonths <- c( "Dec", "Nov", "Oct", "Sep", "Aug", "Jul", "Jun",
                  "May", "Apr", "Mar", "Feb", "Jan")
  
  jnkdf <- data.frame(EngMonths, Num = 12:1)
  rownames(jnkdf) <- GerMonths
  jnkdf2 <- data.frame(EngMonths, Num = 12:1)
  rownames(jnkdf2) <- EngMonths
  jnk <- read.csv("ICE_Brent_Expiry_Calendar_hist.csv", sep=";") %>%
    mutate(FCM = jnkdf[str_sub(as.character(EXPIRING.CONTRACT), 1,3),
                       "EngMonths"],
           FCM = ifelse(is.na(FCM),
                        str_sub(as.character(EXPIRING.CONTRACT), 1,3), 
                        as.character(FCM)),
           FCY = gsub("[^0-9]", "", as.character(EXPIRING.CONTRACT)),
           Num = jnkdf2[FCM,"Num"],
           Contract = paste(FCM, FCY, sep =" "),
           EndDate = as.Date(EndDate, format="%d.%m.%Y"),
           FC = paste(FCcodes[Num], FCY, sep=""),
           FCnum = as.numeric(FCY) + (Num - 1)/12+2000)
  
  write.csv(jnk[,c("FC", "FCnum", "EndDate")], file = "ICE_Brent_Expiry_Calendar.csv",
            row.names=FALSE)
}
