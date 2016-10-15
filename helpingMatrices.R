# Create matricies for calculating maturity dates
tS <- seq(timeDate("2010-01-01"),
          timeDate("2020-05-31"),
          by="day")

# lookup matrix
bizdays <- tS[isBizday(tS, holidays =holidayNYSE())]
bizdaysList <- matrix(1:length(bizdays),ncol=1)
rownames(bizdaysList) <- format(bizdays)

# Futures codes letters indicating the month of the year as used by exchanges
FCcodes <- c("F", "G", "H", "J", "K", "M", "N",	"Q", "U",	"V", "X",	"Z")

# Data.frame of all trading days 
bizdaysdf <- data.frame(as.Date(tS[isBizday(tS, holidayLONDON())]))
colnames(bizdaysdf) <-  "Date"


# Creates a data.frame with the expiry dates of the naphtha contract
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
