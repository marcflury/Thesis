# Create matricies for calculating maturity dates


tS <- seq(timeDate("2010-01-01"),
          timeDate("2020-05-31"),
          by="day")

# lookup matrix
bizdays <- tS[isBizday(tS, holidayLONDON())]
bizdaysList <- matrix(1:length(bizdays),ncol=1)
rownames(bizdaysList) <- as.character(bizdays)

############ FC and EndDate matrix
FCcodes <- c("F", "G", "H", "J", "K", "M", "N",	"Q", "U",	"V", "X",	"Z")

bizdaysdf <- data.frame(as.Date(tS[isBizday(tS, holidayNYSE())])) 
colnames(bizdaysdf) <-  "Date"

EndDatesNap <- transform(bizdaysdf,
                         Year = year(Date),
                         Month = month(Date),
                         Day = day(Date)) %>%
  group_by(Year, Month) %>%
  dplyr::summarise(EndDate=last(Date)) %>%
  dplyr::mutate(FCnum=Year + (Month - 1)/12,
                FC1= FCnum2FC(FCnum)) %>%
  dplyr::mutate(FC=paste(FCcodes[Month], Year, sep="")) %>%
  ungroup()


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

