rm(list=ls(all=TRUE))
# Read excel file
library(XLConnect)
library(ggplot2)
library(dplyr)
library(tidyr)
library(lubridate)
wb <-  loadWorkbook("ThesisData.xlsx")

# Reads part of the Excelfile
# creates a single variable for every continuous contract
myload <- function(num, Name){
  rawdata <- readWorksheet(wb, Name,
                           startRow=2,
                           startCol=1+(num-1)*3,
                           endCol=2+(num-1)*3
  )
  
  fullName <- paste(Name, num, sep="")
  colnames(rawdata) <- c("Date", fullName)
  assign(fullName, rawdata, envir = .GlobalEnv)
  fullName
}

lapply(1:24,myload,Name="Brent") # Run on Brent sheet
lapply(1:24,myload,Name="Crack") # Run on Crack sheet

# Merge all Brent files
allBrent <- merge(Brent1, Brent2, all=TRUE)
for(i in 3:24){
  allBrent <- merge(allBrent,get(paste("Brent", i, sep="")) , all=TRUE)
}

# Merge all Crack files
allCrack <- merge(Crack1, Crack2, all=TRUE)
for(i in 3:24){
  allCrack <- merge(allCrack, get(paste("Crack", i, sep="")), all=TRUE)
}

# find NA and plot them
data.frame(Date=allCrack$Date,
           Count=rowSums(apply(allCrack,c(1,2),is.na))
) %>%
  ggplot(aes(Date, Count))+geom_line()


# Create matrix indicating which columns have NAs
find_NA <- function(x){
  data.frame(Date=x$Date,
             apply(x,c(1,2),is.na)) %>%
    gather(Column, isNA, -Date) %>%
    dplyr::filter(isNA==TRUE) %>%
    ggplot(aes(Date, Column, col=isNA))+geom_point()+
    scale_colour_manual(values=c("orchid"))
}
# Plots depicting at which date and tenor the NAs occur
find_NA(allCrack)

# find_NA(allBrent) breaks as only th last day has NA

# Function that adds Brent and Crack for specific datestamp
# makes sure that the price for the same dates are added
createNap <- function(datestamp){
  napraw <- allBrent[which(allBrent$Date==datestamp), paste("Brent",1:24,sep="")] +
    allCrack[which(allCrack$Date==datestamp), paste("Crack",1:24,sep="")]
  ret <- data.frame(datestamp, napraw)
  colnames(ret) <- c("Date", paste("Naphtha", 1:24, sep=""))
  return(ret)
}

# Find Dates that are in Brent and Crack 
allDates <- allBrent$Date[which(allBrent$Date %in% allCrack$Date) ]

# Which dates are missing?
NADates <- allBrent$Date[which(!(allBrent$Date %in% na.omit(allCrack)$Date)) ]

# Create the Naphtha data.frame 
allNaphtha <- do.call("rbind", lapply(allDates, createNap)) %>%
  t2maturity( EndDatesNap, bizdaysList)


# NA for merge
find_NA(merge(allBrent, allCrack, by="Date", all=TRUE))

# Cut sample so that year(Date) >= 2012
rawData <- merge(allBrent, allNaphtha, by="Date") %>%
  dplyr::filter(year(Date) >= 2012)

Data  <- as.matrix(rawData[,-1])

colnames(Data) <- gsub("[^1234567890]","",colnames(Data))
rownames(Data) <- as.character(rawData$Date)
find_NA(rawData)

save(Data, file="Data.RData")
