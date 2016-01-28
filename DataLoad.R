rm(list=ls(all=TRUE))
# Read excel file
library(XLConnect)
wb <-  loadWorkbook("CrackData.xlsx")

# Reads part of the Excelfile
# creates a single variable for every continuous contract
myload <- function(num, Name){
  rawdata <- readWorksheet(wb, Name,
                           startRow=2,
                           startCol=1+(num-1)*4,
                           endCol=2+(num-1)*4
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
allNaphtha <- do.call("rbind", lapply(allDates, createNap))

rawData <- merge(allBrent, allNaphtha, by="Date")
Data  <- as.matrix(rawData[,-1])
colnames(Data) <- gsub("[^1234567890]","",colnames(Data))
rownames(Data) <- as.character(rawData$Date)


save(Data, file="MF_MFDLM/Data.RData")
