library(dplyr)
library(ggplot2)
library(tidyr)
library(RColorBrewer)
library(rgl)

betadf <- data.frame(as.Date(rownames(Y)), Beta)
colnames(betadf) <- c("Date",
                      paste("BrentB",1:4, sep=""),
                      paste("NapB",1:4, sep=""))

gather(betadf, variable, value, -Date) %>%
  mutate(Facet = as.numeric(variable) %% 4) %>%
  ggplot(aes(Date, value, col=variable)) + geom_line() +
  facet_grid(Facet~., scale="free_y")

spbasis <- Phit%*%d

brent <- Phit%*%d %*% t(Beta[,1:4])
nap <- Phit%*%d %*% t(Beta[,5:8])

diffs <- t(rbind(brent, nap)) -Y

df <- data.frame(Date=as.Date(rownames(Y)), diffs) %>%
  gather(variable, value, -Date)

ggplot(df, aes(x=Date, y=variable, colour=value))+geom_point()

co <- apply(brent, 1, cumsum)
na <- apply(nap, 1, cumsum)

comms <- rbind(data.frame(Date=as.Date(rownames(Y)), co, com="Brent"),
               data.frame(Date=as.Date(rownames(Y)), na, com="Naphtha")) %>%
  gather(variable, value, -Date, -com) %>%
  mutate(variable = as.numeric(gsub("X", "", as.character(variable))))

dplyr::filter(comms, com=="Brent") %>%
  ggplot(aes(x=Date, y=value, col=factor(variable)))+geom_line()

dplyr::filter(comms, com=="Naphtha") %>%
  ggplot(aes(x=Date, y=value, col=factor(variable)))+geom_line()

commsBeta <- rbind(data.frame(Date=as.Date(rownames(Y)), 
                              apply(Beta[, 1:4], 2, cumsum)
                              , com="Brent"),
                   data.frame(Date=as.Date(rownames(Y)),
                              apply(Beta[, 5:8], 2, cumsum),
                              com="Naphtha")) %>%
  gather(variable, value, -Date, -com) %>%
  mutate(variable = as.numeric(gsub("X", "", as.character(variable))))


ggplot(commsBeta,aes(x=Date, y=value, col=com))+
  geom_line()+ facet_grid(variable~., scale="free")

commsEt <- rbind(data.frame(Date=as.Date(rownames(Y)), 
                            postEt[,1]
                              , com="Brent"),
                   data.frame(Date=as.Date(rownames(Y)),
                              postEt[, 2],
                              com="Naphtha")) %>%
  gather(variable, value, -Date, -com) %>%
  mutate(variable = as.numeric(gsub("X", "", as.character(variable))))


ggplot(commsht,aes(x=Date, y=value, col=com))+
  geom_line()+ facet_grid(variable~., scale="free")



zlim <- range(co)
zlen <- zlim[2] - zlim[1]+1

colorlut <- terrain.colors(zlen) # height color lookup table

col <- colorlut[ co - zlim[1] + 1 ] # assign colors to heights for each point


surface3d(1:1026, 1:24, co, color=col)
surface3d(1:1026, 1:24, na+20, color=col)

axis3d('x--', at = seq(0,1300,50),labels=seq(0,1300,50),color="black")
axis3d('z--', at = seq(-100,100,10),
       labels=seq(-100,100,10),color="black")
axis3d('y--', at = seq(1,24,4), labels=seq(1,24,4),color="black")



#################### Spread
z <- co-na
zlim <- range(z)
zlen <- zlim[2] - zlim[1]+1

colorlut <- terrain.colors(zlen) # height color lookup table

col <- colorlut[ z - zlim[1] + 1 ] # assign colors to heights for each point


surface3d(1:1026, 1:24, z, color=col)
surface3d(1:1026, 1:24, z, color=col)

axis3d('x--', at = seq(0,1300,50),labels=seq(0,1300,50),color="black")
axis3d(c('z++','z','z+'), at = seq(-10,10,10),
       labels=seq(-10,10,10),color="black")
axis3d('y++', at = seq(1,24,4), labels=seq(1,24,4),color="black")



realdata <- apply(Y, 2, cumsum)

#################### Spread
z <- co-na
zlim <- range(z)
zlen <- zlim[2] - zlim[1]+1

colorlut <- terrain.colors(zlen) # height color lookup table

col <- colorlut[ z - zlim[1] + 1 ] # assign colors to heights for each point

surface3d(1:1026, 1:24, z, color=col)

axis3d('x--', at = seq(0,1300,50),labels=seq(0,1300,50),color="black")
axis3d(c('z++','z','z+'), at = seq(-2,2,.1),
       labels=seq(-2,2,.1),color="black")
axis3d('y++', at = seq(1,24,4), labels=seq(1,24,4),color="black")
