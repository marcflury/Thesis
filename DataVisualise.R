library(dplyr)
library(ggplot2)
library(tidyr)
spbasis <- Phit%*%d
brentb <- matrix(c(rep(1, 4),
                   c(1.2,1,1,1),
                   rep(1, 4),
                   c(1,2,1,1),
                   rep(1, 4),
                 c(1,1,2,1),
                 c(1,1,1,1),
                 c(1,1,1,1.2)),
                   ncol=4, byrow = TRUE)
napb <- matrix(rep(1, 8*4), ncol=4)+  t(t(brentb-1) * gammaSlopes[5:8]  )


brent <- brentb %*% t(spbasis)
nap <- napb %*% t(spbasis)


B <- gather(data.frame(Date=1:8,brent), Tenor, value,-Date) %>%
  mutate(Com="Brent", Facet=(1+Date) %/% 2, Tenor=as.numeric(Tenor))

N <- gather(data.frame(Date=1:8,nap), Tenor, value,-Date) %>%
  mutate(Com="Nap", Facet=(1+Date) %/% 2, Tenor=as.numeric(Tenor))

rbind(B,N) %>%
  ggplot(aes(x=Tenor, y=value, col=as.factor(Date), linetype=Com))+
  geom_line()+facet_grid(Facet~.)
