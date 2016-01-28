library(rgl)

# Plot the Joint loading curves:
taugrid = seq(allTaus0[1], allTaus0[m], length.out=24)
Phit = splineInfo$basisPhi(seq(splineInfo$a, splineInfo$b, length.out=length(taugrid)))
# Brent Curve
dfB <- t(Phit%*%d%*%t(Beta[,1:4]))
dfN <- t(Phit%*%d%*%t(Beta[,5:8]))
#Diffs to real values
bdiffs <- Y[,1:24]-dfB
ndiffs <- Y[,25:48]-dfN

z <-bdiffs
x <- 1:nrow(z)  # 10 meter spacing (S to N)
y <- 1:ncol(z)   # 10 meter spacing (E to W)

zlim <- range(z)
zlen <- zlim[2] - zlim[1]+1

colorlut <- terrain.colors(zlen) # height color lookup table

col <- colorlut[ z - zlim[1] + 1 ] # assign colors to heights for each point


surface3d(x, y, z, color=col)
axis3d('y--', at = seq(0,1300,50),labels=seq(0,200,10),color="black")
axis3d('z--', at = seq(-13.5,13.5,0.5),
       labels=seq(-13.5,13.5,0.5),color="black")
axis3d('x--', at = y, labels=seq(1,24,1),color="black")
 