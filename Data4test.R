# rm(list=ls(all=TRUE))
Tend <- 200
Ht <- 2
Qt <- 1.2
u <- rnorm(Tend,0,Ht)
eta <- rnorm(Tend,0,Qt)
phi0 <- 3
phi1 <- 0.9
alpha0 <- 0.5
alpha <- rep(0,Tend)
y <- rep(0,Tend)
for(t in 1:Tend){
  if(t == 1){
    alpha[t] <- phi1*alpha0 + phi0
  } else {
    alpha[t] <- phi1*alpha[t-1] + phi0 + eta[t-1]
  }
  y[t] <- alpha[t] #+ u[t]
}

# MCMC size controls
D <-	1e5	# number of final draws
B <-.1*D		# burn-in period
S <- B+D		# burn-in plus draws that are kept


## hyperparmeters for the priors
# kill the influence of priors by settign them to zero
V0_inv	<- 0
beta0		<- rep(0,2)
a0			<- 0
b0			<- 0	

prior<-list(inv_V0	= 0,
            beta0 = rep(0,2),
            a0 = 0,
            b0 = 0)

# storage for Gibbs draws
theta_draws		<- matrix(rep(0,S*3),ncol=3)     

# initialise parameters for Gibbs draws
theta_draws[1,]<- t(c(1, 0.5, 2))
# alpha_draws		<- zeros(T,S)

# constant to be passed to nig_draws
C <- rep(1,Tend-1)

## MAIN GIBBS Sampling Loop
#t0	<- tic 
#t00 <- now # timer paramters 
for (i in 2:S){
  # FFBS of log volatilieties given st_i and theta
  alpha   <- ffbs(y,theta_draws[i-1,])
  # alpha_draws(:,i)  = alpha
  
  # Gibbs sampling of state equation parameters
  theta_draws[i,]  <- nig_draws(alpha[2:Tend],cbind(C, alpha[1:(Tend-1)]),prior)

}

pars			<- theta_draws[(B+1):nrow(theta_draws),]
c_draws   <- pars[,1]
phi_draws <- pars[,2]
sig_draws <- sqrt(pars[,3])

plot(c_draws)
plot(phi_draws)

plot(sig_draws)

hist(c_draws)
hist(phi_draws)
hist(sig_draws)

mean(c_draws)
mean(phi_draws)
mean(sig_draws)
