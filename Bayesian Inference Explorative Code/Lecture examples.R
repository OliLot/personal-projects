###############################
###### Table of Contents ######
###############################

# 1. LECTURE 10 MONTE CARLO INTEGRATION
# 2. LECTURE 13 Gibbs Sampling
# 3. LECTURE 14 Metropolis-Hastings
# 4. LECTURE 15 Mixture Normal model


#--------------------------- START OF EXAMPLES ---------------------------------

################################################
###### LECTURE 10 MONTE CARLO INTEGRATION ######
################################################
# target is x^2 integrated between 0 and 3

#1) simulate 4000 samples from U(0,3)
draws=runif(4000,0,3)

#2) compute MC approximation
(3/4000)*sum(draws^2)

#3) estimate of MC error
sd(draws^2)/sqrt(4000)



#-------------------------------------------------------------------------------



###############################################
########## LECTURE 13 Gibbs Sampling ##########
###############################################
#Example Gibbs sampler, for X~N(mu,sigma^2) when mu and sigma are unknown
# normal(phi,tau^2) prior on mu and IG on sigma^2(alpha,beta)

#Uncomment this line if you want to get the same answer each time you run it
#set.seed(1293820)

#Priors
phi<-0
tau2<-1
alpha<-0.1
beta<-0.01

#Simulate some data
n<-10
x<-rnorm(n,0,1)
xbar<-mean(x)

#Number of samples to take in the chain
T<-10000
#Create space to store the Markov chain in 
mu<-sigma2<-numeric(T)

#Set starting values
mu[1]<-4
sigma2[1]<-0.1

#Run the Gibbs sampler 
for(t in 1:(T-1)){
  mu[t+1]<-rnorm(1,(tau2*n*xbar + sigma2[t]*phi)/(tau2*n + sigma2[t]),
                 sqrt(sigma2[t]*tau2/(tau2*n+sigma2[t])))
  sigma2[t+1]<-1/rgamma(1,shape=n/2 + alpha,
                        rate=1/2*sum((x-mu[t+1])^2) + beta)
}

#Produce trace plots of the outputs
par(mfrow=c(1,2))
plot(1:T,mu,xlab="Iteration",ylab="mu",type="l")
plot(1:T,sigma2,xlab="Iteration",ylab="sigma2",type="l")
par(mfrow=c(1,1))

#Produce density plots of the outputs
#(Note - you'd want to remove some initial samples as burn-in
# before using the samples for inference)
par(mfrow=c(1,2))
plot(density(mu))
plot(density(sigma2))
par(mfrow=c(1,1))

#Produce bivariate plot of samples
plot(mu,sigma2,type="p")

posteriormean_mu=mean(mu)
left95interval=quantile(mu,c(0.025,0.975))[1]
right95interval=quantile(mu,c(0.025,0.975))[2]
cat("posterior mean for mu=",posteriormean_mu, sep=c(""))
cat(" ", sep=c("\n"))
cat("95% credible interval=(",left95interval,",",right95interval,")", sep=c(""))

posteriormean_sigma2=mean(sigma2)
left95interval=quantile(sigma2,c(0.025,0.975))[1]
right95interval=quantile(sigma2,c(0.025,0.975))[2]
cat("posterior mean for sigma^2=",posteriormean_sigma2, sep=c(""))
cat(" ", sep=c("\n"))
cat("95% credible interval=(",left95interval,",",right95interval,")", sep=c(""))

# Of course, if you were to write your own code not as a very basic example, 
# but to properly perform some analysis,
# you would have to allow for a burn-in, and to also obtain additional inferences and 
# diagnostic plots 



#-------------------------------------------------------------------------------



################################################
######## LECTURE 14 Metropolis-Hastings ########
################################################
# Example of Metropolis sampling from section 2.3 of the lecture notes.
# Sampling from the standard Normal distribution. 

T<-500
#Variability of the normal distribution proposal
std<-sqrt(100) #0.1,1,100 

#store samples
theta<-numeric(T)
#store number of acceptances
n.accept<-0

#starting value

theta[1]<-0 #Set start value at 0 when exploring proposals
for (i in 1:(T-1)){
  phi<-rnorm(1,theta[i],std)  
  alpha<-min(c(1,exp(-0.5*(phi^2-theta[i]^2))))
  if(runif(1)<alpha) {
    theta[i+1]<-phi
    n.accept<-n.accept+1
  } else {
    theta[i+1]<-theta[i]
  }
}
p.accept<-n.accept/T
p.accept
#
#
#Calculate effective sample size
library(coda)
chain<-mcmc(theta)
ess<-effectiveSize(chain)
ess

#Plot using home-grown code
#
windows()
par(mfrow=c(2,2))
plot(1:T,theta,type="l",main=paste("p.accept=",format(p.accept,digits=2),", ESS=",format(ess,digits=1),sep=""))
qqnorm(theta)
hist(theta)
acf(theta)
summary(theta)
sd(theta)

# For an analysis (rather than an illustration of the MH sampler properties)
#  remember to discard the burn-in sample. 



#-------------------------------------------------------------------------------



##################################################
######## LECTURE 15 Mixture Normal model #########
##################################################

library(nimble)
library(igraph)
library(coda)
library(R6)
#
# HEIGHTS DATA (lecture 15)
#
# Generate the data #####
set.seed(15)
HeightData=c(rnorm(550,185,8),rnorm(700,170,8))# heights in cm.
# #######################
#
# Specify the statistical model
MixNormCode <- nimbleCode({
  
  # Specify the likelihood:
  for (i in 1:N){
    x[i] ~ dnorm(mu[t[i]],s)
    t[i] <- z[i]+1
    z[i] ~ dbin(pi,1)# auxiliary variable
  }
  # Prior specification:
  mu[1] ~ dnorm(phi,lambda)
  mu[2] ~ dnorm(phi,lambda)
  pi ~ dbeta(1,1)
  
})
#
# Values for some constants in the model
MixNormConsts <- list(N = 700+550, s=1/8^2, phi=177, lambda=1/15^2 ) 
#
# The data values
MixNormData <- list(x=HeightData)

# one set of initial values before building the model                 
MixNormInits <- list(pi=0.5,mu=c(177,177),z=rbinom(1250,1,0.5)) 

# to build the model
MixNorm <- nimbleModel(code = MixNormCode, name = "MixNorm", constants = MixNormConsts,
                       data = MixNormData, inits= MixNormInits)
#
# To compile the model
CMixNorm <- compileNimble(MixNorm)

# set up the monitored quantities. Default is all of the random quantities
MixNormConf <- configureMCMC(MixNorm, monitors = c('pi','mu','z'), print = TRUE)
# binary sampler= gibbs sampler for binary-valued obs
#
# build the MCMC algorithm
MixNormMCMC <- buildMCMC(MixNormConf)
# compile the MCMC chain 
CMixNormMCMC <- compileNimble(MixNormMCMC, project = MixNorm)
#
set.seed(15)
#set.seed(100)
MixNormInits <- list(list(pi = 0.5, mu=c(20,500),z=rbinom(1250,1,0.5)),
                     list(pi=0.5, mu=c(500,20),z=rbinom(1250,1,0.5)))
#
posterior <- runMCMC(CMixNormMCMC, niter = 10000, thin=1, nburnin=1, 
                     summary = TRUE, WAIC = FALSE, samples = TRUE, nchains=2, 
                     samplesAsCodaMCMC=TRUE, inits = MixNormInits) 
#
#
combinedchains <- mcmc.list(posterior$samples$chain1[,c("pi","mu[1]","mu[2]")],
                            posterior$samples$chain2[,c("pi","mu[1]","mu[2]")])
windows()# or quartz() on a MacOS
plot(combinedchains)
windows();gelman.plot(combinedchains)# 
autocorr.plot(posterior$samples$chain1[,c("pi","mu[1]","mu[2]")])
autocorr.plot(posterior$samples$chain2[,c("pi","mu[1]","mu[2]")])
effectiveSize(combinedchains)
posterior$summary$all.chains
summary(combinedchains)# returns Time-series SE
#####################################
# However...
#####################################

set.seed(100)# by just changing the seed..
#
MixNormInits <- list(list(pi = 0.5, mu=c(20,500),z=rbinom(1250,1,0.5)),
                     list(pi=0.5, mu=c(500,20),z=rbinom(1250,1,0.5)))
#
posterior <- runMCMC(CMixNormMCMC, niter = 10000, thin=1, nburnin=1, 
                     summary = TRUE, WAIC = FALSE, samples = TRUE, nchains=2, 
                     samplesAsCodaMCMC=TRUE, inits = MixNormInits) 
#
#
combinedchains <- mcmc.list(posterior$samples$chain1[,c("pi","mu[1]","mu[2]")],
                            posterior$samples$chain2[,c("pi","mu[1]","mu[2]")])
windows()# or quartz() on a MacOS
plot(combinedchains)
#
# The posterior distributions for mu and pi have actually two modes (we saw just
# one mode in the first run -line 60- only by chance). 
# The marginal distributions for mu[1] and mu[2] are identical unless we relabel
# the samples (manually or using postprocessing algorithms) or put constraints 
# to the priors to make the separate mixture components distinguishible (not
# covered in this course) 