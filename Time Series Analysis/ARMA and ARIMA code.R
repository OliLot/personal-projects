require(tseries)
require(zoo)
require(forecast)

################################################################################
########### KEY COMMANDS ###############
################################################################################
# 1. acf, pacf = estimated ACF/PACF
# 2. ARMAacf (with pacf = False standard, or True) = THEORETICAL ACF/PACF
# 3. arima(x, order=(p,d,q)) = fit an arima model of AR(p), differencing d, MA(q)
# 4. Arima(x, order=(p,d,q)) = wrapper for arima function that allows a drift term
# 5. auto.arima(x, ...) = first best ARIMA model according to AIC, AICc, or BIC - requires search specification
#                         forward selection based (try different ones always keeping best, not thorough but good indicator)
# 6. predict(model fit, h steps) = prediction of h next steps in future (with standard error)
# 7. forecast(model fit, h steps) = forecast of h next steps in future (with confidence intervals) - may be slightly different


?auto.arima


### // Simulate MA(1) theta = 0.8 data of length 101
set.seed(31)
ma1 <- arima.sim(n=101, list(ma=0.8)) #standard R command
?arima.sim


### Plot the simulated ARMA(0, 0.8) model
plot(as.zoo(ma1),type="o", lwd=2, pch=16, ylab="MA (1)", cex.axis=1.7, cex.lab=1.7, cex.main=1.7)
abline(h=0,lty=2, lwd=2) # add a straight line through plot at h = 0 y value.


### Calculate ACF of the simulated data
par(mfrow=c(1,2))

# acf function: estimates of ACVF or ACF (pacf function for PACF, ccf for cross-correlation) with zoo class data
acf(as.zoo(ma1), ylim=c(-1,1), main= "MA(1)", cex.axis=1.7, cex.lab=1.7, cex.main=1.7)

# ARMAacf function: theoretical ACF values for a given ARMA process
lines(0:20,ARMAacf(ma=0.8,lag.max=20), type="p", col="blue",lwd=3, pch=16)

# adding a legend
legend(0,-0.5,
       legend=c("Sample ACF", "Theoretical ACF"), lty=c(1,NA), pch=c(NA,16),
       col=c("black","blue"), cex = 1 )

# implementing pacf function for estimate PACF of simulated MA(1) 0.8 process
pacf(as.zoo(ma1), ylim = c(-1,1), main= "MA(1)", cex.axis=1.7, cex.lab=1.7, cex.main=1.7)
# comparing to theoretical pacf values
lines(1:20, ARMAacf(ma=0.8,lag.max=20, pacf=TRUE), type="p", col="blue",lwd=3, pch=16)
legend(0, -0.5,
       legend=c("Sample PACF", "Theoretical PACF"), lty=c(1,NA), pch=c(NA,16),
       col=c("black","blue"), cex = 1 )
par(mfrow=c(1,1))


#-------------------------------------------------------------------------------


### // Simulate  MA(1) theta = -0.8 data of length 101 
# (same estimate vs theoretical comparisons below)
set.seed(13)
ma2=arima.sim(n=101, list(ma=-0.8))

plot(as.zoo(ma2),type="o", lwd=2, pch=16, ylab="MA (1)", cex.axis=1.7, cex.lab=1.7, cex.main=1.7)
abline(h=0,lty=2, lwd=2)

par(mfrow=c(1,2))
acf(ma2, ylim=c(-1,1), main = expression(paste("MA(1), ", theta, "=-0.8")), cex.axis=1.7, cex.lab=1.7, cex.main=1.7)
lines(0:20,ARMAacf(ma=-0.8,lag.max=20), type="p", col="blue",lwd=3, pch=16)
legend(10.7,-0.5,
       legend=c("Sample ACF", "Theoretical ACF"), lty=c(1,NA), pch=c(NA,16),
       col=c("black","blue"), cex = 1.5 )
pacf(ma2, ylim = c(-1,1), main = expression(paste("MA(1), ", theta, "=-0.8")),  cex.axis=1.7, cex.lab=1.7, cex.main=1.7)
lines(1:20, ARMAacf(ma=-0.8,lag.max=20, pacf=TRUE), type="p", col="blue",lwd=3, pch=16)
legend(10.7,-0.5,
       legend=c("Sample PACF", "Theoretical PACF"), lty=c(1,NA), pch=c(NA,16),
       col=c("black","blue"), cex = 1.5 )
par(mfrow=c(1,1))




################################################################################
################################################################################
################################################################################
# ACF and PACF of a simulated AR(1) model

ar1a=function(k,h)  k^h

### // Simulate AR(1) phi = 0.8 data of length 101
set.seed(22173)
ar1=arima.sim(n=101,list(ar=0.8))

### first 21 coefficients phi^(0->20) -> theoretical ACF values (lags 0-21)
aa1=c(1,sapply(1:20, function(x) ar1a(0.8,x)))
### 21 coefficients (0,1 -> 19 rest are zero) -> theoretical PACF values (lags 0-21)
pa1=c(1,0.8,rep(0,19))

### Plot simulated AR(1) 0.8 data
plot(as.zoo(ar1),type="o", lwd=2, pch=16, ylab="AR(1)", xlab="Time", cex.axis=1.7, cex.lab=1.7, cex.main=1.7)
abline(h=0,lty=2, lwd=2)

par(mfrow=c(1,2))
### Calculate ACF and PACF estimates
acf(ar1, ylim=c(-1,1), cex.axis=1.7, cex.lab=1.7, cex.main=1.7)
lines(0:20,aa1, type="p", col="blue",lwd=3, pch=16)
legend(12,-0.5,
       legend=c("Sample ACF", "Theoretical ACF"), lty=c(1,NA), pch=c(NA,16),
       col=c("black","blue"), cex = 1.5 )
pacf(ar1, ylim = c(-1,1), cex.axis=1.7, cex.lab=1.7, cex.main=1.7)
lines(0:20,pa1, type="p", col="blue",lwd=3, pch=16)
legend(12,-0.5,
       legend=c("Sample PACF", "Theoretical PACF"), lty=c(1,NA), pch=c(NA,16),
       col=c("black","blue"), cex = 1.5 )

par(mfrow=c(1,1))


# ------------------------------------------------------------------------------


### // Simulate data from AR(1) phi = -0.8 data of length 101
# same comparisons of ACF and PACF estimates to theoretical values derived from
# the coefficient powers
set.seed(22173)
ar2=arima.sim(n=101, list(ar=-0.8))

aa2=c(1,sapply(1:20, function(x) ar1a(-0.8,x)))
pa2=c(1,-0.8,rep(0,19))

plot(as.zoo(ar2),type="o", lwd=2, pch=16, ylab="AR(1)", xlab="Time", cex.axis=1.7, cex.lab=1.7, cex.main=1.7)
abline(h=0,lty=2, lwd=2)

par(mfrow=c(1,2))
acf(ar2, ylim=c(-1,1), cex.axis=1.7, cex.lab=1.7, cex.main=1.7)
lines(0:20,aa2, type="p", col="blue",lwd=3, pch=16)
legend(12,-0.5,
       legend=c("Sample ACF", "Theoretical ACF"), lty=c(1,NA), pch=c(NA,16),
       col=c("black","blue"), cex = 1.5 )
pacf(ar2, ylim = c(-1,1), cex.axis=1.7, cex.lab=1.7, cex.main=1.7)
lines(0:20,pa2, type="p", col="blue",lwd=3, pch=16)
legend(12,-0.5,
       legend=c("Sample PACF", "Theoretical PACF"), lty=c(1,NA), pch=c(NA,16),
       col=c("black","blue"), cex = 1.5 )

par(mfrow=c(1,1))


#-------------------------------------------------------------------------------


### // Theoretical ACFs of 4 given AR(2) processes
par(mfrow=c(2,2))
plot(ARMAacf(ar=c(1/6,1/6),lag.max = 15),type="h",ylim=c(-0.35,1), ylab="ACF", xlab="Lag",main=expression(paste(phi[1], " = 1/6, ",phi[2], " = 1/6" )), cex.axis=1.7, cex.lab=1.7, cex.main=1.7)
abline(h=0)
plot(ARMAacf(ar=c(-1/6,1/6),lag.max = 15),type="h",ylim=c(-0.35,1), ylab="ACF", xlab="Lag",main=expression(paste(phi[1], " = -1/6, ",phi[2], " = 1/6" )), cex.axis=1.7, cex.lab=1.7, cex.main=1.7)
abline(h=0)
plot(ARMAacf(ar=c(0.2,0.35),lag.max = 15),type="h",ylim=c(-0.35,1), ylab="ACF", xlab="Lag",main=expression(paste(phi[1], " = 0.2, ",phi[2], " = 0.35" )), cex.axis=1.7, cex.lab=1.7, cex.main=1.7)
abline(h=0)
plot(ARMAacf(ar=c(0.6,-0.4),lag.max = 15),type="h",ylim=c(-0.35,1), ylab="ACF", xlab="Lag",main=expression(paste(phi[1], " = 0.6, ",phi[2], " = -0.4" )), cex.axis=1.7, cex.lab=1.7, cex.main=1.7)
abline(h=0)
par(mfrow=c(1,1))

################################################################################
################################################################################
################################################################################
### // ARMA(1,1) phi = 0.6 & theta = 0.8 simulation
set.seed(91816)
arma11=arima.sim(n=101,list(ar=0.6,ma= 0.8))

plot(as.zoo(arma11),type="o", lwd=2, pch=16, ylab="ARMA(1,1)", cex.axis=1.7, cex.lab=1.7)
abline(h=0,lty=2, lwd=2)

# fn to compute theoretical acf directly from formula on slides
# acf11<-function(a,b,h) {   
#   # a is phi
#   # b is theta
#   # h is the lag
#   erg<-(1+a*b)*(a+b)/(1+2*a*b+b^2)
#   for (i in 2:h) {
#     r<- a^{i-1}*erg[1]
#     erg<-c(erg,r)
#   }
#   return(erg) 
# }
# alternatively use the ARMAacf function

acft <-ARMAacf(0.6,0.8,lag.max=20)
pacft = c(1,ARMAacf(0.6,0.8,lag.max=20,pacf=TRUE))

par(mfrow=c(1,2))
acf(arma11,ylim=c(-1,1), cex.axis=1.7, cex.lab=1.7, cex.main=1.7)
lines(0:20,acft, type="p", col="blue",lwd=3, pch=16)
legend(2,-0.5,
       legend=c("Sample ACF", "Theoretical ACF"), lty=c(1,NA), pch=c(NA,16),
       col=c("black","blue"), cex = 1 )
pacf(arma11,ylim=c(-1,1), cex.axis=1.7, cex.lab=1.7, cex.main=1.7)
lines(0:20,pacft, type="p", col="blue",lwd=3, pch=16)
legend(2,-0.5,
       legend=c("Sample PACF", "Theoretical PACF"), lty=c(1,NA), pch=c(NA,16),
       col=c("black","blue"), cex = 1 )
par(mfrow=c(1,1))

set.seed(991816)
armaM11=arima.sim(n=101,list(ar=-0.7,ma= -0.6))

plot(as.zoo(armaM11),type="o", lwd=2, pch=16, ylab="ARMA(1,1)", cex.axis=1.7, cex.lab=1.7, cex.main=1.7)
abline(h=0,lty=2, lwd=2)

par(mfrow=c(1,2))
acf(armaM11,ylim=c(-1,1), cex.axis=1.7, cex.lab=1.7, cex.main=1.7)
lines(0:20, ARMAacf(-0.7,-0.6,lag.max=20), type="p", col="blue",lwd=3, pch=16)
legend(2.7,-0.5,
       legend=c("Sample ACF", "Theoretical ACF"), lty=c(1,NA), pch=c(NA,16),
       col=c("black","blue"), cex = 1.5 )

pacf(armaM11, ylim=c(-1,1), cex.axis=1.7, cex.lab=1.7, cex.main=1.7)
lines(0:20,c(1,ARMAacf(-0.7,-0.6,lag.max=20,pacf=TRUE)), type="p", col="blue",lwd=3, pch=16)
legend(2.7,-0.5,
       legend=c("Sample PACF", "Theoretical PACF"), lty=c(1,NA), pch=c(NA,16),
       col=c("black","blue"), cex = 1 )
par(mfrow=c(1,1))
################################################################################
##### Fitting a model to simulated data ########################################
################################################################################

xt = arima.sim(list(ar=c(0.5,-0.6),ma=c(0.6,-0.8)),n=500)
### FITING AN ARIMA model to the simulated data (order = p, d, q) in AR(p), differencing degree, MA(q)
arima(xt,order=c(2,0,2))
plot(as.zoo(xt),type="o", lwd=2, pch=16, ylab="ARMA(1,1)", cex.axis=1.7, cex.lab=1.7, cex.main=1.7)
abline(h=0,lty=2, lwd=2)
?arima
xt = 5+arima.sim(list(ar=c(0.5,-0.6),ma=c(0.6,-0.8)),n=500)
arima(xt,order=c(2,0,2),include.mean = TRUE)

#############################################################
## CONSOLIDATION OF R COMMANDS FOR ARMA MODELS         ######
#############################################################
require(tseries) # required for the jarque.bera.test() function - normality test
# plus other ARIMA related tests
require(zoo) # sometimes useful for plotting time series
require(forecast) # for obtaining AICc and the BIC for ARIMA models
# for performing 'automatic' model selection
# and for forecasting after fitting a model

# To obtain the theoretical autocor. and partial autocor. for an ARMA model
# Does not require an R package
ARMAacf(ar=c(1/6,1/6),ma=c(1/2,1/5),lag.max = 15)
plot(ARMAacf(ar=c(1/6,1/6),lag.max = 15),type="h",ylim=c(-0.35,1), ylab="ACF", xlab="Lag",main=expression(paste(phi[1], " = 1/6, ",phi[2], " = 1/6" )), cex.axis=1.7, cex.lab=1.7, cex.main=1.7)
ARMAacf(ar=c(1/6,1/6),ma=c(1/2,1/5),lag.max = 15,pacf='true')

# To simulate observations from an ARMA (or ARIMA) process
# Does not require an R package
arma_2_2 <- arima.sim(n=1001,list(ar=c(1/6,1/6),ma=c(1/2,1/5)))
plot(arma_2_2,type="o", lwd=3, xlab="ARMA_2_2", ylab = "Observation", cex.axis=1.7, cex.lab=1.7, cex.main=1.7)
# check if the sample auto cor. are similar to the theoretical ones
acf(arma_2_2,ylim=c(-0.4,1), ci.type="ma", cex.axis=1.7, cex.lab=1.7, cex.main=1.7)
sample_autocor <- acf(arma_2_2,ylim=c(-0.4,1), ci.type="ma", cex.axis=1.7, cex.lab=1.7, cex.main=1.7)



################################################################################
# For a real data analysis
dat2<-read.table("dgnp82.txt")
# get data vector
gnp<-as.vector(dat2[[1]])
plot(gnp,type="o", lwd=3, ylab="GNP", xlab = "Observation number", cex.axis=1.7, cex.lab=1.7, cex.main=1.7)
abline(h=0, lty=2, lwd=3)

# acf and pacf for sample autocor. and partial autocor.
# do not require an R package
par(mfrow=c(1,2))
acf(gnp,ylim=c(-0.4,1), ci.type="ma", cex.axis=1.7, cex.lab=1.7, cex.main=1.7)
pacf(gnp, ylim=c(-0.4,1), ci.type="ma", cex.axis=1.7, cex.lab=1.7, cex.main=1.7)
par(mfrow=c(1,1))
# to obtain the actual numbers
sample_autocor <- acf(gnp,ylim=c(-0.4,1), ci.type="ma", cex.axis=1.7, cex.lab=1.7, cex.main=1.7)
sample_autocor

# To fit an ARMA (or ARIMA model) to data
# Does not require an R package
arima(gnp,order=c(2,0,2), include.mean= FALSE)
arima(gnp,order=c(2,0,2),include.mean = TRUE) # default
gnp_fit <- arima(gnp,order=c(2,0,2),include.mean = TRUE) # default
# check gnp_fit object. For example, 
# gnp_fit$aic, gnp_estimates$coef  
gnp_fit$residuals

# Alternative to fit an arma (or arima model) to data
# 'forecast' R package required - Arima with capital A
gnp_fit2 <- Arima(gnp,order=c(2,0,2),include.mean=TRUE)
gnp_fit2

tabaicc<-matrix(NA,5,5)
tabbic<-matrix(NA,5,5)

for(i in 0:4){
  for(j in 0:4){
    fit<-Arima(gnp,order=c(i,0,j),include.mean=TRUE)
    tabaicc[(i+1),(j+1)]<-fit$aicc
    tabbic[(i+1),(j+1)]<-fit$bic
  }
}

knitr::kable(
  cbind(0:4,tabbic), booktabs = TRUE, col.names = c("p/q",0,1,2,3,4)
)

knitr::kable(
  cbind(0:4,tabaicc), booktabs = TRUE, col.names = c("p/q",0,1,2,3,4)
)

# Heads-up!! There is an automatic way of finding the 'best' model using
# the auto.arima() function from the 'forecasting' package
# (This may result in an ARIMA model than an ARMA model, so we will 
# look at this function in more detail in Chapter 5.
# Caution! Function admits many more important arguments!!)
auto.arima(gnp, ic='aic')
auto.arima(gnp, ic='aicc')
auto.arima(gnp, ic='bic')

################################################################################
#### PREDICTING/FORECASTING FROM SELECTED MODEL FIT ####
################################################################################


# Forecasting after fitting a model with the predict() function
# Does not require an R package
gnp_fit <- arima(gnp,order=c(2,0,2),include.mean = TRUE)
predict(gnp_fit,n.ahead=10)

# Forecasting after fitting a model with the forecast() function
# Requires the 'forecasting' R package
gnp_fit2 <- Arima(gnp,order=c(2,0,2),include.mean = TRUE)
forecast(gnp_fit2,h=10)