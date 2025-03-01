################## STAT0011 ICA: Copula Modelling for VaR ######################
#                                                                              #
#                                   Group 5                                    #
#                                                                              #
#     THIS IS A DRAFT OF STAT0011 GROUP ICA; DO NOT USE IT FOR SUBMISSION!     #
#                                                                              #
#                              INTERNAL USE ONLY                               #
#                                                                              #
#                                                                              #
################################################################################


################################### Set up #####################################

## Install and load the packages required.

library(stats)
install.packages("VineCopula")
library(VineCopula)
install.packages("fGarch")
library(fGarch)
require(fBasics)
install.packages("KScorrect")
library(KScorrect)
install.packages("ADGofTest")
library(ADGofTest)


## Import the weekly financial data for S&P500 and SSE Composite Index, each
## including information about 1214 trading weeks recorded.

SP500 <- read.csv("SP500.csv")
SSE <- read.csv("SSE.csv")

# For coding convenience, all codes related to S&P500 are now labelled with 
# the number '1' (e.g. ret1, model1), and those related to SSE Composite Index 
# are labelled with the number '2' (e.g. ret2, model2) from now on.

## Now compute the weekly log-returns using the weekly prices.

price1 <- rev(SP500$Price)
price2 <- rev(SSE$Price)

ret1 <- log(price1[-1]/price1[-1214])
ret2 <- log(price2[-1]/price2[-1214])


####################### Modelling of log-returns and PIT #######################

## Identification

plot(ret1, col="blue", main = "S&P500", ylab = "log-returns", 
     xlab = "trading week", type = "l")
plot(ret2, col="blue", main = "SSE Composite Index", ylab = "log-returns", 
     xlab = "trading week", type = "l")

jarqueberaTest(ret1)
jarqueberaTest(ret2)

# Both plots exhibit an evident feature of volatility clustering. Meanwhile,
# the p-values for both Jarque-Bera tests are significantly smaller than 0.05, 
# which suggests that the log-returns of neither investment likely follow a 
# normal distribution. Subsequently, the autocorrelation effects in these time 
# series data are examined to fit time series models (autoregressive models),
# with ARCH/GARCH models accounting for time-varying conditional variance.

par(mfrow=c(2,2))
acf(ret1, col="green", lwd=2)   
pacf(ret1, col="green", lwd=2)  
acf(ret1^2, col="red", lwd=2)  
par(mfrow=c(1,1))

# For S&P500, the sample PACF cuts to 0 after lag 1. Thus, an AR(1) model is 
# fitted. Since volatility clustering is detected in the sample ACF of
# log-returns squared, as well as in the previous plot, GARCH is introduced.

par(mfrow=c(2,2))
acf(ret2, col="green", lwd=2)   
pacf(ret2, col="green", lwd=2)  
acf(ret2^2, col="red", lwd=2)    
par(mfrow=c(1,1))

# For SSE Composite Index, the sample PACF approximately cuts to 0 after lag 2. 
# Thus, an AR(2) model is fitted. Since volatility clustering is detected in 
# the sample ACF of log-returns squared, as well as in the previous plot, 
# GARCH is introduced.

## Estimation and Model Checking

# With the initial assumption that the conditional distribution of errors is
# normal, AR(1) models are applied to the weekly log-returns of S&P500, where
# the order 1 is specified by previous examination of the PACF plot, with an 
# ARCH or GARCH component accounting for volatility clustering.

model1_norm_1=garchFit(formula=~arma(1,0)+garch(1,0),data=ret1,trace=F,
                cond.dist="norm")
model1_norm_1
model1_norm_2=garchFit(formula=~arma(1,0)+garch(1,1),data=ret1,trace=F,
                       cond.dist="norm")
model1_norm_2

# Various checking methods are then applied to examine the independency and
# normality of the model errors, as well as the model fit. Since normality is
# assumed, AIC and BIC give a valid evaluation of the model fit.

model1_norm_1@fit$ics
model1_norm_2@fit$ics

# The second model achieves a better fit. We then check if model assumptions
# are satisfied.

res1_norm_1 <- residuals(model1_norm_1, standardize=TRUE)
par(mfrow=c(2,1))
acf(res1_norm_1, col="green", lwd=2)
acf(res1_norm_1^2, col="red", lwd=2)
par(mfrow=c(1,1))
Box.test(res1_norm_1, lag = 10, type = c("Ljung-Box"), fitdf = 1)
Box.test(res1_norm_1^2, lag = 10, type = c("Ljung-Box"), fitdf = 1)

# Although the ACF plots suggest that no autocorrelation seems to exist in the 
# residuals, both Ljung-Box tests suggest that there are statistically 
# significant autocorrelation effects between residuals, which implies that
# the model assumption of error independency is violated. For further testing,
# PIT, Kolmogorov-Smirnov test and Anderson-Darling test are implemented.

u1_norm_1<-pnorm(res1_norm_1) 
hist(u1_norm_1, freq=F, main = "Histogram of u1 (normality assumed)",
     xlab = "u1")
#Kolmogorov-Smirnov test
KStest1_norm_1<-LcKS(u1_norm_1, cdf = "punif")
KStest1_norm_1$p.value
#Anderson-Darling test
ADtest1_norm_1<-ad.test(u1_norm_1, null="punif")
ADtest1_norm_1$p.value 

# The histogram does not display the shape of a U(0,1) distribution. Besides, 
# the p-values of both the Kolmogorov-Smirnov and the Anderson-Darling tests
# are much smaller than 0.05. We reject that the transformed residuals are from 
# a U(0,1) distribution.

# Similarly, check for the second model for S&P500 with normality assumption 
# and we obtain the same conclusion that the assumptions of error independency 
# and conditional normality are violated.

res1_norm_2 <- residuals(model1_norm_2, standardize=TRUE)
par(mfrow=c(2,1))
acf(res1_norm_2, col="green", lwd=2)
acf(res1_norm_2^2, col="red", lwd=2)
par(mfrow=c(1,1))
Box.test(res1_norm_2, lag = 10, type = c("Ljung-Box"), fitdf = 1)
Box.test(res1_norm_2^2, lag = 10, type = c("Ljung-Box"), fitdf = 1)

u1_norm_2<-pnorm(res1_norm_2) 
hist(u1_norm_2, freq=F, main = "Histogram of u1 (normality assumed)",
     xlab = "u1")
#Kolmogorov-Smirnov test
KStest1_norm_2<-LcKS(u1_norm_2, cdf = "punif")
KStest1_norm_2$p.value
#Anderson-Darling test
ADtest1_norm_2<-ad.test(u1_norm_2, null="punif")
ADtest1_norm_2$p.value 

# New conditional distributions of errors are thus attempted, and ultimately, 
# the following model for S&P500 passes most tests for error autocorrelation 
# and post-PIT uniformity and achieves the largest p-value for the rejected
# test:

model1=garchFit(formula=~arma(1,0)+garch(1,1),data=ret1,trace=F, 
                cond.dist="sstd")
model1

# Since the skewness is rather limited, and a student-t distribution tends to
# approach a normal distribution as the degree of freedom gets large, we know
# that AIC and BIC are still valid criteria to evaluate the model fit.

model1@fit$ics

# The AIC and BIC of this model are the greatest among all models with the
# assumed error distribution being normal, skewed normal, student-t and skewed
# student-t. Outcomes of model checking and PIT of residuals for copula 
# modelling are shown below:

# Error independence
res1 <- residuals(model1, standardize=TRUE)
par(mfrow=c(2,1))
acf(res1, col="green", lwd=2)
acf(res1^2, col="red", lwd=2)
par(mfrow=c(1,1))
Box.test(res1, lag = 10, type = c("Ljung-Box"), fitdf = 1)
Box.test(res1^2, lag = 10, type = c("Ljung-Box"), fitdf = 1)

# No autocorrelation is observed for any lag greater than 0, and the p-values 
# of both Ljung-Box tests are significantly greater than 0.05, which suggests 
# that the model assumption of independent error terms is satisfied.

# PIT and uniformity
u1<-psstd(res1, nu = model1@fit$par["shape"], xi = model1@fit$par["skew"])
hist(u1, freq=F)
#Kolmogorov-Smirnov test
KStest1<-LcKS(u1, cdf = "punif")
KStest1$p.value
#Anderson-Darling test
ADtest1<-ad.test(u1, null="punif")
ADtest1$p.value # does not pass the test but fairly close to 0.05

# The histogram roughly displays the shape of a U(0,1) distribution. Besides, 
# the p-value of the Kolmogorov-Smirnov test is above 0.05, and although the
# p-value of the Anderson-Darling test is smaller than 0.05, it is still close. 
# This suggests that it can unlikely be rejected that the samples are from a 
# U(0,1) distribution.

# The identical procedure is applied to SSE Composite Index with an ultimate 
# AR(2)-GARCH(1,1) model found to achieve a good fit, pass all tests for error 
# autocorrelation and post-PIT uniformity. The details of this model are shown
# below:

model2=garchFit(formula=~arma(2,0)+garch(1,1),data=ret2,trace=F,
                cond.dist="sstd")
model2

# Information criteria
model2@fit$ics

# Error independence
res2 <- residuals(model2, standardize=TRUE)
par(mfrow=c(2,1))
acf(res2, col="green", lwd=2)
acf(res2^2, col="red", lwd=2)
par(mfrow=c(1,1))
Box.test(res2, lag = 10, type = c("Ljung-Box"), fitdf = 2)
Box.test(res2^2, lag = 10, type = c("Ljung-Box"), fitdf = 2)

# PIT and uniformity
u2<-psstd(res2, nu = model2@fit$par["shape"], xi = model2@fit$par["skew"])
hist(u2, freq=F)
#Kolmogorov-Smirnov test
KStest2<-LcKS(u2, cdf = "punif")
KStest2$p.value
#Anderson-Darling test
ADtest2<-ad.test(u2, null="punif")
ADtest2$p.value


############################## Copula modelling ################################

model=BiCopSelect(u1, u2, familyset=NA, selectioncrit="AIC", indeptest=TRUE, 
                  level=0.05,se = TRUE)
model


################## Value-at-Risk using Monte Carlo simulation ##################

# Use the explicit method to generate 100000 possible realisations of the
# 1214th error term from the copula model for each of the investments, 
# in order for the simulation of their respective 1214th weekly log-return. 
# As in both AR(1)-GARCH(1,1) and AR(2)-GARCH(1,1) models, the conditional 
# distribution is a skewed student-t, the realisations of the one-step-ahead 
# error term for both investments are obtained by computing the quantiles of 
# the skewed student-t distribution using the simulated u1 and u2.

N=100000
set.seed(0123) # to make all simulated data consistent
u_sim=BiCopSim(N, family=model$family, model$par, model$par2)

res1_sim <- qsstd(u_sim[,1], nu = model1@fit$par["shape"], 
                  xi = model1@fit$par["skew"])
res2_sim <- qsstd(u_sim[,2], nu = model2@fit$par["shape"],
                  xi = model2@fit$par["skew"])

# Simulate the log-return one week ahead for S&P500:

# Now we reintroduce GARCH(1,1) effects to S&P500, starting from the estimated
# unconditional variance, sigma squared 0 = omegahat/(1-alphahat-betahat), and 
# use residuals of model 1 (which are estimates of the 1st to 1213th error 
# terms) to estimate/update the conditional variances sigma squared 1 hat to 
# sigma squared 1213 hat. The 1214th conditional variance then acquires 100000 
# possible realisations, computed via realisations of the 1214th error term, 
# epsilon_1^(1) to epsilon epsilon_1^(100000), as well as the estimated 1213th 
# conditional variance and other model parameters estimated. 

t <- length(ret1)
sp_sigma2hat <- numeric(t+1)

sp_omegahat <- as.numeric(model1@fit$par["omega"])
sp_alphahat <- as.numeric(model1@fit$par["alpha1"])
sp_betahat <- as.numeric(model1@fit$par["beta1"])
sp_errorhat <- model1@residuals

sp_sigma2hat[1] <- sp_omegahat/(1-sp_alphahat-sp_betahat) # sigma squared 0

for(i in 2:(t+1)){
  sp_sigma2hat[i] <- sp_omegahat + 
    sp_alphahat*sp_sigma2hat[i-1]*(sp_errorhat[i-1]^2) + 
    sp_betahat*sp_sigma2hat[i-1]
} # each conditional variance (1st to 1213th) is estimated 
  # using its predecessor. 

sp_sigma2_sim <- sp_omegahat + sp_alphahat*sp_sigma2hat[t+1]*(res1_sim^2) + 
  sp_betahat*sp_sigma2hat[t+1]

# Next, we reintroduce AR(1) autocorrelation effects for S&P500 using the 
# simulated conditional variance and error term, estimated model parameters 
# and the 1213th weekly log-return of S&P500. Consequently, 100000 realisations 
# of the log-return one week ahead are simulated.

sp_chat <- as.numeric(model1@fit$par["mu"])
sp_aralphahat <- as.numeric(model1@fit$par["ar1"])

ret1_sim <- sp_chat + sp_aralphahat*ret1[length(ret1)] + 
  sqrt(sp_sigma2_sim)*res1_sim

# Similarly, reintroduce the GARCH(1,1) and AR(2) autocorrelation effects to
# simulate the weekly log-return for SSE Composite Index:

sse_sigma2hat <- numeric(t+1)

sse_omegahat <- as.numeric(model2@fit$par["omega"])
sse_alphahat <- as.numeric(model2@fit$par["alpha1"])
sse_betahat <- as.numeric(model2@fit$par["beta1"])
sse_errorhat <- model2@residuals

sse_sigma2hat[1] <- sse_omegahat/(1-sse_alphahat-sse_betahat)

for(i in 2:(t+1)){
  sse_sigma2hat[i] <- sse_omegahat + 
    sse_alphahat*sse_sigma2hat[i-1]*(sse_errorhat[i-1]^2) + 
    sse_betahat*sse_sigma2hat[i-1]
}

sse_sigma2_sim <- sse_omegahat + sse_alphahat*sse_sigma2hat[t+1]*(res2_sim^2) + 
  sse_betahat*sse_sigma2hat[t+1]

sse_chat <- as.numeric(model2@fit$par["mu"])
sse_aralpha1hat <- as.numeric(model2@fit$par["ar1"])
sse_aralpha2hat <- as.numeric(model2@fit$par["ar2"])

ret2_sim <- sse_chat + sse_aralpha1hat*ret2[length(ret2)] + 
  sse_aralpha2hat*ret2[length(ret2)-1] + sqrt(sse_sigma2_sim)*res2_sim

# Finally, estimate the 95% and 99% 1-Week Value-at-Risk of the equally 
# weighted portfolio according to the simulation:

port_sim <- log(1+((exp(ret1_sim)-1)+(exp(ret2_sim)-1))*(1/2))
negvar_sim <- quantile(port_sim,c(0.01,0.05))

var_sim <- setNames(c(-as.numeric(negvar_sim[2]), -as.numeric(negvar_sim[1])), 
                      c("95%", "99%"))
var_sim

# Thus the estimated 95% 1-Week VaR of the portfolio is 1.53% and the
# estimated 99% 1-Week VaR is 3.24% according to results stored by the seed 
# 0123, visualised in a histogram of simulated weekly portfolio loss.

hist(-port_sim, breaks = 250, freq = FALSE,
     main = "Distribution of simulated weekly portfolio loss",
     xlab = "loss",col = "lightblue", ylim = c(0, 50))
lines(density(-port_sim), col = "darkblue", lwd = 2)
abline(v = var_sim[1], col = "red",  lwd = 2, lty = 2)
abline(v = var_sim[2], col = "blue", lwd = 2, lty = 2)
legend("topleft",legend = c("95% VaR","99% VaR"),
       col=c("red","blue"),lty = 2, lwd= 2, bty = "n")

