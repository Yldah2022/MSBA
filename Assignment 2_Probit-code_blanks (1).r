#4
library(readr)

LP.data <- read.csv("C:/Users/SAMSUNG/Downloads/CreditCard_LatePayment_data.csv")

# Fit the probit regression model by MLE
probit_LP <- glm(Latepay ~ Usage + Balance, 
                    data = LP.data, 
                    family = binomial(link = "probit"))

# View the summary of the model
summary(probit_LP)


#5
library(truncnorm)
library(mnormt)
#Bayesian estimation for probit regression
#stage 1. read data into R and create columns for censored data


#stage 1. subset the data for Latepay = 1 and =0 
LP.X0 = cbind(1, as.matrix(LP.data[LP.data$Latepay==0, 3:4]))
LP.X1 = cbind(1, as.matrix(LP.data[LP.data$Latepay==1, 3:4])) 
LP.X = cbind(1, as.matrix(LP.data[, 3:4]))
LP.X2 = t(LP.X)%*%LP.X
n0 = dim(LP.X0)[1]
n1 = dim(LP.X2)[1]
nObs = dim(LP.data)[1]
LP.Y = rep(0, nObs)

#stage 2. Initial Setup for the algorithm
NIT = 10000       #num of interations
nBurn = 2000      #num of burn-ins  
NIT.eff = NIT - nBurn    #effective sample size
thin.step = 10           #thinning  
NIT.thin = floor(NIT.eff/thin.step)   #effective sample size after thinning

#stage 3. Record Posterior samples
beta.dim = 3
beta.pos = matrix(0, NIT.thin, beta.dim)
  
#stage 4. priors
#for Beta: mNormal(mu.beta, sigma.beta)
mu.beta = rep(0,beta.dim) 
sigma.beta = 1E6 * diag(beta.dim)  
iSigma.beta = 1E-6 * diag(beta.dim)  #inverse prior covariance matrix 

#stage 5. Gibbs sampler

#initialize the loop
curBeta = c(0.1, 0, 0) #initial (current) regression coeff beta
g = 1

#main loop
for (m in 1:NIT){
	#step 1. sample the latent variable > 0 if Latepay=1, <0 if Latepay=0 
	#use the corresponding truncated normal distribution given curbeta and X variables
	#Please fill in the code 
  # For Latepay = 0
  LP.Y[LP.data$Latepay == 0] <- rtruncnorm(n0, a = -Inf, b = 0, mean = LP.X0 %*% curBeta, sd = 1)
  
  # For Latepay = 1
  LP.Y[LP.data$Latepay == 1] <- rtruncnorm(n1, a = 0, b = Inf, mean = LP.X1 %*% curBeta, sd = 1)
	
	
	#step 2 sample curbeta (same as the linear regression code assuming the error's variance is 1)
	#Please fill in the code 
  # Calculate the posterior mean and variance for beta
  postVarBeta <- solve(LP.X2 + iSigma.beta)
  postMeanBeta <- postVarBeta %*% (t(LP.X) %*% LP.Y + iSigma.beta %*% mu.beta)
 
  curBeta <- mvrnorm(1, mu = postMeanBeta, Sigma = postVarBeta)
	
	
	#save thinned sampled beta after burn-ins
	if ((m > nBurn) & (m%%thin.step == 0)) {
		beta.pos[g,] = curBeta
		g = g+1
	}
}

# Time-series plot for beta0, beta1, beta2
beta0_ts <- ts(beta.pos[, 1])
beta1_ts <- ts(beta.pos[, 2])
beta2_ts <- ts(beta.pos[, 3])

ts.plot(beta0_ts, beta1_ts, beta2_ts, col = c("navy", "red", "green"), 
        main = "Posterior chains for beta coefficients", xlab = "Iteration", ylab = "Beta Coefficients")

legend("topright", inset=.05, legend = c("beta0", "beta1", "beta2"), 
       col = c("navy", "red", "green"), lwd = 1, bty = "n")

# Histograms for beta0, beta1, beta2
hist(beta.pos[, 1], main = "Posterior distribution of beta0", xlab = "beta0")
hist(beta.pos[, 2], main = "Posterior distribution of beta1", xlab = "beta1")
hist(beta.pos[, 3], main = "Posterior distribution of beta2", xlab = "beta2")

# Calculate the 95% posterior intervals for beta0, beta1, beta2
beta0_CI <- quantile(beta.pos[, 1], probs = c(0.025, 0.975))
beta1_CI <- quantile(beta.pos[, 2], probs = c(0.025, 0.975))
beta2_CI <- quantile(beta.pos[, 3], probs = c(0.025, 0.975))

# Print the 95% posterior intervals
print(paste0("95% CI for beta0: [", beta0_CI[1], ", ", beta0_CI[2], "]"))
print(paste0("95% CI for beta1: [", beta1_CI[1], ", ", beta1_CI[2], "]"))
print(paste0("95% CI for beta2: [", beta2_CI[1], ", ", beta2_CI[2], "]"))
