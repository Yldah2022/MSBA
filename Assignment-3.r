#1
library(Matrix)
papertowel.data = read.csv("C:/Users/SAMSUNG/Downloads/Papertowel_repurchase.csv", header=TRUE)
papertowel.data$season = as.factor(ceiling(papertowel.data$week / 13))

interval = c()
for (i in 1:500) {
  papertowel.i = papertowel.data[papertowel.data$consumerID == i,]
  interval.i = rep(0, 52)
  sincePurchase = 0
  for (t in 1:52) {
    sincePurchase = sincePurchase + 1
    interval.i[t] = sincePurchase
    if (papertowel.i$papertowel[t] == 1) sincePurchase = 0
  }
  interval = c(interval, interval.i)
}

papertowel.data$interval = interval

#2
# Fit the logistic regression model
hazard1 <- glm(papertowel ~ interval + season + price + feature + famsize, 
                    family = binomial(link = "logit"), data = papertowel.data)

summary(hazard1)

hazard1_aic <- AIC(hazard1)
print(hazard1_aic)

#3
# Fit the cloglog regression model
cloglogm <- glm(formula = papertowel ~ interval + season + price + feature + famsize,
                     family = binomial(link = "cloglog"), data = papertowel.data)

summary(cloglogm)

clog_aic <- AIC(cloglogm)
print(clog_aic)


#4
library(lme4)
# Fit the cloglog regression model with a random effect for the intercept
cloglog_re_model <- glmer(papertowel ~ interval + season + price + feature + famsize + 
                            (1 | consumerID), 
                          data = papertowel.data, 
                          family = binomial(link = "cloglog"))

summary(cloglog_re_model)

# Calculate the AIC of the model
cloglog_re_model_aic <- AIC(cloglog_re_model)
print(cloglog_re_model_aic)


#2(1)
install.packages("censReg")
library(censReg)
cc_data <- read.csv("C:/Users/SAMSUNG/Downloads/CreditCard_SOW_data2.csv", header = TRUE)

# Fit the censored regression model
cc_creg <- censReg(SOW ~ Balance + Promotion, 
                   left = 0, right = 1, 
                   data = cc_data)

summary(cc_creg)


library(truncnorm)
library(mnormt)
library(MASS)
#Bayesian estimation for truncated regression
#stage 1. read data into R and create columns for censored data
DataFile = "C:/Users/SAMSUNG/Downloads/CreditCard_SOW_data2.csv"
sow.data = read.csv(DataFile, header=T)
sow.data$Cens0 = (sow.data$SOW==0)*1
sow.data$Cens1 = (sow.data$SOW==1)*1

#extract right and left censored data
sow.Y = sow.data$SOW
sow.XLC = cbind(1, as.matrix(sow.data[sow.data$Cens0==1, 3:4]))
sow.XRC = cbind(1, as.matrix(sow.data[sow.data$Cens1==1, 3:4])) 
sow.X = cbind(1, as.matrix(sow.data[, 3:4]))
sow.X2 = t(sow.X)%*%sow.X
nRC = dim(sow.XRC)[1]
nLC = dim(sow.XLC)[1]
nObs = dim(sow.data)[1]

#stage 2. Initial Setup for the algorithm
NIT = 10000       #num of interations
nBurn = 2000      #num of burn-ins  
NIT.eff = NIT - nBurn    #effective sample size
thin.step = 10           #thinning  
NIT.thin = floor(NIT.eff/thin.step)   #effective sample size after thinning

#stage 3. Record Posterior samples
beta.dim = 3
beta.pos = matrix(0, NIT.thin, beta.dim)
tau.pos = rep(0, NIT.thin)

#stage 4. priors
#for Beta: mNormal(mu.beta, sigma.beta)
mu.beta = rep(0,beta.dim) 
sigma.beta = 1E6 * diag(beta.dim)  
iSigma.beta = 1E-6 * diag(beta.dim)  #inverse prior covariance matrix 

#prior for precision: Gamma(a.tau, b.tau)
a.tau = 1/2
b.tau = 1/2

#stage 5. Gibbs sampler

#initialize the loop
curBeta = c(0.5, 0, 0) 
curTau = 4
g = 1

#main loop
for (m in 1:NIT){
  #step 1. sample latent SOW 
  #step 1.a. sample SOW censored at 0
  #Please fill in the blank below the code for sampling the latent SOW when the observed SOW==0
  #Please name your sampled latent SOW as curYLC
  curYLC =  rtruncnorm(nLC, a = 0, b = Inf, mean = sow.XLC %*% curBeta, sd = sqrt(1/curTau))
  
  #step 1.b sample SOW censored at 1
  #Please fill in the blank below the code for sampling the latent SOW when the observed SOW==1
  #Please name your sampled latent SOW curYRC
  curYRC =  rtruncnorm(nRC, a = -Inf, b = 1, mean = sow.XRC %*% curBeta, sd = sqrt(1/curTau))
  
  #step 2 sample beta
  #step 2.a impute the latent variables
  sow.Y[sow.data$Cens0==1] = curYLC
  sow.Y[sow.data$Cens1==1] = curYRC
  
  #step 2.b sample beta
  #Please fill in the blanks below
  sigma.hat =  solve(sow.X2*curTau + iSigma.beta)
  betaPos.mn =  sigma.hat %*% (t(sow.X) %*% (sow.Y*curTau) + iSigma.beta %*% mu.beta)
  curBeta = mvrnorm(1, mu = betaPos.mn, Sigma = sigma.hat)
  
  #step 3 sample tau (precision = 1/sigma^2)
  #Please fill in the blanks below
  sowE.hat =  sow.Y - sow.X %*% curBeta
  a.tau.new = a.tau + nObs/2
  b.tau.new = b.tau + sum(sowE.hat^2)/2
  curTau = rgamma(1, shape = a.tau.new, rate = b.tau.new)
  
  #save thinned samples after burn-ins
  if ((m > nBurn) & (m%%thin.step == 0)) {
    beta.pos[g,] = curBeta
    tau.pos[g] = curTau
    g = g+1
  }
}

par(mfrow = c(2, 2))
plot(beta.pos[, 1], type = 'l', main = expression(beta[0]))
hist(beta.pos[, 1], main = expression(beta[0]), xlab = expression(beta[0]), breaks = 50)
plot(beta.pos[, 2], type = 'l', main = expression(beta[1]))
hist(beta.pos[, 2], main = expression(beta[1]), xlab = expression(beta[1]), breaks = 50)
plot(beta.pos[, 3], type = 'l', main = expression(beta[2]))
hist(beta.pos[, 3], main = expression(beta[2]), xlab = expression(beta[2]), breaks = 50)
plot(tau.pos, type = 'l', main = expression(tau))
hist(tau.pos, main = expression(tau), xlab = expression(tau), breaks =50)

# Calculate the 95% posterior intervals for each beta coefficient
beta0_95CI <- quantile(beta.pos[, 1], probs = c(0.025, 0.975))
beta1_95CI <- quantile(beta.pos[, 2], probs = c(0.025, 0.975))
beta2_95CI <- quantile(beta.pos[, 3], probs = c(0.025, 0.975))

# Calculate the 95% posterior interval for tau
tau_95CI <- quantile(tau.pos, probs = c(0.025, 0.975))

# Print the intervals
cat("95% posterior interval for beta0:", beta0_95CI, "\n")
cat("95% posterior interval for beta1:", beta1_95CI, "\n")
cat("95% posterior interval for beta2:", beta2_95CI, "\n")
cat("95% posterior interval for tau:", tau_95CI, "\n")
