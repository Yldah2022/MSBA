
library(mvtnorm)
Y = rp.data$logPurchase
X = cbind(1, as.matrix(rp.data[,c("Age","logPrice","Feature")]))
X2 = t(X)%*%X
n.obs = length(Y)
#Initial Setup for the algorithm
NIT = 1020       #num of interations
nBurn = 20      #num of burn-ins  

#stage 3. Record Posterior samples
beta.dim = 4
beta.pos = matrix(0, NIT, beta.dim)
tau.pos = rep(0, NIT)
  
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
curBeta = c(1, 0, 0, 0) 
curTau = 4

#main loop
for (m in 1:NIT){
	#step 1 sample beta
	beta.pos.var = solve(curTau*X2 + iSigma.beta)
	beta.pos.mean = beta.pos.var%*%(curTau*t(X)%*%Y + iSigma.beta%*%mu.beta)
	curBeta = as.vector(rmvnorm(1, mean=beta.pos.mean, sigma=beta.pos.var)) 
	
	#step 2 sample tau (precision = 1/sigma^2)
	Y.res = as.vector(Y-X%*%curBeta)
	curTau = rgamma(1, a.tau+0.5*n.obs, b.tau+0.5*(t(Y.res)%*%Y.res))
	
	if(m > nBurn){
		beta.pos[m-nBurn,] = curBeta
		tau.pos[m-nBurn] = curTau
	}
}

ts.plot(beta.pos)
ts.plot(tau.pos)

hist(beta.pos[,3], freq=F, nclass=30)
lines(density(beta.pos[,3]))

hist(1/tau.pos, freq=F, nclass=30)
lines(density(1/tau.pos))

colMeans(beta.pos)
mean(sqrt(1/tau.pos))



