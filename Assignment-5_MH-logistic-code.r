#this is an example of using discrete-valued parameters for MH algorithm

#Import the data
DataFile = "C:/Users/SAMSUNG/Downloads/CreditCard_LatePayment_data (1).csv"
LP.data = read.csv(DataFile, header=T)

#Run the logistic regression
Logitmodel <- glm(Latepay ~ Usage + Balance, data = LP.data, family = binomial(link = "logit"))
summary(Logitmodel)
mle_ <- summary(Logitmodel)$coefficients[, "Std. Error"]

#Extract the binary default data column
yy = LP.data[,2]
xx <- data.matrix(cbind(1, LP.data[, c("Usage", "Balance")]))

#Logistic likelihood calculation
YYloglik = function(y, x, b) {
	exb = exp(x%*%b)
	p = exb/(1+exb)
	y.loglik = sum(dbinom(y, 1, p, log=T))
	y.loglik
}

#initial current b (regression coefficients) value 
b.cur = as.numeric(rep(0, ncol(xx)))
#initial likelihood evaluation at the current p value
yy.llik.cur = YYloglik(yy, xx, b.cur)

#posterior 
n.step = 5000
b.post = matrix(0, n.step, dim(xx)[2])

#MH steps
for(gg in 1:n.step) {
	#Sample proposal regression coefficients:
	#Use the current coefficient values b.cur as the mean and 1/10 of the MLE estimates as the SD. Sample the proposal b.cur from a normal distribution  
	#Please blanks below 
	b.prop = rnorm(length(b.cur),mean = b.cur, sd= mle_ / 10)
	
	#evaluate likelihood and decide on posterior sample
	#Evaluate the likelihood at the new proposal p value using the function YYloglik()
	yy.llik.prop = YYloglik(yy, xx, b.prop)
	
	#Use the Metropolis-Hasting algorithm criterion, determine the posterior sample of the regression coefficients 
	#Hint: use the formula in our lecture slide to compute the ratio of the likelihood at the proposal value and the current value
	#Then determine whether the proposal value or the current value should be accepted as the posterior sample
	alpha_ = exp(min(0, yy.llik.prop - yy.llik.cur))  # MH acceptance ratio
	
	if (runif(1) < alpha_) {
	  b.cur = b.prop  # Accept the new proposal
	  yy.llik.cur = yy.llik.prop  # Update the current likelihood
	}
	
	b.post[gg,] = b.cur
}

#plot the posterior and calculate the posterior mean
plot(b.post[,1], type="l")
plot(b.post[,2], type="l")
plot(b.post[,3], type="l")

hist(b.post[1001:5000,1])
hist(b.post[1001:5000,2])
hist(b.post[1001:5000,3])

