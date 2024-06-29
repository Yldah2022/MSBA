#this is an example of using discrete-valued parameters for MH algorithm

#Import the data
DataFile = "C:/Users/SAMSUNG/Downloads/CreditCard_LatePayment_data (1).csv"
LP.data = read.csv(DataFile, header=T)

#Extract the binary default data column
yy = LP.data[,2]

#estimate pp using MH
#discrete prior 
max.ind = 40
p.priorv = c(0:max.ind)/max.ind
p.priorp = rep(1/(max.ind+1), max.ind+1)

#likelihood calculation
YYloglik = function(y, p) {
	y.loglik = sum(dbinom(y, 1, p, log=T))
	y.loglik
}

#initial current p value 
ind.cur = 11
p.cur = p.priorv[ind.cur] 
#initial likelihood evaluation at the current p value
yy.llik.cur = YYloglik(yy, p.cur)

#posterior 
n.step = 5000
p.post = rep(0, n.step)

#MH steps
for(gg in 1:n.step) {
	#sample proposal parameters
  #For a random walk MH algorithm, please blanks below 
  #When the current p value is 0, the random walk can only jump to the right with probability = 0.5; it stay at 0 with probability = 0.5
  if (ind.cur == 0){
    prop.ind = sample(c(ind.cur,ind.cur+1), 1, prob = c(0.5,0.5))}
  #When the current p value is 1, the random walk can only jump to the left with probability = 0.5; it stay at 1 with probability = 0.5
  else if (ind.cur == max.ind) {
    prop.ind = sample(c(ind.cur-1, ind.cur), 1, prob = c(0.5, 0.5))
  } else {
    prop.ind = sample(c(ind.cur-1, ind.cur+1), size = 1, prob = c(0.5, 0.5))
  }
  
  p.prop = p.priorv[prop.ind]
  
	#Hint: use the index the discrete p values (p.priorv)to simulate a random walk, 
	#with probability = 0.5, the random walk jumps to the right of the current p value, 
	#with probability = 0.5, the random walk jumps to the left of the current p value,

	#When the current p value is 1, the random walk can only jump to the left with probability = 0.5; it stay at 1 with probability = 0.5
	#The new p value the random walk jumps to is the proposal p value
	
	
	
	
	#evaluate likelihood and decide on posterior sample
	#Hint: evaluate the likelihood at the new proposal p value using the function YYloglik()
  yy.llik.prop = YYloglik(yy, p.prop)
	
	#Use the Metropolis-Hasting algorithm criterion, determine the posterior sample of p value
	#Hint: use the formula in our lecture slide to compute the ratio of the likelihood at the proposal p value and the current p value
	#Then determine whether the proposal p value or the current p value should be accepted as the posterior sample
  alpha_accept = exp(yy.llik.prop - yy.llik.cur)
  if (runif(1) < alpha_accept) {
    p.cur = p.prop
    ind.cur = prop.ind
    yy.llik.cur = yy.llik.prop
  }
	p.post[gg]= p.cur
	
}

#plot the posterior and calculate the posterior mean
dev.new(width=5, height=5)
plot(p.post, type="l")
postmean = mean(p.post[1001:n.step])
postmean
