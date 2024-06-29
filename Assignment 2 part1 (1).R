#1
library(readr)

bankd <- read.csv("C:/Users/SAMSUNG/Downloads/Bank_Retention_Data.csv")

head(bankd)

# Convert TractID into a factor
bankd$TractID <- as.factor(bankd$TractID)

# Logistic Regression Model
logit_model1 <- glm(Churn ~ Age + Income + HomeVal + Tenure + DirectDeposit + Loan + Dist + MktShare, 
                    data = bankd, family = binomial(link = "logit"))

# Probit Regression Model
probit_model1 <- glm(Churn ~ Age + Income + HomeVal + Tenure + DirectDeposit + Loan + Dist + MktShare, 
                     data = bankd, family = binomial(link = "probit"))



aic_logit1 <- AIC(logit_model1)
aic_probit1 <- AIC(probit_model1)


summary(logit_model1)
summary(probit_model1)

aic_logit1
aic_probit1


#2
library(lme4)
library(Matrix)
# Fit the logistic regression model with a random intercept for TractID
glmer_model <- glmer(Churn ~ Age + Income + HomeVal + Tenure + DirectDeposit + Loan + Dist + MktShare +
                       (1 | TractID), 
                     data = bankd, 
                     family = binomial(link = "logit"),
                     glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=10000)))

aic_glmer_model <- AIC(glmer_model)
summary(glmer_model)
aic_glmer_model

#3
library(MCMCpack)
bayes_logit = MCMChlogit(fixed = Churn ~ Age + Income + HomeVal + Tenure + DirectDeposit + Loan + Dist + MktShare, 
                           random= ~1, 
                           group="TractID", 
                           data=bankd, 
                           mcmc=20000, 
                           r=2, 
                           R=1, 
                           burnin = 10000, 
                           thin=20,
                           FixOD = 1)


summary(bayes_logit$mcmc[,1:9])


# Plot the posterior sampling chains and densities for b2 and b5 
plot(bayes_logit$mcmc[, "beta.Income"])
plot(bayes_logit$mcmc[, "beta.DirectDeposit"])