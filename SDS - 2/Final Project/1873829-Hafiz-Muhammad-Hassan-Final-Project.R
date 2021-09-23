
# data taken form 
# https://raw.githubusercontent.com/GTPB/PSLS20/master/data/fev.txt
# we need to convert f 

require(tidyverse)
require(magrittr)
require(R2jags)
require(mcmcse)
require(bayesplot)
require(TeachingDemos)

# FEV data from Rosner; 
# 5th edition;   p. 479Use only the subset of data from individuals aged 10-19.   
# Age 10 was the youngest smokern=345y=FEV= Forced expiratory volume = a measure of lung function; 
# it ranges from about 2 to 5 in the data set;    higher is betterx=(Age, Smoke)Smoke is binary (yes/no)Ages range from 10 to 19 


FEVdata <- read.table(file="fev.txt",header=T, sep="\t")
head(ifelse(FEVdata$gender=="m", 1, 0))

FEVdata$gender = ifelse(FEVdata$gender=="m", 1, 0)

# plotting between age and fev
plot(x=FEVdata$age,y=FEVdata$fev,xlab = "Age",ylab = "FEV", col= c("darkgreen"))

# plotting between height and fev

plot(x=FEVdata$height,y=FEVdata$fev,xlab = "Height",ylab = "FEV", col= c("darkblue"))

# plotting between gender and fev

plot(x=FEVdata$gender,y=FEVdata$fev,xlab = "Male - Female",ylab = "FEV", col= c("darkblue"))

# plotting between gender and fev

plot(x=FEVdata$smoking,y=FEVdata$fev,xlab = " Non-Smoking - Smoking",ylab = "FEV", col= c("darkblue"))

#Now we plot the distribution of FEV
hist(FEVdata$fev,breaks = 50,xlab = "FEV Distribution")

#  Data Analysis and Frequentist inference

?coef
?lm

smoke<-FEVdata$smoking
fev<- FEVdata$fev
age<- FEVdata$age
height <- FEVdata$height
n <- length(fev)
lr <- lm(fev~age+smoke+height+age:smoke) # linear regression

# should display summary of the model 
summary(lr)


beta <- coef(lr)
X <- cbind(rep(1,n),age,smoke,age*smoke)

# building getting data from model itself
rstudsmoke <- rstudent(lr)[smoke==1]
rstudnonsmoke <- rstudent(lr)[smoke==0]


#We can observer that we get following estimates of the parameters. Now we analyize our linear line with the dataset from which the model is generated

?fitted
?resid
par(mfrow=c(2,2))
plot(fitted(lr),resid(lr),xlab="Fitted Values", ylab="Residuals",main="", )
lines(lowess(fitted(lr),resid(lr)),col= c("red"))
qqnorm(resid(lr),main="")
qqline(resid(lr), col=c("red"))


#Now we analyize our linear line with the dataset for non-smokers only


par(mfrow=c(2,2))
plot(fitted(lr)[smoke==0],rstudnonsmoke,xlab="Fitted Values", ylab="Studentized Residuals",main="Nonsmokers \n Residual Plot")
lines(lowess(fitted(lr)[smoke==0],rstudnonsmoke), col= c("red"))
qqnorm(rstudnonsmoke,main="Nonsmokers \n Normal Plot")
qqline(rstudnonsmoke,col= c("red"))


#Now we analyize our linear line with the dataset Age wise
plot(age, rstudent(lr),xlab="Age", ylab="Residuals",main="")
lines(lowess(age,rstudent(lr)), col= c("red"))


#We observe that FEV follows normal distribution so we try to develop its distribution
#Lets only include Age and Smoke features and analyize them. 



#Now we write the following model in Jags

model1.jags <- function() {
  # Likelihood
  for(i in 1:n){
    fev[i] ~ dnorm(mu[i],tau)
    mu[i] <- beta[1] + beta[2]*age[i] + beta[3]*smoke[i] + beta[4]*age[i]*smoke[i]
  }
  beta[1] ~ dnorm(0,0.001)  # Diffuse prior
  beta[2] ~ dnorm(0,0.001)
  beta[3] ~ dnorm(0,0.001)
  beta[4] ~ dnorm(0,0.001)
  tau ~ dgamma(0.001,0.001)
}

# <h4>The Parameter $\tau$ represents the variation in the data of FEV </h4>
#   <h4>The Parameter $\beta_1$ is linear transformation parameter used to balance out the mean of data</h4>
#   <h4>The Parameter $\beta_2$ gives the weight to Age variable, the larger this parameter means FEV depends more on value to Age </h4>
#   <h4>The Parameter $\beta_3$ gives the weight to Smoke variabee, the larger this parameter means FEV depends more on value to Smoke </h4>
#   <h4>The Parameter $\beta_4$ gives the weight to product of Smoke and Age variables, the larger this parameters means FEV depends more on product of Smoke and Age</h4>
# 

# Preparing data for JAGS
n <- length(FEVdata$age)

fev <- FEVdata$fev
smoke <- FEVdata$smoking
age <- FEVdata$age

#Lets only include Age and Smoke features and analyize them
dat.jags <- list("n","fev","smoke", "age")


# Defining parameters of interest
mod.params <- c("beta","tau") 

# Starting values
mod.inits <- function(){
  list("tau" = 1, "beta" = c(0,0,0,0))
}


# Run JAGS
set.seed(1873829)
mod1.fit <- jags(data = dat.jags,                                        # DATA
                model.file = model1.jags, inits = mod.inits,                  # MODEL
                parameters.to.save = mod.params,                  
                n.chains = 3, n.iter = 9000, n.burnin = 1000, n.thin=10)# MCMC
mod1.fit
chainArray <- mod.fit$BUGSoutput$sims.array



# Now we get the data of smokers and non-smokers from simulation and analyze it.
# Following is the jags model for simulating MCMC for smokers and non-smoke with Age

# https://oregonstate.app.box.com/s/b8niipwqjaviuujssb5hod6fnoevs1t4
# this model has been used in above page

model2.jags <- function()  {
  # Likelihood
  for(i in 1:n){
    fev[i] ~ dnorm(mu[i],tau)
    mu[i] <- beta[1] + beta[2]* age[i] + beta[3]*smoke[i]+ beta[4]*age[i]*smoke[i]
  }
  
  beta[1] ~ dnorm(0, 0.001)  # Diffuse prior
  beta[2] ~ dnorm(0, 0.001)
  beta[3] ~ dnorm(0, 0.001)
  beta[4] ~ dnorm(0, 0.001)
  tau ~ dgamma(0.001, 0.001)
                  
  ## Estimate mean FEV for smokers and nonsmokers who are 10, ..., 19 years old
  for(i in 1:10){
    meanFEVs[i] <-  beta[1] + (beta[2]+beta[4])*(i+9) + beta[3]
    meanFEVns[i] <-  beta[1] + beta[2]*(i+9)
  }
  
  meanFEVnsZeroError <- ifelse(meanFEVns[9] != 0, meanFEVns[9], 1)
  RM <- meanFEVns[9]/meanFEVnsZeroError ## RM comparing 18 year old smoker to 18 year old nonsmoker
  MD <- meanFEVs[9] - meanFEVns[5] ## MD comparing 18 year old smoker to 13 year old nonsmoker
  ## Easy to estimate relative means and mean differences as well
  
  
  ## Predict the FEV for a 20 year old smoker and nonsmoker
  FEV20s ~ dnorm(mu20s, tau)
  FEV20ns ~ dnorm(mu20ns, tau)
  mu20s <-  beta[1] + (beta[2]+beta[4])*20 + beta[3]
  mu20ns <-  beta[1] + beta[2]*20 
}


fev<- FEVdata$fev
age<- FEVdata$age
smoke<- FEVdata$smoking
height<- FEVdata$height
n <- length(FEVdata$age)

# data that jags will use
dat.jags <- list("n","fev","age","smoke")

# parameters of intrests
mod.params  <- c("beta","tau","meanFEVs","meanFEVns","RM","MD","FEV20s","FEV20ns") 

# Starting values
mod.inits <- function(){
  list("tau" = 1, "beta" = c(0,0,0,0))
}

set.seed(1873829)

mod2.fit <- jags(data = dat.jags,                                        # DATA
                model.file = model2.jags, inits = mod.inits,                  # MODEL
                parameters.to.save = mod.params,                  
                n.chains = 3, n.iter = 9000, n.burnin = 1000, n.thin=10) # MCMC
mod2.fit

# Now lets plot the data of smokers and nonsmokers with their age.

BayesSmoker <- mod.fit$BUGSoutput$mean$meanFEVs
BayesNonsmoker <- mod.fit$BUGSoutput$mean$meanFEVns

Age1 <- seq(10,19,1)
plot(age,fev,xlab="Age", ylab="FEV", lwd=1.5)
lines(Age1,BayesSmoker,type='l', lty=1, lwd=2,cex=1.2)
lines(Age1,BayesNonsmoker,type='l', lty=2, lwd=2,cex=1.2)
legend("bottomright", c("Smoker", "Nonsmoker"), lty=c(1,2), lwd=2)


#<h4>Now we observe that with the age smokers have let FEV compared to non-smokers. </h4>
# <h2>Comparison</h2>
# <h4>Lets compare baysean with frequentist</h4>

summary(lr)
mod1.fit$BUGSoutput$summary


# Diagnostic
# plot(mod1.fit)
# traceplot(mod1.fit)

# We can get better diagnostics with other packages

# Plots with BayesPlot
chainArray <- mod1.fit$BUGSoutput$sims.array
bayesplot::mcmc_combo(chainArray) # combo of density plot 
bayesplot::mcmc_acf(chainArray)


# Diagnostic with coda
coda.fit <- coda::as.mcmc(mod1.fit)

coda::acfplot(coda.fit)
coda::densplot(coda.fit)

#diagnostic
coda::raftery.diag(coda.fit) # checking how many iteration we need. 

coda::geweke.diag(coda.fit) # 
coda::geweke.plot(coda.fit)


coda::gelman.diag(coda.fit)
coda::gelman.plot(coda.fit)

coda::heidel.diag(coda.fit)



# Joining chains --------------------------------------------------------------

chainMat <- mod1.fit$BUGSoutput$sims.matrix






# Intervals
cred <- 0.95
# Point estimates
(p.hat.jags <- colMeans(chainMat))

# Intervals
(p.ET.jags <- apply(chainMat, 2, quantile, 
                    prob=c((1-cred)/2, 1-(1-cred)/2)))

# What about the HPD?
(p.HPD.jags <- coda::HPDinterval(as.mcmc(chainMat)))




p.hat.jags[2]
p.ET.jags
p.HPD.jags
