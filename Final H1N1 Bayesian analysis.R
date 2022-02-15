#Final H1N1 Analysis

library(R2jags)
library(coda)
library(lmtest)
library(foreign)
library(plyr)
library(dplyr)
library(varhandle)

setwd("/Users/sewraith/Box/Research/Homotypic/Data/tempwork/Repeat")

### Analyses of interest: H1N1 ----

# Bring in the data
data_orig <- read.csv2("faustodata_updated.csv", header = T, sep = ",") 


# H1N1 2011 to 2013
try1 <- data_orig %>%
  filter(at2011==1 & at2013==1) %>%
  dplyr::select(Obs, H1N1pdm_2013, H1N1pdm_2011)
try1$H1N1pdm_2013 <- as.numeric(try1$H1N1pdm_2013)
try1$H1N1pdm_2011 <- as.numeric(try1$H1N1pdm_2011)

table(try1$H1N1pdm_2013, try1$H1N1pdm_2011)

# Putting exposure and outcome variables in the environment 
H1N1pdm_2013 <- try1$H1N1pdm_2013
H1N1pdm_2011 <- try1$H1N1pdm_2011

### Bayesian logistic regression
logistic.model <- function() {
  # SAMPLING DISTRIBUTION
  for (i in 1:N) {
    logit(pi[i]) <- b[1] + b[2]*H1N1pdm_2011[i];
    H1N1pdm_2013[i] ~ dbin(pi[i],1);  
  }
  
  # PRIORS ON BETAS
  b[1:Nx] ~ dmnorm(mu[1:Nx], tau[1:Nx,1:Nx]) # multivariate normal prior
  
  # Calculate ORs:
  for (l in 1:Nx) {
    OR[l] <- exp(b[l]);
  }
}

N <- length(H1N1pdm_2013) # Number of observations
Nx <- 2 # Number of variables: intercept + exposure 

# Data, parameter list and starting values
mu <- rep(0, Nx)
tau <- diag(0.001, Nx)
data.logistic <- list("N", "Nx", "H1N1pdm_2011", "H1N1pdm_2013", "mu", "tau")
parameters.logistic <-c("b","OR") # Parameters to keep track of
inits.logistic <- function() {list (b= rep(0, Nx))}

set.seed(114011)
logistic.sim<-jags(data=data.logistic, inits=inits.logistic, parameters.logistic, n.iter=100000,
                   n.burn=25000, model.file=logistic.model, n.thin=5, n.chains = 3)

print(logistic.sim, 2)
OR_H1N1_2011_2013_Bayes <- round(logistic.sim$BUGSoutput$summary[2,c(5,3,7)], 2)
OR_H1N1_2011_2013_Bayes


# H1N1 2011 to 2015
try50 <- data_orig %>%
  filter(at2011==1 & at2015==1 & H1N1pdm_2013==0 ) %>%
  dplyr::select(Obs, H1N1pdm_2011, H1N1pdm_2015)
try50$H1N1pdm_2015 <- as.numeric(try50$H1N1pdm_2015)
try50$H1N1pdm_2011 <- as.numeric(try50$H1N1pdm_2011)

table(try50$H1N1pdm_2015, try50$H1N1pdm_2011)

# Bayesian analyses
# Putting exposure and outcome variables in the environment 
H1N1pdm_2015 <- try50$H1N1pdm_2015
H1N1pdm_2011 <- try50$H1N1pdm_2011

### Bayesian logistic regression
logistic.model <- function() {
  # SAMPLING DISTRIBUTION
  for (i in 1:N) {
    logit(pi[i]) <- b[1] + b[2]*H1N1pdm_2011[i];
    H1N1pdm_2015[i] ~ dbin(pi[i],1);  
  }
  
  # PRIORS ON BETAS
  b[1:Nx] ~ dmnorm(mu[1:Nx], tau[1:Nx,1:Nx]) # multivariate normal prior
  
  # Calculate ORs:
  for (l in 1:Nx) {
    OR[l] <- exp(b[l]);
  }
}

N <- length(H1N1pdm_2015) # Number of observations
Nx <- 2 # Number of variables: intercept + exposure 

# Data, parameter list and starting values
mu <- rep(0, Nx)
tau <- diag(0.001, Nx)
data.logistic <- list("N", "Nx", "H1N1pdm_2011", "H1N1pdm_2015", "mu", "tau")
parameters.logistic <-c("b","OR") # Parameters to keep track of
inits.logistic <- function() {list (b= rep(0, Nx))}

set.seed(114011)
logistic.sim<-jags(data=data.logistic, inits=inits.logistic, parameters.logistic, n.iter=100000,
                   n.burn=25000, model.file=logistic.model, n.thin=5, n.chains = 3)

print(logistic.sim, 2)
OR_H1N1_2011_2015_Bayes_SA <- round(logistic.sim$BUGSoutput$summary[2,c(5,3,7)], 2)
OR_H1N1_2011_2015_Bayes_SA



# H1N1 2011 to 2018
try50 <- data_orig %>%
  filter(at2011==1 & at2018==1 & H1N1pdm_2013==0 & H1N1pdm_2015==0 & H1N1pdm_2016==0 & H1N1pdm_2017==0) %>%
  dplyr::select(Obs, H1N1pdm_2011, H1N1pdm_2018)
try50$H1N1pdm_2018 <- as.numeric(try50$H1N1pdm_2018)
try50$H1N1pdm_2011 <- as.numeric(try50$H1N1pdm_2011)

table(try50$H1N1pdm_2018, try50$H1N1pdm_2011)

# Bayesian analyses
# Putting exposure and outcome variables in the environment 
H1N1pdm_2018 <- try50$H1N1pdm_2018
H1N1pdm_2011 <- try50$H1N1pdm_2011

### Bayesian logistic regression
logistic.model <- function() {
  # SAMPLING DISTRIBUTION
  for (i in 1:N) {
    logit(pi[i]) <- b[1] + b[2]*H1N1pdm_2011[i];
    H1N1pdm_2018[i] ~ dbin(pi[i],1);  
  }
  
  # PRIORS ON BETAS
  b[1:Nx] ~ dmnorm(mu[1:Nx], tau[1:Nx,1:Nx]) # multivariate normal prior
  
  # Calculate ORs:
  for (l in 1:Nx) {
    OR[l] <- exp(b[l]);
  }
}

N <- length(H1N1pdm_2018) # Number of observations
Nx <- 2 # Number of variables: intercept + exposure 

# Data, parameter list and starting values
mu <- rep(0, Nx)
tau <- diag(0.001, Nx)
data.logistic <- list("N", "Nx", "H1N1pdm_2011", "H1N1pdm_2018", "mu", "tau")
parameters.logistic <-c("b","OR") # Parameters to keep track of
inits.logistic <- function() {list (b= rep(0, Nx))}

set.seed(114011)
logistic.sim<-jags(data=data.logistic, inits=inits.logistic, parameters.logistic, n.iter=100000,
                   n.burn=25000, model.file=logistic.model, n.thin=5, n.chains = 3)

print(logistic.sim, 2)
OR_H1N1_2011_2018_Bayes_SA <- round(logistic.sim$BUGSoutput$summary[2,c(5,3,7)], 2)
OR_H1N1_2011_2018_Bayes_SA


# H1N1 2013 to 2015
try2 <- data_orig %>%
  filter(at2013==1 & at2015==1) %>%
  dplyr::select(Obs, H1N1pdm_2013, H1N1pdm_2015)
try2$H1N1pdm_2013 <- as.numeric(try2$H1N1pdm_2013)
try2$H1N1pdm_2015 <- as.numeric(try2$H1N1pdm_2015)

table(try2$H1N1pdm_2013, try2$H1N1pdm_2015)

# Bayesian analyses
# Putting exposure and outcome variables in the environment 
H1N1pdm_2013 <- try2$H1N1pdm_2013
H1N1pdm_2015 <- try2$H1N1pdm_2015

### Bayesian logistic regression
logistic.model <- function() {
  # SAMPLING DISTRIBUTION
  for (i in 1:N) {
    logit(pi[i]) <- b[1] + b[2]*H1N1pdm_2013[i];
    H1N1pdm_2015[i] ~ dbin(pi[i],1);  
  }
  
  # PRIORS ON BETAS
  b[1:Nx] ~ dmnorm(mu[1:Nx], tau[1:Nx,1:Nx]) # multivariate normal prior
  
  # Calculate ORs:
  for (l in 1:Nx) {
    OR[l] <- exp(b[l]);
  }
}

N <- length(H1N1pdm_2015) # Number of observations
Nx <- 2 # Number of variables: intercept + exposure 

# Data, parameter list and starting values
mu <- rep(0, Nx)
tau <- diag(0.001, Nx)
data.logistic <- list("N", "Nx", "H1N1pdm_2013", "H1N1pdm_2015", "mu", "tau")
parameters.logistic <-c("b","OR") # Parameters to keep track of
inits.logistic <- function() {list (b= rep(0, Nx))}

set.seed(114011)
logistic.sim<-jags(data=data.logistic, inits=inits.logistic, parameters.logistic, n.iter=100000,
                   n.burn=25000, model.file=logistic.model, n.thin=5, n.chains = 3)

print(logistic.sim, 2)
OR_H1N1_2013_2015_Bayes <- round(logistic.sim$BUGSoutput$summary[2,c(5,3,7)], 2)


# H1N1 2013 to 2018
try50 <- data_orig %>%
  filter(at2013==1 & at2018==1 & H1N1pdm_2015==0 & H1N1pdm_2016==0 & H1N1pdm_2017==0) %>%
  dplyr::select(Obs, H1N1pdm_2013, H1N1pdm_2018)
try50$H1N1pdm_2018 <- as.numeric(try50$H1N1pdm_2018)
try50$H1N1pdm_2013 <- as.numeric(try50$H1N1pdm_2013)

table(try50$H1N1pdm_2018, try50$H1N1pdm_2013)

# Bayesian analyses
# Putting exposure and outcome variables in the environment 
H1N1pdm_2018 <- try50$H1N1pdm_2018
H1N1pdm_2013 <- try50$H1N1pdm_2013

### Bayesian logistic regression
logistic.model <- function() {
  # SAMPLING DISTRIBUTION
  for (i in 1:N) {
    logit(pi[i]) <- b[1] + b[2]*H1N1pdm_2013[i];
    H1N1pdm_2018[i] ~ dbin(pi[i],1);  
  }
  
  # PRIORS ON BETAS
  b[1:Nx] ~ dmnorm(mu[1:Nx], tau[1:Nx,1:Nx]) # multivariate normal prior
  
  # Calculate ORs:
  for (l in 1:Nx) {
    OR[l] <- exp(b[l]);
  }
}

N <- length(H1N1pdm_2018) # Number of observations
Nx <- 2 # Number of variables: intercept + exposure 

# Data, parameter list and starting values
mu <- rep(0, Nx)
tau <- diag(0.001, Nx)
data.logistic <- list("N", "Nx", "H1N1pdm_2013", "H1N1pdm_2018", "mu", "tau")
parameters.logistic <-c("b","OR") # Parameters to keep track of
inits.logistic <- function() {list (b= rep(0, Nx))}

set.seed(114011)
logistic.sim<-jags(data=data.logistic, inits=inits.logistic, parameters.logistic, n.iter=100000,
                   n.burn=25000, model.file=logistic.model, n.thin=5, n.chains = 3)

print(logistic.sim, 2)
OR_H1N1_2013_2018_Bayes_SA <- round(logistic.sim$BUGSoutput$summary[2,c(5,3,7)], 2)
OR_H1N1_2013_2018_Bayes_SA



# H1N1 2015 to 2018
try50 <- data_orig %>%
  filter(at2015==1 & at2018==1 & H1N1pdm_2016==0 & H1N1pdm_2017==0) %>%
  dplyr::select(Obs, H1N1pdm_2015, H1N1pdm_2018)
try50$H1N1pdm_2018 <- as.numeric(try50$H1N1pdm_2018)
try50$H1N1pdm_2015 <- as.numeric(try50$H1N1pdm_2015)

table(try50$H1N1pdm_2018, try50$H1N1pdm_2015)

# Bayesian analyses
# Putting exposure and outcome variables in the environment 
H1N1pdm_2018 <- try50$H1N1pdm_2018
H1N1pdm_2015 <- try50$H1N1pdm_2015

### Bayesian logistic regression
logistic.model <- function() {
  # SAMPLING DISTRIBUTION
  for (i in 1:N) {
    logit(pi[i]) <- b[1] + b[2]*H1N1pdm_2015[i];
    H1N1pdm_2018[i] ~ dbin(pi[i],1);  
  }
  
  # PRIORS ON BETAS
  b[1:Nx] ~ dmnorm(mu[1:Nx], tau[1:Nx,1:Nx]) # multivariate normal prior
  
  # Calculate ORs:
  for (l in 1:Nx) {
    OR[l] <- exp(b[l]);
  }
}

N <- length(H1N1pdm_2018) # Number of observations
Nx <- 2 # Number of variables: intercept + exposure 

# Data, parameter list and starting values
mu <- rep(0, Nx)
tau <- diag(0.001, Nx)
data.logistic <- list("N", "Nx", "H1N1pdm_2015", "H1N1pdm_2018", "mu", "tau")
parameters.logistic <-c("b","OR") # Parameters to keep track of
inits.logistic <- function() {list (b= rep(0, Nx))}

set.seed(114011)
logistic.sim<-jags(data=data.logistic, inits=inits.logistic, parameters.logistic, n.iter=100000,
                   n.burn=25000, model.file=logistic.model, n.thin=5, n.chains = 3)

print(logistic.sim, 2)
OR_H1N1_2015_2018_Bayes_SA <- round(logistic.sim$BUGSoutput$summary[2,c(5,3,7)], 2)
OR_H1N1_2015_2018_Bayes_SA

#---- Analyses of interest: without dropped infections ----

# H1N1 2011 to 2013 - nothing was excluded

# H1N1 2011 to 2015
try50 <- data_orig %>%
  filter(at2011==1 & at2015==1) %>%
  dplyr::select(Obs, H1N1pdm_2011, H1N1pdm_2015)
try50$H1N1pdm_2015 <- as.numeric(try50$H1N1pdm_2015)
try50$H1N1pdm_2011 <- as.numeric(try50$H1N1pdm_2011)

table(try50$H1N1pdm_2015, try50$H1N1pdm_2011)

# Bayesian analyses
# Putting exposure and outcome variables in the environment 
H1N1pdm_2015 <- try50$H1N1pdm_2015
H1N1pdm_2011 <- try50$H1N1pdm_2011

### Bayesian logistic regression
logistic.model <- function() {
  # SAMPLING DISTRIBUTION
  for (i in 1:N) {
    logit(pi[i]) <- b[1] + b[2]*H1N1pdm_2011[i];
    H1N1pdm_2015[i] ~ dbin(pi[i],1);  
  }
  
  # PRIORS ON BETAS
  b[1:Nx] ~ dmnorm(mu[1:Nx], tau[1:Nx,1:Nx]) # multivariate normal prior
  
  # Calculate ORs:
  for (l in 1:Nx) {
    OR[l] <- exp(b[l]);
  }
}

N <- length(H1N1pdm_2015) # Number of observations
Nx <- 2 # Number of variables: intercept + exposure 

# Data, parameter list and starting values
mu <- rep(0, Nx)
tau <- diag(0.001, Nx)
data.logistic <- list("N", "Nx", "H1N1pdm_2011", "H1N1pdm_2015", "mu", "tau")
parameters.logistic <-c("b","OR") # Parameters to keep track of
inits.logistic <- function() {list (b= rep(0, Nx))}

set.seed(114011)
logistic.sim<-jags(data=data.logistic, inits=inits.logistic, parameters.logistic, n.iter=100000,
                   n.burn=25000, model.file=logistic.model, n.thin=5, n.chains = 3)

print(logistic.sim, 2)
OR_H1N1_2011_2015_Bayes_SA <- round(logistic.sim$BUGSoutput$summary[2,c(5,3,7)], 2)
OR_H1N1_2011_2015_Bayes_SA


# H1N1 2011 to 2018
try50 <- data_orig %>%
  filter(at2011==1 & at2018==1 ) %>%
  dplyr::select(Obs, H1N1pdm_2011, H1N1pdm_2018)
try50$H1N1pdm_2018 <- as.numeric(try50$H1N1pdm_2018)
try50$H1N1pdm_2011 <- as.numeric(try50$H1N1pdm_2011)

table(try50$H1N1pdm_2018, try50$H1N1pdm_2011)

# Bayesian analyses
# Putting exposure and outcome variables in the environment 
H1N1pdm_2018 <- try50$H1N1pdm_2018
H1N1pdm_2011 <- try50$H1N1pdm_2011

### Bayesian logistic regression
logistic.model <- function() {
  # SAMPLING DISTRIBUTION
  for (i in 1:N) {
    logit(pi[i]) <- b[1] + b[2]*H1N1pdm_2011[i];
    H1N1pdm_2018[i] ~ dbin(pi[i],1);  
  }
  
  # PRIORS ON BETAS
  b[1:Nx] ~ dmnorm(mu[1:Nx], tau[1:Nx,1:Nx]) # multivariate normal prior
  
  # Calculate ORs:
  for (l in 1:Nx) {
    OR[l] <- exp(b[l]);
  }
}

N <- length(H1N1pdm_2018) # Number of observations
Nx <- 2 # Number of variables: intercept + exposure 

# Data, parameter list and starting values
mu <- rep(0, Nx)
tau <- diag(0.001, Nx)
data.logistic <- list("N", "Nx", "H1N1pdm_2011", "H1N1pdm_2018", "mu", "tau")
parameters.logistic <-c("b","OR") # Parameters to keep track of
inits.logistic <- function() {list (b= rep(0, Nx))}

set.seed(114011)
logistic.sim<-jags(data=data.logistic, inits=inits.logistic, parameters.logistic, n.iter=100000,
                   n.burn=25000, model.file=logistic.model, n.thin=5, n.chains = 3)

print(logistic.sim, 2)
OR_H1N1_2011_2018_Bayes_SA <- round(logistic.sim$BUGSoutput$summary[2,c(5,3,7)], 2)
OR_H1N1_2011_2018_Bayes_SA

#2013 to 2015 - nothing was excluded

# H1N1 2013 to 2018
try50 <- data_orig %>%
  filter(at2013==1 & at2018==1 ) %>%
  dplyr::select(Obs, H1N1pdm_2013, H1N1pdm_2018)
try50$H1N1pdm_2018 <- as.numeric(try50$H1N1pdm_2018)
try50$H1N1pdm_2013 <- as.numeric(try50$H1N1pdm_2013)

table(try50$H1N1pdm_2018, try50$H1N1pdm_2013)

# Bayesian analyses
# Putting exposure and outcome variables in the environment 
H1N1pdm_2018 <- try50$H1N1pdm_2018
H1N1pdm_2013 <- try50$H1N1pdm_2013

### Bayesian logistic regression
logistic.model <- function() {
  # SAMPLING DISTRIBUTION
  for (i in 1:N) {
    logit(pi[i]) <- b[1] + b[2]*H1N1pdm_2013[i];
    H1N1pdm_2018[i] ~ dbin(pi[i],1);  
  }
  
  # PRIORS ON BETAS
  b[1:Nx] ~ dmnorm(mu[1:Nx], tau[1:Nx,1:Nx]) # multivariate normal prior
  
  # Calculate ORs:
  for (l in 1:Nx) {
    OR[l] <- exp(b[l]);
  }
}

N <- length(H1N1pdm_2018) # Number of observations
Nx <- 2 # Number of variables: intercept + exposure 

# Data, parameter list and starting values
mu <- rep(0, Nx)
tau <- diag(0.001, Nx)
data.logistic <- list("N", "Nx", "H1N1pdm_2013", "H1N1pdm_2018", "mu", "tau")
parameters.logistic <-c("b","OR") # Parameters to keep track of
inits.logistic <- function() {list (b= rep(0, Nx))}

set.seed(114011)
logistic.sim<-jags(data=data.logistic, inits=inits.logistic, parameters.logistic, n.iter=100000,
                   n.burn=25000, model.file=logistic.model, n.thin=5, n.chains = 3)

print(logistic.sim, 2)
OR_H1N1_2013_2018_Bayes_SA <- round(logistic.sim$BUGSoutput$summary[2,c(5,3,7)], 2)
OR_H1N1_2013_2018_Bayes_SA



# H1N1 2015 to 2018
try50 <- data_orig %>%
  filter(at2015==1 & at2018==1 ) %>%
  dplyr::select(Obs, H1N1pdm_2015, H1N1pdm_2018)
try50$H1N1pdm_2018 <- as.numeric(try50$H1N1pdm_2018)
try50$H1N1pdm_2015 <- as.numeric(try50$H1N1pdm_2015)

table(try50$H1N1pdm_2018, try50$H1N1pdm_2015)

# Bayesian analyses
# Putting exposure and outcome variables in the environment 
H1N1pdm_2018 <- try50$H1N1pdm_2018
H1N1pdm_2015 <- try50$H1N1pdm_2015

### Bayesian logistic regression
logistic.model <- function() {
  # SAMPLING DISTRIBUTION
  for (i in 1:N) {
    logit(pi[i]) <- b[1] + b[2]*H1N1pdm_2015[i];
    H1N1pdm_2018[i] ~ dbin(pi[i],1);  
  }
  
  # PRIORS ON BETAS
  b[1:Nx] ~ dmnorm(mu[1:Nx], tau[1:Nx,1:Nx]) # multivariate normal prior
  
  # Calculate ORs:
  for (l in 1:Nx) {
    OR[l] <- exp(b[l]);
  }
}

N <- length(H1N1pdm_2018) # Number of observations
Nx <- 2 # Number of variables: intercept + exposure 

# Data, parameter list and starting values
mu <- rep(0, Nx)
tau <- diag(0.001, Nx)
data.logistic <- list("N", "Nx", "H1N1pdm_2015", "H1N1pdm_2018", "mu", "tau")
parameters.logistic <-c("b","OR") # Parameters to keep track of
inits.logistic <- function() {list (b= rep(0, Nx))}

set.seed(114011)
logistic.sim<-jags(data=data.logistic, inits=inits.logistic, parameters.logistic, n.iter=100000,
                   n.burn=25000, model.file=logistic.model, n.thin=5, n.chains = 3)

print(logistic.sim, 2)
OR_H1N1_2015_2018_Bayes_SA <- round(logistic.sim$BUGSoutput$summary[2,c(5,3,7)], 2)
OR_H1N1_2015_2018_Bayes_SA





#---- Analyses of interest: Age-stratified  ----

#2011 to 2013

try7 <- data_orig %>%
  filter(at2011==1 & at2013==1) %>%
  dplyr::select(Obs, H1N1pdm_2011, H1N1pdm_2013, age_2011_cat)

try2$H1N1pdm_2011 <- unfactor(try2$H1N1pdm_2011)
try2$H1N1pdm_2013 <- unfactor(try2$H1N1pdm_2013)

# 2x2 table
table(try7[try7$age_2011_cat==0,]$H1N1pdm_2011, try7[try7$age_2011_cat==0,]$H1N1pdm_2013) 
table(try7[try7$age_2011_cat==1,]$H1N1pdm_2011, try7[try7$age_2011_cat==1,]$H1N1pdm_2013) 

# Freq analysis of H1N1 2011 on H1N1 2013, under 5s
freq.fit <- glm(H1N1pdm_2013 ~ H1N1pdm_2011, family=binomial, data=try7[try7$age_2011_cat==0,])
summary(freq.fit) 
OR.freq <- exp(cbind(coef(freq.fit), confint(freq.fit)))
colnames(OR.freq) <- c("OR","95% LL","95% UL")
round(OR.freq, digits=2) 
OR_H1N1_2011_2013_cat0_freq <- round(OR.freq, digits=2)[2,]

# Freq analysis of H1N1 2011 on H1N1 2013, over 5s
freq.fit <- glm(H1N1pdm_2013 ~ H1N1pdm_2011, family=binomial, data=try7[try7$age_2011_cat==1,])
summary(freq.fit) 
OR.freq <- exp(cbind(coef(freq.fit), confint(freq.fit)))
colnames(OR.freq) <- c("OR","95% LL","95% UL")
round(OR.freq, digits=2) 
OR_H1N1_2011_2013_cat1_freq <- round(OR.freq, digits=2)[2,]

# Bayesian analysis of H1N1 2011 on H1N1 2013, under 5s
H1N1pdm_2011 <- try7[try7$age_2011_cat==0,]$H1N1pdm_2011
H1N1pdm_2013 <- try7[try7$age_2011_cat==0,]$H1N1pdm_2013

logistic.model <- function() {
  # SAMPLING DISTRIBUTION
  for (i in 1:N) {
    logit(pi[i]) <- b[1] + b[2]*H1N1pdm_2011[i];
    H1N1pdm_2013[i] ~ dbin(pi[i],1);  
  }
  
  # PRIORS ON BETAS
  b[1:Nx] ~ dmnorm(mu[1:Nx], tau[1:Nx,1:Nx]) # multivariate normal prior
  
  # Calculate ORs:
  for (l in 1:Nx) {
    OR[l] <- exp(b[l]);
  }
}

N <- length(H1N1pdm_2013) 
Nx <- 2 

mu <- rep(0, Nx)
tau <- diag(0.001, Nx)
data.logistic <- list("N", "Nx", "H1N1pdm_2011", "H1N1pdm_2013", "mu", "tau")
parameters.logistic <-c("b","OR") # Parameters to keep track of
inits.logistic <- function() {list (b= rep(0, Nx))}

set.seed(114011)
logistic.sim_2011_2013_cat0<-jags(data=data.logistic, inits=inits.logistic, parameters.logistic, n.iter=100000,
                                  n.burn=25000, model.file=logistic.model, n.thin=5, n.chains = 3)

print(logistic.sim_2011_2013_cat0, 2)
OR_H1N1_2011_2013_cat0_Bayes <- round(logistic.sim_2011_2013_cat0$BUGSoutput$summary[2,c(5,3,7)], 2)

OR_H1N1_2011_2013_cat0_freq
OR_H1N1_2011_2013_cat0_Bayes

# Bayesian analysis of H1N1 2011 on H1N1 2013, over 5s
H1N1pdm_2011 <- try7[try7$age_2011_cat==1,]$H1N1pdm_2011
H1N1pdm_2013 <- try7[try7$age_2011_cat==1,]$H1N1pdm_2013

logistic.model <- function() {
  # SAMPLING DISTRIBUTION
  for (i in 1:N) {
    logit(pi[i]) <- b[1] + b[2]*H1N1pdm_2011[i];
    H1N1pdm_2013[i] ~ dbin(pi[i],1);  
  }
  
  # PRIORS ON BETAS
  b[1:Nx] ~ dmnorm(mu[1:Nx], tau[1:Nx,1:Nx]) # multivariate normal prior
  
  # Calculate ORs:
  for (l in 1:Nx) {
    OR[l] <- exp(b[l]);
  }
}

N <- length(H1N1pdm_2013) 
Nx <- 2 

mu <- rep(0, Nx)
tau <- diag(0.001, Nx)
data.logistic <- list("N", "Nx", "H1N1pdm_2011", "H1N1pdm_2013", "mu", "tau")
parameters.logistic <-c("b","OR") # Parameters to keep track of
inits.logistic <- function() {list (b= rep(0, Nx))}

set.seed(114011)
logistic.sim_2011_2013_cat1<-jags(data=data.logistic, inits=inits.logistic, parameters.logistic, n.iter=100000,
                                  n.burn=25000, model.file=logistic.model, n.thin=5, n.chains = 3)

print(logistic.sim_2011_2013_cat1, 2)
OR_H1N1_2011_2013_cat1_Bayes <- round(logistic.sim_2011_2013_cat1$BUGSoutput$summary[2,c(5,3,7)], 2)

OR_H1N1_2011_2013_cat1_freq
OR_H1N1_2011_2013_cat1_Bayes


### Analyses of interest: Age-stratified H1N1 2013 to 2015 ----
try8 <- data_orig %>%
  filter(at2013==1 & at2015==1) %>%
  dplyr::select(Obs, H1N1pdm_2013, H1N1pdm_2015, age_2013_cat)

try8$H1N1pdm_2013 <- unfactor(try8$H1N1pdm_2013)
try8$H1N1pdm_2015 <- unfactor(try8$H1N1pdm_2015)

# 2x2 table
table(try8[try8$age_2013_cat==0,]$H1N1pdm_2013, try8[try8$age_2013_cat==0,]$H1N1pdm_2015) 
table(try8[try8$age_2013_cat==1,]$H1N1pdm_2013, try8[try8$age_2013_cat==1,]$H1N1pdm_2015) 

# Freq analysis of H1N1 2013 on H1N1 2015, under 5s
freq.fit <- glm(H1N1pdm_2015 ~ H1N1pdm_2013, family=binomial, data=try8[try8$age_2013_cat==0,])
summary(freq.fit) 
OR.freq <- exp(cbind(coef(freq.fit), confint(freq.fit)))
colnames(OR.freq) <- c("OR","95% LL","95% UL")
round(OR.freq, digits=2) 
OR_H1N1_2013_2015_cat0_freq <- round(OR.freq, digits=2)[2,]

# Freq analysis of H1N1 2013 on H1N1 2015, over 5s
freq.fit <- glm(H1N1pdm_2015 ~ H1N1pdm_2013, family=binomial, data=try8[try8$age_2013_cat==1,])
summary(freq.fit) 
OR.freq <- exp(cbind(coef(freq.fit), confint(freq.fit)))
colnames(OR.freq) <- c("OR","95% LL","95% UL")
round(OR.freq, digits=2) 
OR_H1N1_2013_2015_cat1_freq <- round(OR.freq, digits=2)[2,]

# Bayesian analysis of H1N1 2013 on H1N1 2015, under 5s
H1N1pdm_2013 <- try8[try8$age_2013_cat==0,]$H1N1pdm_2013
H1N1pdm_2015 <- try8[try8$age_2013_cat==0,]$H1N1pdm_2015

logistic.model <- function() {
  # SAMPLING DISTRIBUTION
  for (i in 1:N) {
    logit(pi[i]) <- b[1] + b[2]*H1N1pdm_2013[i];
    H1N1pdm_2015[i] ~ dbin(pi[i],1);  
  }
  
  # PRIORS ON BETAS
  b[1:Nx] ~ dmnorm(mu[1:Nx], tau[1:Nx,1:Nx]) # multivariate normal prior
  
  # Calculate ORs:
  for (l in 1:Nx) {
    OR[l] <- exp(b[l]);
  }
}

N <- length(H1N1pdm_2015) 
Nx <- 2 

mu <- rep(0, Nx)
tau <- diag(0.001, Nx)
data.logistic <- list("N", "Nx", "H1N1pdm_2013", "H1N1pdm_2015", "mu", "tau")
parameters.logistic <-c("b","OR") # Parameters to keep track of
inits.logistic <- function() {list (b= rep(0, Nx))}

set.seed(114011)
logistic.sim_2013_2015_cat0<-jags(data=data.logistic, inits=inits.logistic, parameters.logistic, n.iter=100000,
                                  n.burn=25000, model.file=logistic.model, n.thin=5, n.chains = 3)

print(logistic.sim_2013_2015_cat0, 2)
OR_H1N1_2013_2015_cat0_Bayes <- round(logistic.sim_2013_2015_cat0$BUGSoutput$summary[2,c(5,3,7)], 2)

OR_H1N1_2013_2015_cat0_freq
OR_H1N1_2013_2015_cat0_Bayes


# Bayesian analysis of H1N1 2013 on H1N1 2015, over 5s
H1N1pdm_2013 <- try8[try8$age_2013_cat==1,]$H1N1pdm_2013
H1N1pdm_2015 <- try8[try8$age_2013_cat==1,]$H1N1pdm_2015

logistic.model <- function() {
  # SAMPLING DISTRIBUTION
  for (i in 1:N) {
    logit(pi[i]) <- b[1] + b[2]*H1N1pdm_2013[i];
    H1N1pdm_2015[i] ~ dbin(pi[i],1);  
  }
  
  # PRIORS ON BETAS
  b[1:Nx] ~ dmnorm(mu[1:Nx], tau[1:Nx,1:Nx]) # multivariate normal prior
  
  # Calculate ORs:
  for (l in 1:Nx) {
    OR[l] <- exp(b[l]);
  }
}

N <- length(H1N1pdm_2015) 
Nx <- 2 

mu <- rep(0, Nx)
tau <- diag(0.001, Nx)
data.logistic <- list("N", "Nx", "H1N1pdm_2013", "H1N1pdm_2015", "mu", "tau")
parameters.logistic <-c("b","OR") # Parameters to keep track of
inits.logistic <- function() {list (b= rep(0, Nx))}

set.seed(114011)
logistic.sim_2013_2015_cat1<-jags(data=data.logistic, inits=inits.logistic, parameters.logistic, n.iter=100000,
                                  n.burn=25000, model.file=logistic.model, n.thin=5, n.chains = 3)

print(logistic.sim_2013_2015_cat1, 2)
OR_H1N1_2013_2015_cat1_Bayes <- round(logistic.sim_2013_2015_cat1$BUGSoutput$summary[2,c(5,3,7)], 2)

OR_H1N1_2013_2015_cat1_freq
OR_H1N1_2013_2015_cat1_Bayes


### Analyses of interest: Age-stratified H1N1 2015 to 2018 ----
try9 <- data_orig %>%
  filter(at2015==1 & at2018==1 & H1N1pdm_2016==0 & H1N1pdm_2017==0) %>%
  dplyr::select(Obs, H1N1pdm_2015, H1N1pdm_2018, age_2015_cat)

try9$H1N1pdm_2015 <- as.numeric(try9$H1N1pdm_2015)
try9$H1N1pdm_2018 <- as.numeric(try9$H1N1pdm_2018)

# 2x2 table
table(try9[try9$age_2015_cat==0,]$H1N1pdm_2015, try9[try9$age_2015_cat==0,]$H1N1pdm_2018) 
table(try9[try9$age_2015_cat==1,]$H1N1pdm_2015, try9[try9$age_2015_cat==1,]$H1N1pdm_2018) 


# Bayesian analysis of H1N1 2015 on H1N1 2018, under 5s
H1N1pdm_2015 <- try9[try9$age_2015_cat==0,]$H1N1pdm_2015
H1N1pdm_2018 <- try9[try9$age_2015_cat==0,]$H1N1pdm_2018

logistic.model <- function() {
  # SAMPLING DISTRIBUTION
  for (i in 1:N) {
    logit(pi[i]) <- b[1] + b[2]*H1N1pdm_2015[i];
    H1N1pdm_2018[i] ~ dbin(pi[i],1);  
  }
  
  # PRIORS ON BETAS
  b[1:Nx] ~ dmnorm(mu[1:Nx], tau[1:Nx,1:Nx]) # multivariate normal prior
  
  # Calculate ORs:
  for (l in 1:Nx) {
    OR[l] <- exp(b[l]);
  }
}

N <- length(H1N1pdm_2018) 
Nx <- 2 

mu <- rep(0, Nx)
tau <- diag(0.001, Nx)
data.logistic <- list("N", "Nx", "H1N1pdm_2015", "H1N1pdm_2018", "mu", "tau")
parameters.logistic <-c("b","OR") # Parameters to keep track of
inits.logistic <- function() {list (b= rep(0, Nx))}

set.seed(114011)
logistic.sim_2015_2018_cat0<-jags(data=data.logistic, inits=inits.logistic, parameters.logistic, n.iter=100000,
                                  n.burn=25000, model.file=logistic.model, n.thin=5, n.chains = 3)

print(logistic.sim_2015_2018_cat0, 2)
OR_H1N1_2015_2018_cat0_Bayes <- round(logistic.sim_2015_2018_cat0$BUGSoutput$summary[2,c(5,3,7)], 2)
OR_H1N1_2015_2018_cat0_Bayes

# Bayesian analysis of H1N1 2015 on H1N1 2018, over 5s
H1N1pdm_2015 <- try9[try9$age_2015_cat==1,]$H1N1pdm_2015
H1N1pdm_2018 <- try9[try9$age_2015_cat==1,]$H1N1pdm_2018

logistic.model <- function() {
  # SAMPLING DISTRIBUTION
  for (i in 1:N) {
    logit(pi[i]) <- b[1] + b[2]*H1N1pdm_2015[i];
    H1N1pdm_2018[i] ~ dbin(pi[i],1);  
  }
  
  # PRIORS ON BETAS
  b[1:Nx] ~ dmnorm(mu[1:Nx], tau[1:Nx,1:Nx]) # multivariate normal prior
  
  # Calculate ORs:
  for (l in 1:Nx) {
    OR[l] <- exp(b[l]);
  }
}

N <- length(H1N1pdm_2018) 
Nx <- 2 

mu <- rep(0, Nx)
tau <- diag(0.001, Nx)
data.logistic <- list("N", "Nx", "H1N1pdm_2015", "H1N1pdm_2018", "mu", "tau")
parameters.logistic <-c("b","OR") # Parameters to keep track of
inits.logistic <- function() {list (b= rep(0, Nx))}

set.seed(114011)
logistic.sim_2015_2018_cat1<-jags(data=data.logistic, inits=inits.logistic, parameters.logistic, n.iter=100000,
                                  n.burn=25000, model.file=logistic.model, n.thin=5, n.chains = 3)

print(logistic.sim_2015_2018_cat1, 2)
OR_H1N1_2015_2018_cat1_Bayes <- round(logistic.sim_2015_2018_cat1$BUGSoutput$summary[2,c(5,3,7)], 2)
OR_H1N1_2015_2018_cat1_Bayes

### sensitivity analysis

try9 <- data_orig %>%
  filter(at2015==1 & at2018==1 ) %>%
  dplyr::select(Obs, H1N1pdm_2015, H1N1pdm_2018, age_2015_cat)

try9$H1N1pdm_2015 <- as.numeric(try9$H1N1pdm_2015)
try9$H1N1pdm_2018 <- as.numeric(try9$H1N1pdm_2018)

# 2x2 table
table(try9[try9$age_2015_cat==0,]$H1N1pdm_2015, try9[try9$age_2015_cat==0,]$H1N1pdm_2018) 
table(try9[try9$age_2015_cat==1,]$H1N1pdm_2015, try9[try9$age_2015_cat==1,]$H1N1pdm_2018) 


# Bayesian analysis of H1N1 2015 on H1N1 2018, under 5s
H1N1pdm_2015 <- try9[try9$age_2015_cat==0,]$H1N1pdm_2015
H1N1pdm_2018 <- try9[try9$age_2015_cat==0,]$H1N1pdm_2018

logistic.model <- function() {
  # SAMPLING DISTRIBUTION
  for (i in 1:N) {
    logit(pi[i]) <- b[1] + b[2]*H1N1pdm_2015[i];
    H1N1pdm_2018[i] ~ dbin(pi[i],1);  
  }
  
  # PRIORS ON BETAS
  b[1:Nx] ~ dmnorm(mu[1:Nx], tau[1:Nx,1:Nx]) # multivariate normal prior
  
  # Calculate ORs:
  for (l in 1:Nx) {
    OR[l] <- exp(b[l]);
  }
}

N <- length(H1N1pdm_2018) 
Nx <- 2 

mu <- rep(0, Nx)
tau <- diag(0.001, Nx)
data.logistic <- list("N", "Nx", "H1N1pdm_2015", "H1N1pdm_2018", "mu", "tau")
parameters.logistic <-c("b","OR") # Parameters to keep track of
inits.logistic <- function() {list (b= rep(0, Nx))}

set.seed(114011)
logistic.sim_2015_2018_cat0<-jags(data=data.logistic, inits=inits.logistic, parameters.logistic, n.iter=100000,
                                  n.burn=25000, model.file=logistic.model, n.thin=5, n.chains = 3)

print(logistic.sim_2015_2018_cat0, 2)
OR_H1N1_2015_2018_cat0_Bayes <- round(logistic.sim_2015_2018_cat0$BUGSoutput$summary[2,c(5,3,7)], 2)
OR_H1N1_2015_2018_cat0_Bayes

# Bayesian analysis of H1N1 2015 on H1N1 2018, over 5s
H1N1pdm_2015 <- try9[try9$age_2015_cat==1,]$H1N1pdm_2015
H1N1pdm_2018 <- try9[try9$age_2015_cat==1,]$H1N1pdm_2018

logistic.model <- function() {
  # SAMPLING DISTRIBUTION
  for (i in 1:N) {
    logit(pi[i]) <- b[1] + b[2]*H1N1pdm_2015[i];
    H1N1pdm_2018[i] ~ dbin(pi[i],1);  
  }
  
  # PRIORS ON BETAS
  b[1:Nx] ~ dmnorm(mu[1:Nx], tau[1:Nx,1:Nx]) # multivariate normal prior
  
  # Calculate ORs:
  for (l in 1:Nx) {
    OR[l] <- exp(b[l]);
  }
}

N <- length(H1N1pdm_2018) 
Nx <- 2 

mu <- rep(0, Nx)
tau <- diag(0.001, Nx)
data.logistic <- list("N", "Nx", "H1N1pdm_2015", "H1N1pdm_2018", "mu", "tau")
parameters.logistic <-c("b","OR") # Parameters to keep track of
inits.logistic <- function() {list (b= rep(0, Nx))}

set.seed(114011)
logistic.sim_2015_2018_cat1<-jags(data=data.logistic, inits=inits.logistic, parameters.logistic, n.iter=100000,
                                  n.burn=25000, model.file=logistic.model, n.thin=5, n.chains = 3)

print(logistic.sim_2015_2018_cat1, 2)
OR_H1N1_2015_2018_cat1_Bayes <- round(logistic.sim_2015_2018_cat1$BUGSoutput$summary[2,c(5,3,7)], 2)
OR_H1N1_2015_2018_cat1_Bayes


### Analyses of interest: Age-stratified H1N1 2011 to 2015 ----
try9 <- data_orig %>%
  filter(at2011==1 & at2015==1 & H1N1pdm_2013==0) %>%
  dplyr::select(Obs, H1N1pdm_2011, H1N1pdm_2015, age_2011_cat)

try9$H1N1pdm_2011 <- as.numeric(try9$H1N1pdm_2011)
try9$H1N1pdm_2015 <- as.numeric(try9$H1N1pdm_2015)

# 2x2 table
table(try9[try9$age_2011_cat==0,]$H1N1pdm_2011, try9[try9$age_2011_cat==0,]$H1N1pdm_2015) 
table(try9[try9$age_2011_cat==1,]$H1N1pdm_2011, try9[try9$age_2011_cat==1,]$H1N1pdm_2015) 

# Bayesian analysis of H1N1 2015 on H1N1 2018, under 5s
H1N1pdm_2011 <- try9[try9$age_2011_cat==0,]$H1N1pdm_2011
H1N1pdm_2015 <- try9[try9$age_2011_cat==0,]$H1N1pdm_2015

logistic.model <- function() {
  # SAMPLING DISTRIBUTION
  for (i in 1:N) {
    logit(pi[i]) <- b[1] + b[2]*H1N1pdm_2011[i];
    H1N1pdm_2015[i] ~ dbin(pi[i],1);  
  }
  
  # PRIORS ON BETAS
  b[1:Nx] ~ dmnorm(mu[1:Nx], tau[1:Nx,1:Nx]) # multivariate normal prior
  
  # Calculate ORs:
  for (l in 1:Nx) {
    OR[l] <- exp(b[l]);
  }
}

N <- length(H1N1pdm_2015) 
Nx <- 2 

mu <- rep(0, Nx)
tau <- diag(0.001, Nx)
data.logistic <- list("N", "Nx", "H1N1pdm_2011", "H1N1pdm_2015", "mu", "tau")
parameters.logistic <-c("b","OR") # Parameters to keep track of
inits.logistic <- function() {list (b= rep(0, Nx))}

set.seed(114011)
logistic.sim_2011_2015_cat0<-jags(data=data.logistic, inits=inits.logistic, parameters.logistic, n.iter=100000,
                                  n.burn=25000, model.file=logistic.model, n.thin=5, n.chains = 3)

print(logistic.sim_2011_2015_cat0, 2)
OR_H1N1_2011_2015_cat0_Bayes <- round(logistic.sim_2011_2015_cat0$BUGSoutput$summary[2,c(5,3,7)], 2)
OR_H1N1_2011_2015_cat0_Bayes

# Bayesian analysis of H1N1 2015 on H1N1 2018, over 5s
H1N1pdm_2011 <- try9[try9$age_2011_cat==1,]$H1N1pdm_2011
H1N1pdm_2015 <- try9[try9$age_2011_cat==1,]$H1N1pdm_2015

logistic.model <- function() {
  # SAMPLING DISTRIBUTION
  for (i in 1:N) {
    logit(pi[i]) <- b[1] + b[2]*H1N1pdm_2011[i];
    H1N1pdm_2015[i] ~ dbin(pi[i],1);  
  }
  
  # PRIORS ON BETAS
  b[1:Nx] ~ dmnorm(mu[1:Nx], tau[1:Nx,1:Nx]) # multivariate normal prior
  
  # Calculate ORs:
  for (l in 1:Nx) {
    OR[l] <- exp(b[l]);
  }
}

N <- length(H1N1pdm_2015) 
Nx <- 2 

mu <- rep(0, Nx)
tau <- diag(0.001, Nx)
data.logistic <- list("N", "Nx", "H1N1pdm_2011", "H1N1pdm_2015", "mu", "tau")
parameters.logistic <-c("b","OR") # Parameters to keep track of
inits.logistic <- function() {list (b= rep(0, Nx))}

set.seed(114011)
logistic.sim_2011_2015_cat1<-jags(data=data.logistic, inits=inits.logistic, parameters.logistic, n.iter=100000,
                                  n.burn=25000, model.file=logistic.model, n.thin=5, n.chains = 3)

print(logistic.sim_2011_2015_cat1, 2)
OR_H1N1_2011_2015_cat1_Bayes <- round(logistic.sim_2011_2015_cat1$BUGSoutput$summary[2,c(5,3,7)], 2)
OR_H1N1_2011_2015_cat1_Bayes


### SENSITIVITY ANALYSIS

try9 <- data_orig %>%
  filter(at2011==1 & at2015==1) %>%
  dplyr::select(Obs, H1N1pdm_2011, H1N1pdm_2015, age_2011_cat)

try9$H1N1pdm_2011 <- as.numeric(try9$H1N1pdm_2011)
try9$H1N1pdm_2015 <- as.numeric(try9$H1N1pdm_2015)

# 2x2 table
table(try9[try9$age_2011_cat==0,]$H1N1pdm_2011, try9[try9$age_2011_cat==0,]$H1N1pdm_2015) 
table(try9[try9$age_2011_cat==1,]$H1N1pdm_2011, try9[try9$age_2011_cat==1,]$H1N1pdm_2015) 

# Bayesian analysis of H1N1 2015 on H1N1 2018, under 5s
H1N1pdm_2011 <- try9[try9$age_2011_cat==0,]$H1N1pdm_2011
H1N1pdm_2015 <- try9[try9$age_2011_cat==0,]$H1N1pdm_2015

logistic.model <- function() {
  # SAMPLING DISTRIBUTION
  for (i in 1:N) {
    logit(pi[i]) <- b[1] + b[2]*H1N1pdm_2011[i];
    H1N1pdm_2015[i] ~ dbin(pi[i],1);  
  }
  
  # PRIORS ON BETAS
  b[1:Nx] ~ dmnorm(mu[1:Nx], tau[1:Nx,1:Nx]) # multivariate normal prior
  
  # Calculate ORs:
  for (l in 1:Nx) {
    OR[l] <- exp(b[l]);
  }
}

N <- length(H1N1pdm_2015) 
Nx <- 2 

mu <- rep(0, Nx)
tau <- diag(0.001, Nx)
data.logistic <- list("N", "Nx", "H1N1pdm_2011", "H1N1pdm_2015", "mu", "tau")
parameters.logistic <-c("b","OR") # Parameters to keep track of
inits.logistic <- function() {list (b= rep(0, Nx))}

set.seed(114011)
logistic.sim_2011_2015_cat0<-jags(data=data.logistic, inits=inits.logistic, parameters.logistic, n.iter=100000,
                                  n.burn=25000, model.file=logistic.model, n.thin=5, n.chains = 3)

print(logistic.sim_2011_2015_cat0, 2)
OR_H1N1_2011_2015_cat0_Bayes <- round(logistic.sim_2011_2015_cat0$BUGSoutput$summary[2,c(5,3,7)], 2)
OR_H1N1_2011_2015_cat0_Bayes

# Bayesian analysis of H1N1 2015 on H1N1 2018, over 5s
H1N1pdm_2011 <- try9[try9$age_2011_cat==1,]$H1N1pdm_2011
H1N1pdm_2015 <- try9[try9$age_2011_cat==1,]$H1N1pdm_2015

logistic.model <- function() {
  # SAMPLING DISTRIBUTION
  for (i in 1:N) {
    logit(pi[i]) <- b[1] + b[2]*H1N1pdm_2011[i];
    H1N1pdm_2015[i] ~ dbin(pi[i],1);  
  }
  
  # PRIORS ON BETAS
  b[1:Nx] ~ dmnorm(mu[1:Nx], tau[1:Nx,1:Nx]) # multivariate normal prior
  
  # Calculate ORs:
  for (l in 1:Nx) {
    OR[l] <- exp(b[l]);
  }
}

N <- length(H1N1pdm_2015) 
Nx <- 2 

mu <- rep(0, Nx)
tau <- diag(0.001, Nx)
data.logistic <- list("N", "Nx", "H1N1pdm_2011", "H1N1pdm_2015", "mu", "tau")
parameters.logistic <-c("b","OR") # Parameters to keep track of
inits.logistic <- function() {list (b= rep(0, Nx))}

set.seed(114011)
logistic.sim_2011_2015_cat1<-jags(data=data.logistic, inits=inits.logistic, parameters.logistic, n.iter=100000,
                                  n.burn=25000, model.file=logistic.model, n.thin=5, n.chains = 3)

print(logistic.sim_2011_2015_cat1, 2)
OR_H1N1_2011_2015_cat1_Bayes <- round(logistic.sim_2011_2015_cat1$BUGSoutput$summary[2,c(5,3,7)], 2)
OR_H1N1_2011_2015_cat1_Bayes


### Analyses of interest: Age-stratified H1N1 2011 to 2018 ----
try9 <- data_orig %>%
  filter(at2011==1 & at2018==1 & H1N1pdm_2013==0 & H1N1pdm_2015==0 & H1N1pdm_2016==0 & H1N1pdm_2017==0) %>%
  dplyr::select(Obs, H1N1pdm_2011, H1N1pdm_2018, age_2011_cat)

try9$H1N1pdm_2011 <- as.numeric(try9$H1N1pdm_2011)
try9$H1N1pdm_2018 <- as.numeric(try9$H1N1pdm_2018)

# 2x2 table
table(try9[try9$age_2011_cat==0,]$H1N1pdm_2011, try9[try9$age_2011_cat==0,]$H1N1pdm_2018) 
table(try9[try9$age_2011_cat==1,]$H1N1pdm_2011, try9[try9$age_2011_cat==1,]$H1N1pdm_2018) 

# Bayesian analysis of H1N1 2015 on H1N1 2018, under 5s
H1N1pdm_2011 <- try9[try9$age_2011_cat==0,]$H1N1pdm_2011
H1N1pdm_2018 <- try9[try9$age_2011_cat==0,]$H1N1pdm_2018

logistic.model <- function() {
  # SAMPLING DISTRIBUTION
  for (i in 1:N) {
    logit(pi[i]) <- b[1] + b[2]*H1N1pdm_2011[i];
    H1N1pdm_2018[i] ~ dbin(pi[i],1);  
  }
  
  # PRIORS ON BETAS
  b[1:Nx] ~ dmnorm(mu[1:Nx], tau[1:Nx,1:Nx]) # multivariate normal prior
  
  # Calculate ORs:
  for (l in 1:Nx) {
    OR[l] <- exp(b[l]);
  }
}

N <- length(H1N1pdm_2018) 
Nx <- 2 

mu <- rep(0, Nx)
tau <- diag(0.001, Nx)
data.logistic <- list("N", "Nx", "H1N1pdm_2011", "H1N1pdm_2018", "mu", "tau")
parameters.logistic <-c("b","OR") # Parameters to keep track of
inits.logistic <- function() {list (b= rep(0, Nx))}

set.seed(114011)
logistic.sim_2011_2018_cat0<-jags(data=data.logistic, inits=inits.logistic, parameters.logistic, n.iter=100000,
                                  n.burn=25000, model.file=logistic.model, n.thin=5, n.chains = 3)

print(logistic.sim_2011_2018_cat0, 2)
OR_H1N1_2011_2018_cat0_Bayes <- round(logistic.sim_2011_2018_cat0$BUGSoutput$summary[2,c(5,3,7)], 2)
OR_H1N1_2011_2018_cat0_Bayes

# Bayesian analysis of H1N1 2015 on H1N1 2018, over 5s
H1N1pdm_2011 <- try9[try9$age_2011_cat==1,]$H1N1pdm_2011
H1N1pdm_2018 <- try9[try9$age_2011_cat==1,]$H1N1pdm_2018

logistic.model <- function() {
  # SAMPLING DISTRIBUTION
  for (i in 1:N) {
    logit(pi[i]) <- b[1] + b[2]*H1N1pdm_2011[i];
    H1N1pdm_2018[i] ~ dbin(pi[i],1);  
  }
  
  # PRIORS ON BETAS
  b[1:Nx] ~ dmnorm(mu[1:Nx], tau[1:Nx,1:Nx]) # multivariate normal prior
  
  # Calculate ORs:
  for (l in 1:Nx) {
    OR[l] <- exp(b[l]);
  }
}

N <- length(H1N1pdm_2018) 
Nx <- 2 

mu <- rep(0, Nx)
tau <- diag(0.001, Nx)
data.logistic <- list("N", "Nx", "H1N1pdm_2011", "H1N1pdm_2018", "mu", "tau")
parameters.logistic <-c("b","OR") # Parameters to keep track of
inits.logistic <- function() {list (b= rep(0, Nx))}

set.seed(114011)
logistic.sim_2011_2018_cat1<-jags(data=data.logistic, inits=inits.logistic, parameters.logistic, n.iter=100000,
                                  n.burn=25000, model.file=logistic.model, n.thin=5, n.chains = 3)

print(logistic.sim_2011_2018_cat1, 2)
OR_H1N1_2011_2018_cat1_Bayes <- round(logistic.sim_2011_2018_cat1$BUGSoutput$summary[2,c(5,3,7)], 2)
OR_H1N1_2011_2018_cat1_Bayes


### SENSITIVITY ANALYSIS
try9 <- data_orig %>%
  filter(at2011==1 & at2018==1) %>%
  dplyr::select(Obs, H1N1pdm_2011, H1N1pdm_2018, age_2011_cat)

try9$H1N1pdm_2011 <- as.numeric(try9$H1N1pdm_2011)
try9$H1N1pdm_2018 <- as.numeric(try9$H1N1pdm_2018)

# 2x2 table
table(try9[try9$age_2011_cat==0,]$H1N1pdm_2011, try9[try9$age_2011_cat==0,]$H1N1pdm_2018) 
table(try9[try9$age_2011_cat==1,]$H1N1pdm_2011, try9[try9$age_2011_cat==1,]$H1N1pdm_2018) 

# Bayesian analysis of H1N1 2015 on H1N1 2018, under 5s
H1N1pdm_2011 <- try9[try9$age_2011_cat==0,]$H1N1pdm_2011
H1N1pdm_2018 <- try9[try9$age_2011_cat==0,]$H1N1pdm_2018

logistic.model <- function() {
  # SAMPLING DISTRIBUTION
  for (i in 1:N) {
    logit(pi[i]) <- b[1] + b[2]*H1N1pdm_2011[i];
    H1N1pdm_2018[i] ~ dbin(pi[i],1);  
  }
  
  # PRIORS ON BETAS
  b[1:Nx] ~ dmnorm(mu[1:Nx], tau[1:Nx,1:Nx]) # multivariate normal prior
  
  # Calculate ORs:
  for (l in 1:Nx) {
    OR[l] <- exp(b[l]);
  }
}

N <- length(H1N1pdm_2018) 
Nx <- 2 

mu <- rep(0, Nx)
tau <- diag(0.001, Nx)
data.logistic <- list("N", "Nx", "H1N1pdm_2011", "H1N1pdm_2018", "mu", "tau")
parameters.logistic <-c("b","OR") # Parameters to keep track of
inits.logistic <- function() {list (b= rep(0, Nx))}

set.seed(114011)
logistic.sim_2011_2018_cat0<-jags(data=data.logistic, inits=inits.logistic, parameters.logistic, n.iter=100000,
                                  n.burn=25000, model.file=logistic.model, n.thin=5, n.chains = 3)

print(logistic.sim_2011_2018_cat0, 2)
OR_H1N1_2011_2018_cat0_Bayes <- round(logistic.sim_2011_2018_cat0$BUGSoutput$summary[2,c(5,3,7)], 2)
OR_H1N1_2011_2018_cat0_Bayes

# Bayesian analysis of H1N1 2015 on H1N1 2018, over 5s
H1N1pdm_2011 <- try9[try9$age_2011_cat==1,]$H1N1pdm_2011
H1N1pdm_2018 <- try9[try9$age_2011_cat==1,]$H1N1pdm_2018

logistic.model <- function() {
  # SAMPLING DISTRIBUTION
  for (i in 1:N) {
    logit(pi[i]) <- b[1] + b[2]*H1N1pdm_2011[i];
    H1N1pdm_2018[i] ~ dbin(pi[i],1);  
  }
  
  # PRIORS ON BETAS
  b[1:Nx] ~ dmnorm(mu[1:Nx], tau[1:Nx,1:Nx]) # multivariate normal prior
  
  # Calculate ORs:
  for (l in 1:Nx) {
    OR[l] <- exp(b[l]);
  }
}

N <- length(H1N1pdm_2018) 
Nx <- 2 

mu <- rep(0, Nx)
tau <- diag(0.001, Nx)
data.logistic <- list("N", "Nx", "H1N1pdm_2011", "H1N1pdm_2018", "mu", "tau")
parameters.logistic <-c("b","OR") # Parameters to keep track of
inits.logistic <- function() {list (b= rep(0, Nx))}

set.seed(114011)
logistic.sim_2011_2018_cat1<-jags(data=data.logistic, inits=inits.logistic, parameters.logistic, n.iter=100000,
                                  n.burn=25000, model.file=logistic.model, n.thin=5, n.chains = 3)

print(logistic.sim_2011_2018_cat1, 2)
OR_H1N1_2011_2018_cat1_Bayes <- round(logistic.sim_2011_2018_cat1$BUGSoutput$summary[2,c(5,3,7)], 2)
OR_H1N1_2011_2018_cat1_Bayes


### Analyses of interest: Age-stratified H1N1 2013 to 2018 ----
try9 <- data_orig %>%
  filter(at2013==1 & at2018==1 & H1N1pdm_2015==0 & H1N1pdm_2016==0 & H1N1pdm_2017==0) %>%
  dplyr::select(Obs, H1N1pdm_2013, H1N1pdm_2018, age_2013_cat)

try9$H1N1pdm_2013 <- as.numeric(try9$H1N1pdm_2013)
try9$H1N1pdm_2018 <- as.numeric(try9$H1N1pdm_2018)

# 2x2 table
table(try9[try9$age_2013_cat==0,]$H1N1pdm_2013, try9[try9$age_2013_cat==0,]$H1N1pdm_2018) 
table(try9[try9$age_2013_cat==1,]$H1N1pdm_2013, try9[try9$age_2013_cat==1,]$H1N1pdm_2018) 

# Bayesian analysis of H1N1 2015 on H1N1 2018, under 5s
H1N1pdm_2013 <- try9[try9$age_2013_cat==0,]$H1N1pdm_2013
H1N1pdm_2018 <- try9[try9$age_2013_cat==0,]$H1N1pdm_2018

logistic.model <- function() {
  # SAMPLING DISTRIBUTION
  for (i in 1:N) {
    logit(pi[i]) <- b[1] + b[2]*H1N1pdm_2013[i];
    H1N1pdm_2018[i] ~ dbin(pi[i],1);  
  }
  
  # PRIORS ON BETAS
  b[1:Nx] ~ dmnorm(mu[1:Nx], tau[1:Nx,1:Nx]) # multivariate normal prior
  
  # Calculate ORs:
  for (l in 1:Nx) {
    OR[l] <- exp(b[l]);
  }
}

N <- length(H1N1pdm_2018) 
Nx <- 2 

mu <- rep(0, Nx)
tau <- diag(0.001, Nx)
data.logistic <- list("N", "Nx", "H1N1pdm_2013", "H1N1pdm_2018", "mu", "tau")
parameters.logistic <-c("b","OR") # Parameters to keep track of
inits.logistic <- function() {list (b= rep(0, Nx))}

set.seed(114011)
logistic.sim_2013_2018_cat0<-jags(data=data.logistic, inits=inits.logistic, parameters.logistic, n.iter=100000,
                                  n.burn=25000, model.file=logistic.model, n.thin=5, n.chains = 3)

print(logistic.sim_2011_2018_cat0, 2)
OR_H1N1_2013_2018_cat0_Bayes <- round(logistic.sim_2013_2018_cat0$BUGSoutput$summary[2,c(5,3,7)], 2)
OR_H1N1_2013_2018_cat0_Bayes

# Bayesian analysis of H1N1 2015 on H1N1 2018, over 5s
H1N1pdm_2013 <- try9[try9$age_2013_cat==1,]$H1N1pdm_2013
H1N1pdm_2018 <- try9[try9$age_2013_cat==1,]$H1N1pdm_2018

logistic.model <- function() {
  # SAMPLING DISTRIBUTION
  for (i in 1:N) {
    logit(pi[i]) <- b[1] + b[2]*H1N1pdm_2013[i];
    H1N1pdm_2018[i] ~ dbin(pi[i],1);  
  }
  
  # PRIORS ON BETAS
  b[1:Nx] ~ dmnorm(mu[1:Nx], tau[1:Nx,1:Nx]) # multivariate normal prior
  
  # Calculate ORs:
  for (l in 1:Nx) {
    OR[l] <- exp(b[l]);
  }
}

N <- length(H1N1pdm_2018) 
Nx <- 2 

mu <- rep(0, Nx)
tau <- diag(0.001, Nx)
data.logistic <- list("N", "Nx", "H1N1pdm_2013", "H1N1pdm_2018", "mu", "tau")
parameters.logistic <-c("b","OR") # Parameters to keep track of
inits.logistic <- function() {list (b= rep(0, Nx))}

set.seed(114011)
logistic.sim_2013_2018_cat1<-jags(data=data.logistic, inits=inits.logistic, parameters.logistic, n.iter=100000,
                                  n.burn=25000, model.file=logistic.model, n.thin=5, n.chains = 3)

print(logistic.sim_2013_2018_cat1, 2)
OR_H1N1_2013_2018_cat1_Bayes <- round(logistic.sim_2013_2018_cat1$BUGSoutput$summary[2,c(5,3,7)], 2)
OR_H1N1_2013_2018_cat1_Bayes


### SENSITIVITY ANALYSIS

try9 <- data_orig %>%
  filter(at2013==1 & at2018==1) %>%
  dplyr::select(Obs, H1N1pdm_2013, H1N1pdm_2018, age_2013_cat)

try9$H1N1pdm_2013 <- as.numeric(try9$H1N1pdm_2013)
try9$H1N1pdm_2018 <- as.numeric(try9$H1N1pdm_2018)

# 2x2 table
table(try9[try9$age_2013_cat==0,]$H1N1pdm_2013, try9[try9$age_2013_cat==0,]$H1N1pdm_2018) 
table(try9[try9$age_2013_cat==1,]$H1N1pdm_2013, try9[try9$age_2013_cat==1,]$H1N1pdm_2018) 

# Bayesian analysis of H1N1 2015 on H1N1 2018, under 5s
H1N1pdm_2013 <- try9[try9$age_2013_cat==0,]$H1N1pdm_2013
H1N1pdm_2018 <- try9[try9$age_2013_cat==0,]$H1N1pdm_2018

logistic.model <- function() {
  # SAMPLING DISTRIBUTION
  for (i in 1:N) {
    logit(pi[i]) <- b[1] + b[2]*H1N1pdm_2013[i];
    H1N1pdm_2018[i] ~ dbin(pi[i],1);  
  }
  
  # PRIORS ON BETAS
  b[1:Nx] ~ dmnorm(mu[1:Nx], tau[1:Nx,1:Nx]) # multivariate normal prior
  
  # Calculate ORs:
  for (l in 1:Nx) {
    OR[l] <- exp(b[l]);
  }
}

N <- length(H1N1pdm_2018) 
Nx <- 2 

mu <- rep(0, Nx)
tau <- diag(0.001, Nx)
data.logistic <- list("N", "Nx", "H1N1pdm_2013", "H1N1pdm_2018", "mu", "tau")
parameters.logistic <-c("b","OR") # Parameters to keep track of
inits.logistic <- function() {list (b= rep(0, Nx))}

set.seed(114011)
logistic.sim_2013_2018_cat0<-jags(data=data.logistic, inits=inits.logistic, parameters.logistic, n.iter=100000,
                                  n.burn=25000, model.file=logistic.model, n.thin=5, n.chains = 3)

print(logistic.sim_2011_2018_cat0, 2)
OR_H1N1_2013_2018_cat0_Bayes <- round(logistic.sim_2013_2018_cat0$BUGSoutput$summary[2,c(5,3,7)], 2)
OR_H1N1_2013_2018_cat0_Bayes

# Bayesian analysis of H1N1 2015 on H1N1 2018, over 5s
H1N1pdm_2013 <- try9[try9$age_2013_cat==1,]$H1N1pdm_2013
H1N1pdm_2018 <- try9[try9$age_2013_cat==1,]$H1N1pdm_2018

logistic.model <- function() {
  # SAMPLING DISTRIBUTION
  for (i in 1:N) {
    logit(pi[i]) <- b[1] + b[2]*H1N1pdm_2013[i];
    H1N1pdm_2018[i] ~ dbin(pi[i],1);  
  }
  
  # PRIORS ON BETAS
  b[1:Nx] ~ dmnorm(mu[1:Nx], tau[1:Nx,1:Nx]) # multivariate normal prior
  
  # Calculate ORs:
  for (l in 1:Nx) {
    OR[l] <- exp(b[l]);
  }
}

N <- length(H1N1pdm_2018) 
Nx <- 2 

mu <- rep(0, Nx)
tau <- diag(0.001, Nx)
data.logistic <- list("N", "Nx", "H1N1pdm_2013", "H1N1pdm_2018", "mu", "tau")
parameters.logistic <-c("b","OR") # Parameters to keep track of
inits.logistic <- function() {list (b= rep(0, Nx))}

set.seed(114011)
logistic.sim_2013_2018_cat1<-jags(data=data.logistic, inits=inits.logistic, parameters.logistic, n.iter=100000,
                                  n.burn=25000, model.file=logistic.model, n.thin=5, n.chains = 3)

print(logistic.sim_2013_2018_cat1, 2)
OR_H1N1_2013_2018_cat1_Bayes <- round(logistic.sim_2013_2018_cat1$BUGSoutput$summary[2,c(5,3,7)], 2)
OR_H1N1_2013_2018_cat1_Bayes

