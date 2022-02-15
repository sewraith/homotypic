## Analyses for IB by lineage

library(R2jags)
library(coda)
library(lmtest)
library(foreign)
library(plyr)
library(dplyr)
library(varhandle)

setwd("/Users/sewraith/Box/Research/Homotypic/Data/tempwork/Repeat")


### Analyses of interest: IB by Lineage ----

# Bring in the data
data_orig <- read.csv2("faustodata_updated6nov20.csv", header = T, sep = ",") 

# IB 2014 to 2017 - Yamagata
try52 <- data_orig %>%
  filter(at2014==1 & at2017==1 & IB_2016_Yamagata==0) %>%
  dplyr::select(Obs, IB_2014_Yamagata, IB_2017_Yamagata)
try52$IB_2014_Yamagata <- as.numeric(try52$IB_2014_Yamagata)
try52$IB_2017_Yamagata <- as.numeric(try52$IB_2017_Yamagata)

table(try52$IB_2014_Yamagata, try52$IB_2017_Yamagata)

# Bayesian analyses
# Putting exposure and outcome variables in the environment 
IB_2014_Yamagata <- try52$IB_2014_Yamagata
IB_2017_Yamagata <- try52$IB_2017_Yamagata

### Bayesian logistic regression
logistic.model <- function() {
  # SAMPLING DISTRIBUTION
  for (i in 1:N) {
    logit(pi[i]) <- b[1] + b[2]*IB_2014_Yamagata[i];
    IB_2017_Yamagata[i] ~ dbin(pi[i],1);  
  }
  
  # PRIORS ON BETAS
  b[1:Nx] ~ dmnorm(mu[1:Nx], tau[1:Nx,1:Nx]) # multivariate normal prior
  
  # Calculate ORs:
  for (l in 1:Nx) {
    OR[l] <- exp(b[l]);
  }
}

N <- length(IB_2017_Yamagata) # Number of observations
Nx <- 2 # Number of variables: intercept + exposure 

# Data, parameter list and starting values
mu <- rep(0, Nx)
tau <- diag(0.001, Nx)
data.logistic <- list("N", "Nx", "IB_2014_Yamagata", "IB_2017_Yamagata", "mu", "tau")
parameters.logistic <-c("b","OR") # Parameters to keep track of
inits.logistic <- function() {list (b= rep(0, Nx))}

set.seed(114011)
logistic.sim<-jags(data=data.logistic, inits=inits.logistic, parameters.logistic, n.iter=100000,
                   n.burn=25000, model.file=logistic.model, n.thin=5, n.chains = 3)

print(logistic.sim, 2)
OR_IB_2014_2017_Yamagata_Bayes_SA <- round(logistic.sim$BUGSoutput$summary[2,c(5,3,7)], 2)
OR_IB_2014_2017_Yamagata_Bayes_SA

### sensitivity analysis

try52 <- data_orig %>%
  filter(at2014==1 & at2017==1 ) %>%
  dplyr::select(Obs, IB_2014_Yamagata, IB_2017_Yamagata)
try52$IB_2014_Yamagata <- as.numeric(try52$IB_2014_Yamagata)
try52$IB_2017_Yamagata <- as.numeric(try52$IB_2017_Yamagata)

table(try52$IB_2014_Yamagata, try52$IB_2017_Yamagata)

# Bayesian analyses
# Putting exposure and outcome variables in the environment 
IB_2014_Yamagata <- try52$IB_2014_Yamagata
IB_2017_Yamagata <- try52$IB_2017_Yamagata

### Bayesian logistic regression
logistic.model <- function() {
  # SAMPLING DISTRIBUTION
  for (i in 1:N) {
    logit(pi[i]) <- b[1] + b[2]*IB_2014_Yamagata[i];
    IB_2017_Yamagata[i] ~ dbin(pi[i],1);  
  }
  
  # PRIORS ON BETAS
  b[1:Nx] ~ dmnorm(mu[1:Nx], tau[1:Nx,1:Nx]) # multivariate normal prior
  
  # Calculate ORs:
  for (l in 1:Nx) {
    OR[l] <- exp(b[l]);
  }
}

N <- length(IB_2017_Yamagata) # Number of observations
Nx <- 2 # Number of variables: intercept + exposure 

# Data, parameter list and starting values
mu <- rep(0, Nx)
tau <- diag(0.001, Nx)
data.logistic <- list("N", "Nx", "IB_2014_Yamagata", "IB_2017_Yamagata", "mu", "tau")
parameters.logistic <-c("b","OR") # Parameters to keep track of
inits.logistic <- function() {list (b= rep(0, Nx))}

set.seed(114011)
logistic.sim<-jags(data=data.logistic, inits=inits.logistic, parameters.logistic, n.iter=100000,
                   n.burn=25000, model.file=logistic.model, n.thin=5, n.chains = 3)

print(logistic.sim, 2)
OR_IB_2014_2017_Yamagata_Bayes_SA <- round(logistic.sim$BUGSoutput$summary[2,c(5,3,7)], 2)
OR_IB_2014_2017_Yamagata_Bayes_SA



# IB 2014 to 2019 - Yamagata
try52 <- data_orig %>%
  filter(at2014==1 & at2019==1 & IB_2016_Yamagata==0 & IB_2017_Yamagata==0 & IB_2018_Yamagata==0) %>%
  dplyr::select(Obs, IB_2014_Yamagata, IB_2019_Yamagata)
try52$IB_2014_Yamagata <- as.numeric(try52$IB_2014_Yamagata)
try52$IB_2019_Yamagata <- as.numeric(try52$IB_2019_Yamagata)

table(try52$IB_2014_Yamagata, try52$IB_2019_Yamagata)

# Bayesian analyses
# Putting exposure and outcome variables in the environment 
IB_2014_Yamagata <- try52$IB_2014_Yamagata
IB_2019_Yamagata <- try52$IB_2019_Yamagata

### Bayesian logistic regression
logistic.model <- function() {
  # SAMPLING DISTRIBUTION
  for (i in 1:N) {
    logit(pi[i]) <- b[1] + b[2]*IB_2014_Yamagata[i];
    IB_2019_Yamagata[i] ~ dbin(pi[i],1);  
  }
  
  # PRIORS ON BETAS
  b[1:Nx] ~ dmnorm(mu[1:Nx], tau[1:Nx,1:Nx]) # multivariate normal prior
  
  # Calculate ORs:
  for (l in 1:Nx) {
    OR[l] <- exp(b[l]);
  }
}

N <- length(IB_2019_Yamagata) # Number of observations
Nx <- 2 # Number of variables: intercept + exposure 

# Data, parameter list and starting values
mu <- rep(0, Nx)
tau <- diag(0.001, Nx)
data.logistic <- list("N", "Nx", "IB_2014_Yamagata", "IB_2019_Yamagata", "mu", "tau")
parameters.logistic <-c("b","OR") # Parameters to keep track of
inits.logistic <- function() {list (b= rep(0, Nx))}

set.seed(114011)
logistic.sim<-jags(data=data.logistic, inits=inits.logistic, parameters.logistic, n.iter=100000,
                   n.burn=25000, model.file=logistic.model, n.thin=5, n.chains = 3)

print(logistic.sim, 2)
OR_IB_2014_2019_Yamagata_Bayes_SA <- round(logistic.sim$BUGSoutput$summary[2,c(5,3,7)], 2)
OR_IB_2014_2019_Yamagata_Bayes_SA

### sensitivity analysis

try52 <- data_orig %>%
  filter(at2014==1 & at2019==1) %>%
  dplyr::select(Obs, IB_2014_Yamagata, IB_2019_Yamagata)
try52$IB_2014_Yamagata <- as.numeric(try52$IB_2014_Yamagata)
try52$IB_2019_Yamagata <- as.numeric(try52$IB_2019_Yamagata)

table(try52$IB_2014_Yamagata, try52$IB_2019_Yamagata)

# Bayesian analyses
# Putting exposure and outcome variables in the environment 
IB_2014_Yamagata <- try52$IB_2014_Yamagata
IB_2019_Yamagata <- try52$IB_2019_Yamagata

### Bayesian logistic regression
logistic.model <- function() {
  # SAMPLING DISTRIBUTION
  for (i in 1:N) {
    logit(pi[i]) <- b[1] + b[2]*IB_2014_Yamagata[i];
    IB_2019_Yamagata[i] ~ dbin(pi[i],1);  
  }
  
  # PRIORS ON BETAS
  b[1:Nx] ~ dmnorm(mu[1:Nx], tau[1:Nx,1:Nx]) # multivariate normal prior
  
  # Calculate ORs:
  for (l in 1:Nx) {
    OR[l] <- exp(b[l]);
  }
}

N <- length(IB_2019_Yamagata) # Number of observations
Nx <- 2 # Number of variables: intercept + exposure 

# Data, parameter list and starting values
mu <- rep(0, Nx)
tau <- diag(0.001, Nx)
data.logistic <- list("N", "Nx", "IB_2014_Yamagata", "IB_2019_Yamagata", "mu", "tau")
parameters.logistic <-c("b","OR") # Parameters to keep track of
inits.logistic <- function() {list (b= rep(0, Nx))}

set.seed(114011)
logistic.sim<-jags(data=data.logistic, inits=inits.logistic, parameters.logistic, n.iter=100000,
                   n.burn=25000, model.file=logistic.model, n.thin=5, n.chains = 3)

print(logistic.sim, 2)
OR_IB_2014_2019_Yamagata_Bayes_SA <- round(logistic.sim$BUGSoutput$summary[2,c(5,3,7)], 2)
OR_IB_2014_2019_Yamagata_Bayes_SA



# IB 2017 to 2019 - Yamagata
try52 <- data_orig %>%
  filter(at2017==1 & at2019==1 & IB_2018_Yamagata==0) %>%
  dplyr::select(Obs, IB_2017_Yamagata, IB_2019_Yamagata)
try52$IB_2017_Yamagata <- as.numeric(try52$IB_2017_Yamagata)
try52$IB_2019_Yamagata <- as.numeric(try52$IB_2019_Yamagata)

table(try52$IB_2017_Yamagata, try52$IB_2019_Yamagata)

# Bayesian analyses
# Putting exposure and outcome variables in the environment 
IB_2017_Yamagata <- try52$IB_2017_Yamagata
IB_2019_Yamagata <- try52$IB_2019_Yamagata

### Bayesian logistic regression
logistic.model <- function() {
  # SAMPLING DISTRIBUTION
  for (i in 1:N) {
    logit(pi[i]) <- b[1] + b[2]*IB_2017_Yamagata[i];
    IB_2019_Yamagata[i] ~ dbin(pi[i],1);  
  }
  
  # PRIORS ON BETAS
  b[1:Nx] ~ dmnorm(mu[1:Nx], tau[1:Nx,1:Nx]) # multivariate normal prior
  
  # Calculate ORs:
  for (l in 1:Nx) {
    OR[l] <- exp(b[l]);
  }
}

N <- length(IB_2019_Yamagata) # Number of observations
Nx <- 2 # Number of variables: intercept + exposure 

# Data, parameter list and starting values
mu <- rep(0, Nx)
tau <- diag(0.001, Nx)
data.logistic <- list("N", "Nx", "IB_2017_Yamagata", "IB_2019_Yamagata", "mu", "tau")
parameters.logistic <-c("b","OR") # Parameters to keep track of
inits.logistic <- function() {list (b= rep(0, Nx))}

set.seed(114011)
logistic.sim<-jags(data=data.logistic, inits=inits.logistic, parameters.logistic, n.iter=100000,
                   n.burn=25000, model.file=logistic.model, n.thin=5, n.chains = 3)

print(logistic.sim, 2)
OR_IB_2017_2019_Yamagata_Bayes_SA <- round(logistic.sim$BUGSoutput$summary[2,c(5,3,7)], 2)
OR_IB_2017_2019_Yamagata_Bayes_SA

### sensitivity analysis

try52 <- data_orig %>%
  filter(at2017==1 & at2019==1 ) %>%
  dplyr::select(Obs, IB_2017_Yamagata, IB_2019_Yamagata)
try52$IB_2017_Yamagata <- as.numeric(try52$IB_2017_Yamagata)
try52$IB_2019_Yamagata <- as.numeric(try52$IB_2019_Yamagata)

table(try52$IB_2017_Yamagata, try52$IB_2019_Yamagata)

# Bayesian analyses
# Putting exposure and outcome variables in the environment 
IB_2017_Yamagata <- try52$IB_2017_Yamagata
IB_2019_Yamagata <- try52$IB_2019_Yamagata

### Bayesian logistic regression
logistic.model <- function() {
  # SAMPLING DISTRIBUTION
  for (i in 1:N) {
    logit(pi[i]) <- b[1] + b[2]*IB_2017_Yamagata[i];
    IB_2019_Yamagata[i] ~ dbin(pi[i],1);  
  }
  
  # PRIORS ON BETAS
  b[1:Nx] ~ dmnorm(mu[1:Nx], tau[1:Nx,1:Nx]) # multivariate normal prior
  
  # Calculate ORs:
  for (l in 1:Nx) {
    OR[l] <- exp(b[l]);
  }
}

N <- length(IB_2019_Yamagata) # Number of observations
Nx <- 2 # Number of variables: intercept + exposure 

# Data, parameter list and starting values
mu <- rep(0, Nx)
tau <- diag(0.001, Nx)
data.logistic <- list("N", "Nx", "IB_2017_Yamagata", "IB_2019_Yamagata", "mu", "tau")
parameters.logistic <-c("b","OR") # Parameters to keep track of
inits.logistic <- function() {list (b= rep(0, Nx))}

set.seed(114011)
logistic.sim<-jags(data=data.logistic, inits=inits.logistic, parameters.logistic, n.iter=100000,
                   n.burn=25000, model.file=logistic.model, n.thin=5, n.chains = 3)

print(logistic.sim, 2)
OR_IB_2017_2019_Yamagata_Bayes_SA <- round(logistic.sim$BUGSoutput$summary[2,c(5,3,7)], 2)
OR_IB_2017_2019_Yamagata_Bayes_SA

# AGE STRATIFIED

# Yamagata 2017 to 2019

try35 <- data_orig %>%
  filter(at2017==1 & at2019==1 & IB_2018_Yamagata==0) %>%
  dplyr::select(Obs, IB_2017_Yamagata, IB_2019_Yamagata, age_2017_cat)

try35$IB_2017_Yamagata <- as.numeric(try35$IB_2017_Yamagata)
try35$IB_2019_Yamagata <- as.numeric(try35$IB_2019_Yamagata)

# 2x2 table
table(try35[try35$age_2017_cat==0,]$IB_2017_Yamagata, try35[try35$age_2017_cat==0,]$IB_2019_Yamagata) 
table(try35[try35$age_2017_cat==1,]$IB_2017_Yamagata, try35[try35$age_2017_cat==1,]$IB_2019_Yamagata) 

# Bayesian analysis of H3N2 2014 on H3N2 2016, under 5s
IB_2017_Yamagata <- try35[try35$age_2017_cat==0,]$IB_2017_Yamagata
IB_2019_Yamagata <- try35[try35$age_2017_cat==0,]$IB_2019_Yamagata

logistic.model <- function() {
  # SAMPLING DISTRIBUTION
  for (i in 1:N) {
    logit(pi[i]) <- b[1] + b[2]*IB_2017_Yamagata[i];
    IB_2019_Yamagata[i] ~ dbin(pi[i],1);  
  }
  
  # PRIORS ON BETAS
  b[1:Nx] ~ dmnorm(mu[1:Nx], tau[1:Nx,1:Nx]) # multivariate normal prior
  
  # Calculate ORs:
  for (l in 1:Nx) {
    OR[l] <- exp(b[l]);
  }
}

N <- length(IB_2019_Yamagata) 
Nx <- 2 

mu <- rep(0, Nx)
tau <- diag(0.001, Nx)
data.logistic <- list("N", "Nx", "IB_2017_Yamagata", "IB_2019_Yamagata", "mu", "tau")
parameters.logistic <-c("b","OR") # Parameters to keep track of
inits.logistic <- function() {list (b= rep(0, Nx))}

set.seed(114011)
logistic.sim_IB_Yamagata_2017_2019_cat0<-jags(data=data.logistic, inits=inits.logistic, parameters.logistic, n.iter=100000,
                                            n.burn=25000, model.file=logistic.model, n.thin=5, n.chains = 3)

print(logistic.sim_IB_Yamagata_2017_2019_cat0, 2)
OR_IB_Yamagata_2017_2019_cat0 <- round(logistic.sim_IB_Yamagata_2017_2019_cat0$BUGSoutput$summary[2,c(5,3,7)], 2)
OR_IB_Yamagata_2017_2019_cat0

# Bayesian analysis of H3N2 2014 on H3N2 2016, over 5s
IB_2017_Yamagata <- try35[try35$age_2017_cat==1,]$IB_2017_Yamagata
IB_2019_Yamagata <- try35[try35$age_2017_cat==1,]$IB_2019_Yamagata

logistic.model <- function() {
  # SAMPLING DISTRIBUTION
  for (i in 1:N) {
    logit(pi[i]) <- b[1] + b[2]*IB_2017_Yamagata[i];
    IB_2019_Yamagata[i] ~ dbin(pi[i],1);  
  }
  
  # PRIORS ON BETAS
  b[1:Nx] ~ dmnorm(mu[1:Nx], tau[1:Nx,1:Nx]) # multivariate normal prior
  
  # Calculate ORs:
  for (l in 1:Nx) {
    OR[l] <- exp(b[l]);
  }
}

N <- length(IB_2019_Yamagata) 
Nx <- 2 

mu <- rep(0, Nx)
tau <- diag(0.001, Nx)
data.logistic <- list("N", "Nx", "IB_2017_Yamagata", "IB_2019_Yamagata", "mu", "tau")
parameters.logistic <-c("b","OR") # Parameters to keep track of
inits.logistic <- function() {list (b= rep(0, Nx))}

set.seed(114011)
logistic.sim_IB_Yamagata_2017_2019_cat1<-jags(data=data.logistic, inits=inits.logistic, parameters.logistic, n.iter=100000,
                                            n.burn=25000, model.file=logistic.model, n.thin=5, n.chains = 3)

print(logistic.sim_IB_Yamagata_2017_2019_cat1, 2)
OR_IB_Yamagata_2017_2019_cat1 <- round(logistic.sim_IB_Yamagata_2017_2019_cat1$BUGSoutput$summary[2,c(5,3,7)], 2)
OR_IB_Yamagata_2017_2019_cat1

### sensitivity analysis


# Yamagata 2017 to 2019

try35 <- data_orig %>%
  filter(at2017==1 & at2019==1) %>%
  dplyr::select(Obs, IB_2017_Yamagata, IB_2019_Yamagata, age_2017_cat)

try35$IB_2017_Yamagata <- as.numeric(try35$IB_2017_Yamagata)
try35$IB_2019_Yamagata <- as.numeric(try35$IB_2019_Yamagata)

# 2x2 table
table(try35[try35$age_2017_cat==0,]$IB_2017_Yamagata, try35[try35$age_2017_cat==0,]$IB_2019_Yamagata) 
table(try35[try35$age_2017_cat==1,]$IB_2017_Yamagata, try35[try35$age_2017_cat==1,]$IB_2019_Yamagata) 


# Bayesian analysis of H3N2 2014 on H3N2 2016, under 5s
IB_2017_Yamagata <- try35[try35$age_2017_cat==0,]$IB_2017_Yamagata
IB_2019_Yamagata <- try35[try35$age_2017_cat==0,]$IB_2019_Yamagata

logistic.model <- function() {
  # SAMPLING DISTRIBUTION
  for (i in 1:N) {
    logit(pi[i]) <- b[1] + b[2]*IB_2017_Yamagata[i];
    IB_2019_Yamagata[i] ~ dbin(pi[i],1);  
  }
  
  # PRIORS ON BETAS
  b[1:Nx] ~ dmnorm(mu[1:Nx], tau[1:Nx,1:Nx]) # multivariate normal prior
  
  # Calculate ORs:
  for (l in 1:Nx) {
    OR[l] <- exp(b[l]);
  }
}

N <- length(IB_2019_Yamagata) 
Nx <- 2 

mu <- rep(0, Nx)
tau <- diag(0.001, Nx)
data.logistic <- list("N", "Nx", "IB_2017_Yamagata", "IB_2019_Yamagata", "mu", "tau")
parameters.logistic <-c("b","OR") # Parameters to keep track of
inits.logistic <- function() {list (b= rep(0, Nx))}

set.seed(114011)
logistic.sim_IB_Yamagata_2017_2019_cat0<-jags(data=data.logistic, inits=inits.logistic, parameters.logistic, n.iter=100000,
                                              n.burn=25000, model.file=logistic.model, n.thin=5, n.chains = 3)

print(logistic.sim_IB_Yamagata_2017_2019_cat0, 2)
OR_IB_Yamagata_2017_2019_cat0 <- round(logistic.sim_IB_Yamagata_2017_2019_cat0$BUGSoutput$summary[2,c(5,3,7)], 2)
OR_IB_Yamagata_2017_2019_cat0

# Bayesian analysis of H3N2 2014 on H3N2 2016, over 5s
IB_2017_Yamagata <- try35[try35$age_2017_cat==1,]$IB_2017_Yamagata
IB_2019_Yamagata <- try35[try35$age_2017_cat==1,]$IB_2019_Yamagata

logistic.model <- function() {
  # SAMPLING DISTRIBUTION
  for (i in 1:N) {
    logit(pi[i]) <- b[1] + b[2]*IB_2017_Yamagata[i];
    IB_2019_Yamagata[i] ~ dbin(pi[i],1);  
  }
  
  # PRIORS ON BETAS
  b[1:Nx] ~ dmnorm(mu[1:Nx], tau[1:Nx,1:Nx]) # multivariate normal prior
  
  # Calculate ORs:
  for (l in 1:Nx) {
    OR[l] <- exp(b[l]);
  }
}

N <- length(IB_2019_Yamagata) 
Nx <- 2 

mu <- rep(0, Nx)
tau <- diag(0.001, Nx)
data.logistic <- list("N", "Nx", "IB_2017_Yamagata", "IB_2019_Yamagata", "mu", "tau")
parameters.logistic <-c("b","OR") # Parameters to keep track of
inits.logistic <- function() {list (b= rep(0, Nx))}

set.seed(114011)
logistic.sim_IB_Yamagata_2017_2019_cat1<-jags(data=data.logistic, inits=inits.logistic, parameters.logistic, n.iter=100000,
                                            n.burn=25000, model.file=logistic.model, n.thin=5, n.chains = 3)

print(logistic.sim_IB_Yamagata_2017_2019_cat1, 2)
OR_IB_Yamagata_2017_2019_cat1 <- round(logistic.sim_IB_Yamagata_2017_2019_cat1$BUGSoutput$summary[2,c(5,3,7)], 2)
OR_IB_Yamagata_2017_2019_cat1

# Yamagata 2014 to 2019


try35 <- data_orig %>%
  filter(at2014==1 & at2019==1 & IB_2015_Yamagata==0 & IB_2016_Yamagata==0 & IB_2017_Yamagata==0 & IB_2018_Yamagata==0) %>%
  dplyr::select(Obs, IB_2014_Yamagata, IB_2019_Yamagata, age_2014_cat)

try35$IB_2014_Yamagata <- as.numeric(try35$IB_2014_Yamagata)
try35$IB_2019_Yamagata <- as.numeric(try35$IB_2019_Yamagata)

# 2x2 table
table(try35[try35$age_2014_cat==0,]$IB_2014_Yamagata, try35[try35$age_2014_cat==0,]$IB_2019_Yamagata) 
table(try35[try35$age_2014_cat==1,]$IB_2014_Yamagata, try35[try35$age_2014_cat==1,]$IB_2019_Yamagata) 


# Bayesian analysis of H3N2 2014 on H3N2 2016, under 5s
IB_2014_Yamagata <- try35[try35$age_2014_cat==0,]$IB_2014_Yamagata
IB_2019_Yamagata <- try35[try35$age_2014_cat==0,]$IB_2019_Yamagata

logistic.model <- function() {
  # SAMPLING DISTRIBUTION
  for (i in 1:N) {
    logit(pi[i]) <- b[1] + b[2]*IB_2014_Yamagata[i];
    IB_2019_Yamagata[i] ~ dbin(pi[i],1);  
  }
  
  # PRIORS ON BETAS
  b[1:Nx] ~ dmnorm(mu[1:Nx], tau[1:Nx,1:Nx]) # multivariate normal prior
  
  # Calculate ORs:
  for (l in 1:Nx) {
    OR[l] <- exp(b[l]);
  }
}

N <- length(IB_2019_Yamagata) 
Nx <- 2 

mu <- rep(0, Nx)
tau <- diag(0.001, Nx)
data.logistic <- list("N", "Nx", "IB_2014_Yamagata", "IB_2019_Yamagata", "mu", "tau")
parameters.logistic <-c("b","OR") # Parameters to keep track of
inits.logistic <- function() {list (b= rep(0, Nx))}

set.seed(114011)
logistic.sim_IB_Yamagata_2014_2019_cat0<-jags(data=data.logistic, inits=inits.logistic, parameters.logistic, n.iter=100000,
                                              n.burn=25000, model.file=logistic.model, n.thin=5, n.chains = 3)

print(logistic.sim_IB_Yamagata_2014_2019_cat0, 2)
OR_IB_Yamagata_2014_2019_cat0 <- round(logistic.sim_IB_Yamagata_2014_2019_cat0$BUGSoutput$summary[2,c(5,3,7)], 2)
OR_IB_Yamagata_2014_2019_cat0

# Bayesian analysis of H3N2 2014 on H3N2 2016, over 5s
IB_2014_Yamagata <- try35[try35$age_2014_cat==1,]$IB_2014_Yamagata
IB_2019_Yamagata <- try35[try35$age_2014_cat==1,]$IB_2019_Yamagata

logistic.model <- function() {
  # SAMPLING DISTRIBUTION
  for (i in 1:N) {
    logit(pi[i]) <- b[1] + b[2]*IB_2014_Yamagata[i];
    IB_2019_Yamagata[i] ~ dbin(pi[i],1);  
  }
  
  # PRIORS ON BETAS
  b[1:Nx] ~ dmnorm(mu[1:Nx], tau[1:Nx,1:Nx]) # multivariate normal prior
  
  # Calculate ORs:
  for (l in 1:Nx) {
    OR[l] <- exp(b[l]);
  }
}

N <- length(IB_2019_Yamagata) 
Nx <- 2 

mu <- rep(0, Nx)
tau <- diag(0.001, Nx)
data.logistic <- list("N", "Nx", "IB_2014_Yamagata", "IB_2019_Yamagata", "mu", "tau")
parameters.logistic <-c("b","OR") # Parameters to keep track of
inits.logistic <- function() {list (b= rep(0, Nx))}

set.seed(114011)
logistic.sim_IB_Yamagata_2014_2019_cat1<-jags(data=data.logistic, inits=inits.logistic, parameters.logistic, n.iter=100000,
                                            n.burn=25000, model.file=logistic.model, n.thin=5, n.chains = 3)

print(logistic.sim_IB_Yamagata_2014_2019_cat1, 2)
OR_IB_Yamagata_2014_2019_cat1 <- round(logistic.sim_IB_Yamagata_2014_2019_cat1$BUGSoutput$summary[2,c(5,3,7)], 2)
OR_IB_Yamagata_2014_2019_cat1

### sensitivity analysis


try35 <- data_orig %>%
  filter(at2014==1 & at2019==1 ) %>%
  dplyr::select(Obs, IB_2014_Yamagata, IB_2019_Yamagata, age_2014_cat)

try35$IB_2014_Yamagata <- as.numeric(try35$IB_2014_Yamagata)
try35$IB_2019_Yamagata <- as.numeric(try35$IB_2019_Yamagata)

# 2x2 table
table(try35[try35$age_2014_cat==0,]$IB_2014_Yamagata, try35[try35$age_2014_cat==0,]$IB_2019_Yamagata) 
table(try35[try35$age_2014_cat==1,]$IB_2014_Yamagata, try35[try35$age_2014_cat==1,]$IB_2019_Yamagata) 


# Bayesian analysis of H3N2 2014 on H3N2 2016, under 5s
IB_2014_Yamagata <- try35[try35$age_2014_cat==0,]$IB_2014_Yamagata
IB_2019_Yamagata <- try35[try35$age_2014_cat==0,]$IB_2019_Yamagata

logistic.model <- function() {
  # SAMPLING DISTRIBUTION
  for (i in 1:N) {
    logit(pi[i]) <- b[1] + b[2]*IB_2014_Yamagata[i];
    IB_2019_Yamagata[i] ~ dbin(pi[i],1);  
  }
  
  # PRIORS ON BETAS
  b[1:Nx] ~ dmnorm(mu[1:Nx], tau[1:Nx,1:Nx]) # multivariate normal prior
  
  # Calculate ORs:
  for (l in 1:Nx) {
    OR[l] <- exp(b[l]);
  }
}

N <- length(IB_2019_Yamagata) 
Nx <- 2 

mu <- rep(0, Nx)
tau <- diag(0.001, Nx)
data.logistic <- list("N", "Nx", "IB_2014_Yamagata", "IB_2019_Yamagata", "mu", "tau")
parameters.logistic <-c("b","OR") # Parameters to keep track of
inits.logistic <- function() {list (b= rep(0, Nx))}

set.seed(114011)
logistic.sim_IB_Yamagata_2014_2019_cat0<-jags(data=data.logistic, inits=inits.logistic, parameters.logistic, n.iter=100000,
                                              n.burn=25000, model.file=logistic.model, n.thin=5, n.chains = 3)

print(logistic.sim_IB_Yamagata_2014_2019_cat0, 2)
OR_IB_Yamagata_2014_2019_cat0 <- round(logistic.sim_IB_Yamagata_2014_2019_cat0$BUGSoutput$summary[2,c(5,3,7)], 2)
OR_IB_Yamagata_2014_2019_cat0

# Bayesian analysis of H3N2 2014 on H3N2 2016, over 5s
IB_2014_Yamagata <- try35[try35$age_2014_cat==1,]$IB_2014_Yamagata
IB_2019_Yamagata <- try35[try35$age_2014_cat==1,]$IB_2019_Yamagata

logistic.model <- function() {
  # SAMPLING DISTRIBUTION
  for (i in 1:N) {
    logit(pi[i]) <- b[1] + b[2]*IB_2014_Yamagata[i];
    IB_2019_Yamagata[i] ~ dbin(pi[i],1);  
  }
  
  # PRIORS ON BETAS
  b[1:Nx] ~ dmnorm(mu[1:Nx], tau[1:Nx,1:Nx]) # multivariate normal prior
  
  # Calculate ORs:
  for (l in 1:Nx) {
    OR[l] <- exp(b[l]);
  }
}

N <- length(IB_2019_Yamagata) 
Nx <- 2 

mu <- rep(0, Nx)
tau <- diag(0.001, Nx)
data.logistic <- list("N", "Nx", "IB_2014_Yamagata", "IB_2019_Yamagata", "mu", "tau")
parameters.logistic <-c("b","OR") # Parameters to keep track of
inits.logistic <- function() {list (b= rep(0, Nx))}

set.seed(114011)
logistic.sim_IB_Yamagata_2014_2019_cat1<-jags(data=data.logistic, inits=inits.logistic, parameters.logistic, n.iter=100000,
                                              n.burn=25000, model.file=logistic.model, n.thin=5, n.chains = 3)

print(logistic.sim_IB_Yamagata_2014_2019_cat1, 2)
OR_IB_Yamagata_2014_2019_cat1 <- round(logistic.sim_IB_Yamagata_2014_2019_cat1$BUGSoutput$summary[2,c(5,3,7)], 2)
OR_IB_Yamagata_2014_2019_cat1

# Yamagata 2014 to 2017


try35 <- data_orig %>%
  filter(at2014==1 & at2017==1 & IB_2015_Yamagata==0 & IB_2016_Yamagata==0 ) %>%
  dplyr::select(Obs, IB_2014_Yamagata, IB_2017_Yamagata, age_2014_cat)

try35$IB_2014_Yamagata <- as.numeric(try35$IB_2014_Yamagata)
try35$IB_2017_Yamagata <- as.numeric(try35$IB_2017_Yamagata)

# 2x2 table
table(try35[try35$age_2014_cat==0,]$IB_2014_Yamagata, try35[try35$age_2014_cat==0,]$IB_2017_Yamagata) 
table(try35[try35$age_2014_cat==1,]$IB_2014_Yamagata, try35[try35$age_2014_cat==1,]$IB_2017_Yamagata) 


# Bayesian analysis of H3N2 2014 on H3N2 2016, under 5s
IB_2014_Yamagata <- try35[try35$age_2014_cat==0,]$IB_2014_Yamagata
IB_2017_Yamagata <- try35[try35$age_2014_cat==0,]$IB_2017_Yamagata

logistic.model <- function() {
  # SAMPLING DISTRIBUTION
  for (i in 1:N) {
    logit(pi[i]) <- b[1] + b[2]*IB_2014_Yamagata[i];
    IB_2017_Yamagata[i] ~ dbin(pi[i],1);  
  }
  
  # PRIORS ON BETAS
  b[1:Nx] ~ dmnorm(mu[1:Nx], tau[1:Nx,1:Nx]) # multivariate normal prior
  
  # Calculate ORs:
  for (l in 1:Nx) {
    OR[l] <- exp(b[l]);
  }
}

N <- length(IB_2017_Yamagata) 
Nx <- 2 

mu <- rep(0, Nx)
tau <- diag(0.001, Nx)
data.logistic <- list("N", "Nx", "IB_2014_Yamagata", "IB_2017_Yamagata", "mu", "tau")
parameters.logistic <-c("b","OR") # Parameters to keep track of
inits.logistic <- function() {list (b= rep(0, Nx))}

set.seed(114011)
logistic.sim_IB_Yamagata_2014_2017_cat0<-jags(data=data.logistic, inits=inits.logistic, parameters.logistic, n.iter=100000,
                                              n.burn=25000, model.file=logistic.model, n.thin=5, n.chains = 3)

print(logistic.sim_IB_Yamagata_2014_2017_cat0, 2)
OR_IB_Yamagata_2014_2017_cat0 <- round(logistic.sim_IB_Yamagata_2014_2017_cat0$BUGSoutput$summary[2,c(5,3,7)], 2)
OR_IB_Yamagata_2014_2017_cat0

# Bayesian analysis of H3N2 2014 on H3N2 2016, over 5s
IB_2014_Yamagata <- try35[try35$age_2014_cat==1,]$IB_2014_Yamagata
IB_2017_Yamagata <- try35[try35$age_2014_cat==1,]$IB_2017_Yamagata

logistic.model <- function() {
  # SAMPLING DISTRIBUTION
  for (i in 1:N) {
    logit(pi[i]) <- b[1] + b[2]*IB_2014_Yamagata[i];
    IB_2017_Yamagata[i] ~ dbin(pi[i],1);  
  }
  
  # PRIORS ON BETAS
  b[1:Nx] ~ dmnorm(mu[1:Nx], tau[1:Nx,1:Nx]) # multivariate normal prior
  
  # Calculate ORs:
  for (l in 1:Nx) {
    OR[l] <- exp(b[l]);
  }
}

N <- length(IB_2017_Yamagata) 
Nx <- 2 

mu <- rep(0, Nx)
tau <- diag(0.001, Nx)
data.logistic <- list("N", "Nx", "IB_2014_Yamagata", "IB_2017_Yamagata", "mu", "tau")
parameters.logistic <-c("b","OR") # Parameters to keep track of
inits.logistic <- function() {list (b= rep(0, Nx))}

set.seed(114011)
logistic.sim_IB_Yamagata_2014_2017_cat1<-jags(data=data.logistic, inits=inits.logistic, parameters.logistic, n.iter=100000,
                                              n.burn=25000, model.file=logistic.model, n.thin=5, n.chains = 3)

print(logistic.sim_IB_Yamagata_2014_2017_cat1, 2)
OR_IB_Yamagata_2014_2017_cat1 <- round(logistic.sim_IB_Yamagata_2014_2017_cat1$BUGSoutput$summary[2,c(5,3,7)], 2)
OR_IB_Yamagata_2014_2017_cat1

### sensitivity analysis


try35 <- data_orig %>%
  filter(at2014==1 & at2017==1 ) %>%
  dplyr::select(Obs, IB_2014_Yamagata, IB_2017_Yamagata, age_2014_cat)

try35$IB_2014_Yamagata <- as.numeric(try35$IB_2014_Yamagata)
try35$IB_2017_Yamagata <- as.numeric(try35$IB_2017_Yamagata)

# 2x2 table
table(try35[try35$age_2014_cat==0,]$IB_2014_Yamagata, try35[try35$age_2014_cat==0,]$IB_2017_Yamagata) 
table(try35[try35$age_2014_cat==1,]$IB_2014_Yamagata, try35[try35$age_2014_cat==1,]$IB_2017_Yamagata) 


# Bayesian analysis of H3N2 2014 on H3N2 2016, under 5s
IB_2014_Yamagata <- try35[try35$age_2014_cat==0,]$IB_2014_Yamagata
IB_2017_Yamagata <- try35[try35$age_2014_cat==0,]$IB_2017_Yamagata

logistic.model <- function() {
  # SAMPLING DISTRIBUTION
  for (i in 1:N) {
    logit(pi[i]) <- b[1] + b[2]*IB_2014_Yamagata[i];
    IB_2017_Yamagata[i] ~ dbin(pi[i],1);  
  }
  
  # PRIORS ON BETAS
  b[1:Nx] ~ dmnorm(mu[1:Nx], tau[1:Nx,1:Nx]) # multivariate normal prior
  
  # Calculate ORs:
  for (l in 1:Nx) {
    OR[l] <- exp(b[l]);
  }
}

N <- length(IB_2017_Yamagata) 
Nx <- 2 

mu <- rep(0, Nx)
tau <- diag(0.001, Nx)
data.logistic <- list("N", "Nx", "IB_2014_Yamagata", "IB_2017_Yamagata", "mu", "tau")
parameters.logistic <-c("b","OR") # Parameters to keep track of
inits.logistic <- function() {list (b= rep(0, Nx))}

set.seed(114011)
logistic.sim_IB_Yamagata_2014_2017_cat0<-jags(data=data.logistic, inits=inits.logistic, parameters.logistic, n.iter=100000,
                                              n.burn=25000, model.file=logistic.model, n.thin=5, n.chains = 3)

print(logistic.sim_IB_Yamagata_2014_2017_cat0, 2)
OR_IB_Yamagata_2014_2017_cat0 <- round(logistic.sim_IB_Yamagata_2014_2017_cat0$BUGSoutput$summary[2,c(5,3,7)], 2)
OR_IB_Yamagata_2014_2017_cat0

# Bayesian analysis of H3N2 2014 on H3N2 2016, over 5s
IB_2014_Yamagata <- try35[try35$age_2014_cat==1,]$IB_2014_Yamagata
IB_2017_Yamagata <- try35[try35$age_2014_cat==1,]$IB_2017_Yamagata

logistic.model <- function() {
  # SAMPLING DISTRIBUTION
  for (i in 1:N) {
    logit(pi[i]) <- b[1] + b[2]*IB_2014_Yamagata[i];
    IB_2017_Yamagata[i] ~ dbin(pi[i],1);  
  }
  
  # PRIORS ON BETAS
  b[1:Nx] ~ dmnorm(mu[1:Nx], tau[1:Nx,1:Nx]) # multivariate normal prior
  
  # Calculate ORs:
  for (l in 1:Nx) {
    OR[l] <- exp(b[l]);
  }
}

N <- length(IB_2017_Yamagata) 
Nx <- 2 

mu <- rep(0, Nx)
tau <- diag(0.001, Nx)
data.logistic <- list("N", "Nx", "IB_2014_Yamagata", "IB_2017_Yamagata", "mu", "tau")
parameters.logistic <-c("b","OR") # Parameters to keep track of
inits.logistic <- function() {list (b= rep(0, Nx))}

set.seed(114011)
logistic.sim_IB_Yamagata_2014_2017_cat1<-jags(data=data.logistic, inits=inits.logistic, parameters.logistic, n.iter=100000,
                                              n.burn=25000, model.file=logistic.model, n.thin=5, n.chains = 3)

print(logistic.sim_IB_Yamagata_2014_2017_cat1, 2)
OR_IB_Yamagata_2014_2017_cat1 <- round(logistic.sim_IB_Yamagata_2014_2017_cat1$BUGSoutput$summary[2,c(5,3,7)], 2)
OR_IB_Yamagata_2014_2017_cat1

# VICTORIA ANALYSES


# IB 2012 to 2016 - Victoria
try52 <- data_orig %>%
  filter(at2012==1 & at2016==1) %>%
  dplyr::select(Obs, IB_2012_Victoria, IB_2016_Victoria)
try52$IB_2012_Victoria <- as.numeric(try52$IB_2012_Victoria)
try52$IB_2016_Victoria <- as.numeric(try52$IB_2016_Victoria)

table(try52$IB_2012_Victoria, try52$IB_2016_Victoria)

# Bayesian analyses
# Putting exposure and outcome variables in the environment 
IB_2012_Victoria <- try52$IB_2012_Victoria
IB_2016_Victoria <- try52$IB_2016_Victoria

### Bayesian logistic regression
logistic.model <- function() {
  # SAMPLING DISTRIBUTION
  for (i in 1:N) {
    logit(pi[i]) <- b[1] + b[2]*IB_2012_Victoria[i];
    IB_2016_Victoria[i] ~ dbin(pi[i],1);  
  }
  
  # PRIORS ON BETAS
  b[1:Nx] ~ dmnorm(mu[1:Nx], tau[1:Nx,1:Nx]) # multivariate normal prior
  
  # Calculate ORs:
  for (l in 1:Nx) {
    OR[l] <- exp(b[l]);
  }
}

N <- length(IB_2016_Victoria) # Number of observations
Nx <- 2 # Number of variables: intercept + exposure 

# Data, parameter list and starting values
mu <- rep(0, Nx)
tau <- diag(0.001, Nx)
data.logistic <- list("N", "Nx", "IB_2012_Victoria", "IB_2016_Victoria", "mu", "tau")
parameters.logistic <-c("b","OR") # Parameters to keep track of
inits.logistic <- function() {list (b= rep(0, Nx))}

set.seed(114011)
logistic.sim<-jags(data=data.logistic, inits=inits.logistic, parameters.logistic, n.iter=100000,
                   n.burn=25000, model.file=logistic.model, n.thin=5, n.chains = 3)

print(logistic.sim, 2)
OR_IB_2012_2016_Victoria_Bayes_SA <- round(logistic.sim$BUGSoutput$summary[2,c(5,3,7)], 2)
OR_IB_2012_2016_Victoria_Bayes_SA



# IB 2012 to 2017 - Victoria
try52 <- data_orig %>%
  filter(at2012==1 & at2017==1 & IB_2016_Victoria==0) %>%
  dplyr::select(Obs, IB_2012_Victoria, IB_2017_Victoria)
try52$IB_2012_Victoria <- as.numeric(try52$IB_2012_Victoria)
try52$IB_2017_Victoria <- as.numeric(try52$IB_2017_Victoria)

table(try52$IB_2012_Victoria, try52$IB_2017_Victoria)

# Bayesian analyses
# Putting exposure and outcome variables in the environment 
IB_2012_Victoria <- try52$IB_2012_Victoria
IB_2017_Victoria <- try52$IB_2017_Victoria

### Bayesian logistic regression
logistic.model <- function() {
  # SAMPLING DISTRIBUTION
  for (i in 1:N) {
    logit(pi[i]) <- b[1] + b[2]*IB_2012_Victoria[i];
    IB_2017_Victoria[i] ~ dbin(pi[i],1);  
  }
  
  # PRIORS ON BETAS
  b[1:Nx] ~ dmnorm(mu[1:Nx], tau[1:Nx,1:Nx]) # multivariate normal prior
  
  # Calculate ORs:
  for (l in 1:Nx) {
    OR[l] <- exp(b[l]);
  }
}

N <- length(IB_2017_Victoria) # Number of observations
Nx <- 2 # Number of variables: intercept + exposure 

# Data, parameter list and starting values
mu <- rep(0, Nx)
tau <- diag(0.001, Nx)
data.logistic <- list("N", "Nx", "IB_2012_Victoria", "IB_2017_Victoria", "mu", "tau")
parameters.logistic <-c("b","OR") # Parameters to keep track of
inits.logistic <- function() {list (b= rep(0, Nx))}

set.seed(114011)
logistic.sim<-jags(data=data.logistic, inits=inits.logistic, parameters.logistic, n.iter=100000,
                   n.burn=25000, model.file=logistic.model, n.thin=5, n.chains = 3)

print(logistic.sim, 2)
OR_IB_2012_2017_Victoria_Bayes_SA <- round(logistic.sim$BUGSoutput$summary[2,c(5,3,7)], 2)
OR_IB_2012_2017_Victoria_Bayes_SA

# sensitivity analysis

try52 <- data_orig %>%
  filter(at2012==1 & at2017==1 ) %>%
  dplyr::select(Obs, IB_2012_Victoria, IB_2017_Victoria)
try52$IB_2012_Victoria <- as.numeric(try52$IB_2012_Victoria)
try52$IB_2017_Victoria <- as.numeric(try52$IB_2017_Victoria)

table(try52$IB_2012_Victoria, try52$IB_2017_Victoria)

# Bayesian analyses
# Putting exposure and outcome variables in the environment 
IB_2012_Victoria <- try52$IB_2012_Victoria
IB_2017_Victoria <- try52$IB_2017_Victoria

### Bayesian logistic regression
logistic.model <- function() {
  # SAMPLING DISTRIBUTION
  for (i in 1:N) {
    logit(pi[i]) <- b[1] + b[2]*IB_2012_Victoria[i];
    IB_2017_Victoria[i] ~ dbin(pi[i],1);  
  }
  
  # PRIORS ON BETAS
  b[1:Nx] ~ dmnorm(mu[1:Nx], tau[1:Nx,1:Nx]) # multivariate normal prior
  
  # Calculate ORs:
  for (l in 1:Nx) {
    OR[l] <- exp(b[l]);
  }
}

N <- length(IB_2017_Victoria) # Number of observations
Nx <- 2 # Number of variables: intercept + exposure 

# Data, parameter list and starting values
mu <- rep(0, Nx)
tau <- diag(0.001, Nx)
data.logistic <- list("N", "Nx", "IB_2012_Victoria", "IB_2017_Victoria", "mu", "tau")
parameters.logistic <-c("b","OR") # Parameters to keep track of
inits.logistic <- function() {list (b= rep(0, Nx))}

set.seed(114011)
logistic.sim<-jags(data=data.logistic, inits=inits.logistic, parameters.logistic, n.iter=100000,
                   n.burn=25000, model.file=logistic.model, n.thin=5, n.chains = 3)

print(logistic.sim, 2)
OR_IB_2012_2017_Victoria_Bayes_SA <- round(logistic.sim$BUGSoutput$summary[2,c(5,3,7)], 2)
OR_IB_2012_2017_Victoria_Bayes_SA


# IB 2012 to 2018 - Victoria
try52 <- data_orig %>%
  filter(at2012==1 & at2018==1 & IB_2016_Victoria==0 & IB_2017_Victoria==0) %>%
  dplyr::select(Obs, IB_2012_Victoria, IB_2018_Victoria)
try52$IB_2012_Victoria <- as.numeric(try52$IB_2012_Victoria)
try52$IB_2018_Victoria <- as.numeric(try52$IB_2018_Victoria)

table(try52$IB_2012_Victoria, try52$IB_2018_Victoria)

# Bayesian analyses
# Putting exposure and outcome variables in the environment 
IB_2012_Victoria <- try52$IB_2012_Victoria
IB_2018_Victoria <- try52$IB_2018_Victoria

### Bayesian logistic regression
logistic.model <- function() {
  # SAMPLING DISTRIBUTION
  for (i in 1:N) {
    logit(pi[i]) <- b[1] + b[2]*IB_2012_Victoria[i];
    IB_2018_Victoria[i] ~ dbin(pi[i],1);  
  }
  
  # PRIORS ON BETAS
  b[1:Nx] ~ dmnorm(mu[1:Nx], tau[1:Nx,1:Nx]) # multivariate normal prior
  
  # Calculate ORs:
  for (l in 1:Nx) {
    OR[l] <- exp(b[l]);
  }
}

N <- length(IB_2018_Victoria) # Number of observations
Nx <- 2 # Number of variables: intercept + exposure 

# Data, parameter list and starting values
mu <- rep(0, Nx)
tau <- diag(0.001, Nx)
data.logistic <- list("N", "Nx", "IB_2012_Victoria", "IB_2018_Victoria", "mu", "tau")
parameters.logistic <-c("b","OR") # Parameters to keep track of
inits.logistic <- function() {list (b= rep(0, Nx))}

set.seed(114011)
logistic.sim<-jags(data=data.logistic, inits=inits.logistic, parameters.logistic, n.iter=100000,
                   n.burn=25000, model.file=logistic.model, n.thin=5, n.chains = 3)

print(logistic.sim, 2)
OR_IB_2012_2018_Victoria_Bayes_SA <- round(logistic.sim$BUGSoutput$summary[2,c(5,3,7)], 2)
OR_IB_2012_2018_Victoria_Bayes_SA

# sensitivity analysis

try52 <- data_orig %>%
  filter(at2012==1 & at2018==1) %>%
  dplyr::select(Obs, IB_2012_Victoria, IB_2018_Victoria)
try52$IB_2012_Victoria <- as.numeric(try52$IB_2012_Victoria)
try52$IB_2018_Victoria <- as.numeric(try52$IB_2018_Victoria)

table(try52$IB_2012_Victoria, try52$IB_2018_Victoria)

# Bayesian analyses
# Putting exposure and outcome variables in the environment 
IB_2012_Victoria <- try52$IB_2012_Victoria
IB_2018_Victoria <- try52$IB_2018_Victoria

### Bayesian logistic regression
logistic.model <- function() {
  # SAMPLING DISTRIBUTION
  for (i in 1:N) {
    logit(pi[i]) <- b[1] + b[2]*IB_2012_Victoria[i];
    IB_2018_Victoria[i] ~ dbin(pi[i],1);  
  }
  
  # PRIORS ON BETAS
  b[1:Nx] ~ dmnorm(mu[1:Nx], tau[1:Nx,1:Nx]) # multivariate normal prior
  
  # Calculate ORs:
  for (l in 1:Nx) {
    OR[l] <- exp(b[l]);
  }
}

N <- length(IB_2018_Victoria) # Number of observations
Nx <- 2 # Number of variables: intercept + exposure 

# Data, parameter list and starting values
mu <- rep(0, Nx)
tau <- diag(0.001, Nx)
data.logistic <- list("N", "Nx", "IB_2012_Victoria", "IB_2018_Victoria", "mu", "tau")
parameters.logistic <-c("b","OR") # Parameters to keep track of
inits.logistic <- function() {list (b= rep(0, Nx))}

set.seed(114011)
logistic.sim<-jags(data=data.logistic, inits=inits.logistic, parameters.logistic, n.iter=100000,
                   n.burn=25000, model.file=logistic.model, n.thin=5, n.chains = 3)

print(logistic.sim, 2)
OR_IB_2012_2018_Victoria_Bayes_SA <- round(logistic.sim$BUGSoutput$summary[2,c(5,3,7)], 2)
OR_IB_2012_2018_Victoria_Bayes_SA



# IB 2012 to 2019 - Victoria
try52 <- data_orig %>%
  filter(at2012==1 & at2019==1 & IB_2016_Victoria==0 & IB_2017_Victoria==0 & IB_2018_Victoria==0) %>%
  dplyr::select(Obs, IB_2012_Victoria, IB_2019_Victoria)
try52$IB_2012_Victoria <- as.numeric(try52$IB_2012_Victoria)
try52$IB_2019_Victoria <- as.numeric(try52$IB_2019_Victoria)

table(try52$IB_2012_Victoria, try52$IB_2019_Victoria)

# Bayesian analyses
# Putting exposure and outcome variables in the environment 
IB_2012_Victoria <- try52$IB_2012_Victoria
IB_2019_Victoria <- try52$IB_2019_Victoria

### Bayesian logistic regression
logistic.model <- function() {
  # SAMPLING DISTRIBUTION
  for (i in 1:N) {
    logit(pi[i]) <- b[1] + b[2]*IB_2012_Victoria[i];
    IB_2019_Victoria[i] ~ dbin(pi[i],1);  
  }
  
  # PRIORS ON BETAS
  b[1:Nx] ~ dmnorm(mu[1:Nx], tau[1:Nx,1:Nx]) # multivariate normal prior
  
  # Calculate ORs:
  for (l in 1:Nx) {
    OR[l] <- exp(b[l]);
  }
}

N <- length(IB_2019_Victoria) # Number of observations
Nx <- 2 # Number of variables: intercept + exposure 

# Data, parameter list and starting values
mu <- rep(0, Nx)
tau <- diag(0.001, Nx)
data.logistic <- list("N", "Nx", "IB_2012_Victoria", "IB_2019_Victoria", "mu", "tau")
parameters.logistic <-c("b","OR") # Parameters to keep track of
inits.logistic <- function() {list (b= rep(0, Nx))}

set.seed(114011)
logistic.sim<-jags(data=data.logistic, inits=inits.logistic, parameters.logistic, n.iter=100000,
                   n.burn=25000, model.file=logistic.model, n.thin=5, n.chains = 3)

print(logistic.sim, 2)
OR_IB_2012_2019_Victoria_Bayes_SA <- round(logistic.sim$BUGSoutput$summary[2,c(5,3,7)], 2)
OR_IB_2012_2019_Victoria_Bayes_SA

# sensitivity analysis

try52 <- data_orig %>%
  filter(at2012==1 & at2019==1) %>%
  dplyr::select(Obs, IB_2012_Victoria, IB_2019_Victoria)
try52$IB_2012_Victoria <- as.numeric(try52$IB_2012_Victoria)
try52$IB_2019_Victoria <- as.numeric(try52$IB_2019_Victoria)

table(try52$IB_2012_Victoria, try52$IB_2019_Victoria)

# Bayesian analyses
# Putting exposure and outcome variables in the environment 
IB_2012_Victoria <- try52$IB_2012_Victoria
IB_2019_Victoria <- try52$IB_2019_Victoria

### Bayesian logistic regression
logistic.model <- function() {
  # SAMPLING DISTRIBUTION
  for (i in 1:N) {
    logit(pi[i]) <- b[1] + b[2]*IB_2012_Victoria[i];
    IB_2019_Victoria[i] ~ dbin(pi[i],1);  
  }
  
  # PRIORS ON BETAS
  b[1:Nx] ~ dmnorm(mu[1:Nx], tau[1:Nx,1:Nx]) # multivariate normal prior
  
  # Calculate ORs:
  for (l in 1:Nx) {
    OR[l] <- exp(b[l]);
  }
}

N <- length(IB_2019_Victoria) # Number of observations
Nx <- 2 # Number of variables: intercept + exposure 

# Data, parameter list and starting values
mu <- rep(0, Nx)
tau <- diag(0.001, Nx)
data.logistic <- list("N", "Nx", "IB_2012_Victoria", "IB_2019_Victoria", "mu", "tau")
parameters.logistic <-c("b","OR") # Parameters to keep track of
inits.logistic <- function() {list (b= rep(0, Nx))}

set.seed(114011)
logistic.sim<-jags(data=data.logistic, inits=inits.logistic, parameters.logistic, n.iter=100000,
                   n.burn=25000, model.file=logistic.model, n.thin=5, n.chains = 3)

print(logistic.sim, 2)
OR_IB_2012_2019_Victoria_Bayes_SA <- round(logistic.sim$BUGSoutput$summary[2,c(5,3,7)], 2)
OR_IB_2012_2019_Victoria_Bayes_SA



# IB 2016 to 2017 - Victoria
try52 <- data_orig %>%
  filter(at2016==1 & at2017==1 ) %>%
  dplyr::select(Obs, IB_2016_Victoria, IB_2017_Victoria)
try52$IB_2016_Victoria <- as.numeric(try52$IB_2016_Victoria)
try52$IB_2017_Victoria <- as.numeric(try52$IB_2017_Victoria)

table(try52$IB_2016_Victoria, try52$IB_2017_Victoria)

# Bayesian analyses
# Putting exposure and outcome variables in the environment 
IB_2016_Victoria <- try52$IB_2016_Victoria
IB_2017_Victoria <- try52$IB_2017_Victoria

### Bayesian logistic regression
logistic.model <- function() {
  # SAMPLING DISTRIBUTION
  for (i in 1:N) {
    logit(pi[i]) <- b[1] + b[2]*IB_2016_Victoria[i];
    IB_2017_Victoria[i] ~ dbin(pi[i],1);  
  }
  
  # PRIORS ON BETAS
  b[1:Nx] ~ dmnorm(mu[1:Nx], tau[1:Nx,1:Nx]) # multivariate normal prior
  
  # Calculate ORs:
  for (l in 1:Nx) {
    OR[l] <- exp(b[l]);
  }
}

N <- length(IB_2017_Victoria) # Number of observations
Nx <- 2 # Number of variables: intercept + exposure 

# Data, parameter list and starting values
mu <- rep(0, Nx)
tau <- diag(0.001, Nx)
data.logistic <- list("N", "Nx", "IB_2016_Victoria", "IB_2017_Victoria", "mu", "tau")
parameters.logistic <-c("b","OR") # Parameters to keep track of
inits.logistic <- function() {list (b= rep(0, Nx))}

set.seed(114011)
logistic.sim<-jags(data=data.logistic, inits=inits.logistic, parameters.logistic, n.iter=100000,
                   n.burn=25000, model.file=logistic.model, n.thin=5, n.chains = 3)

print(logistic.sim, 2)
OR_IB_2016_2017_Victoria_Bayes_SA <- round(logistic.sim$BUGSoutput$summary[2,c(5,3,7)], 2)
OR_IB_2016_2017_Victoria_Bayes_SA




# IB 2016 to 2018 - Victoria
try52 <- data_orig %>%
  filter(at2016==1 & at2018==1 & IB_2017_Victoria==0) %>%
  dplyr::select(Obs, IB_2016_Victoria, IB_2018_Victoria)
try52$IB_2016_Victoria <- as.numeric(try52$IB_2016_Victoria)
try52$IB_2018_Victoria <- as.numeric(try52$IB_2018_Victoria)

table(try52$IB_2016_Victoria, try52$IB_2018_Victoria)

# Bayesian analyses
# Putting exposure and outcome variables in the environment 
IB_2016_Victoria <- try52$IB_2016_Victoria
IB_2018_Victoria <- try52$IB_2018_Victoria

### Bayesian logistic regression
logistic.model <- function() {
  # SAMPLING DISTRIBUTION
  for (i in 1:N) {
    logit(pi[i]) <- b[1] + b[2]*IB_2016_Victoria[i];
    IB_2018_Victoria[i] ~ dbin(pi[i],1);  
  }
  
  # PRIORS ON BETAS
  b[1:Nx] ~ dmnorm(mu[1:Nx], tau[1:Nx,1:Nx]) # multivariate normal prior
  
  # Calculate ORs:
  for (l in 1:Nx) {
    OR[l] <- exp(b[l]);
  }
}

N <- length(IB_2018_Victoria) # Number of observations
Nx <- 2 # Number of variables: intercept + exposure 

# Data, parameter list and starting values
mu <- rep(0, Nx)
tau <- diag(0.001, Nx)
data.logistic <- list("N", "Nx", "IB_2016_Victoria", "IB_2018_Victoria", "mu", "tau")
parameters.logistic <-c("b","OR") # Parameters to keep track of
inits.logistic <- function() {list (b= rep(0, Nx))}

set.seed(114011)
logistic.sim<-jags(data=data.logistic, inits=inits.logistic, parameters.logistic, n.iter=100000,
                   n.burn=25000, model.file=logistic.model, n.thin=5, n.chains = 3)

print(logistic.sim, 2)
OR_IB_2016_2018_Victoria_Bayes_SA <- round(logistic.sim$BUGSoutput$summary[2,c(5,3,7)], 2)
OR_IB_2016_2018_Victoria_Bayes_SA

# sensitivity analysis

try52 <- data_orig %>%
  filter(at2016==1 & at2018==1 ) %>%
  dplyr::select(Obs, IB_2016_Victoria, IB_2018_Victoria)
try52$IB_2016_Victoria <- as.numeric(try52$IB_2016_Victoria)
try52$IB_2018_Victoria <- as.numeric(try52$IB_2018_Victoria)

table(try52$IB_2016_Victoria, try52$IB_2018_Victoria)

# Bayesian analyses
# Putting exposure and outcome variables in the environment 
IB_2016_Victoria <- try52$IB_2016_Victoria
IB_2018_Victoria <- try52$IB_2018_Victoria

### Bayesian logistic regression
logistic.model <- function() {
  # SAMPLING DISTRIBUTION
  for (i in 1:N) {
    logit(pi[i]) <- b[1] + b[2]*IB_2016_Victoria[i];
    IB_2018_Victoria[i] ~ dbin(pi[i],1);  
  }
  
  # PRIORS ON BETAS
  b[1:Nx] ~ dmnorm(mu[1:Nx], tau[1:Nx,1:Nx]) # multivariate normal prior
  
  # Calculate ORs:
  for (l in 1:Nx) {
    OR[l] <- exp(b[l]);
  }
}

N <- length(IB_2018_Victoria) # Number of observations
Nx <- 2 # Number of variables: intercept + exposure 

# Data, parameter list and starting values
mu <- rep(0, Nx)
tau <- diag(0.001, Nx)
data.logistic <- list("N", "Nx", "IB_2016_Victoria", "IB_2018_Victoria", "mu", "tau")
parameters.logistic <-c("b","OR") # Parameters to keep track of
inits.logistic <- function() {list (b= rep(0, Nx))}

set.seed(114011)
logistic.sim<-jags(data=data.logistic, inits=inits.logistic, parameters.logistic, n.iter=100000,
                   n.burn=25000, model.file=logistic.model, n.thin=5, n.chains = 3)

print(logistic.sim, 2)
OR_IB_2016_2018_Victoria_Bayes_SA <- round(logistic.sim$BUGSoutput$summary[2,c(5,3,7)], 2)
OR_IB_2016_2018_Victoria_Bayes_SA




# IB 2016 to 2019 - Victoria
try52 <- data_orig %>%
  filter(at2016==1 & at2019==1 & IB_2017_Victoria==0 & IB_2018_Victoria==0) %>%
  dplyr::select(Obs, IB_2016_Victoria, IB_2019_Victoria)
try52$IB_2016_Victoria <- as.numeric(try52$IB_2016_Victoria)
try52$IB_2019_Victoria <- as.numeric(try52$IB_2019_Victoria)

table(try52$IB_2016_Victoria, try52$IB_2019_Victoria)

# Bayesian analyses
# Putting exposure and outcome variables in the environment 
IB_2016_Victoria <- try52$IB_2016_Victoria
IB_2019_Victoria <- try52$IB_2019_Victoria

### Bayesian logistic regression
logistic.model <- function() {
  # SAMPLING DISTRIBUTION
  for (i in 1:N) {
    logit(pi[i]) <- b[1] + b[2]*IB_2016_Victoria[i];
    IB_2019_Victoria[i] ~ dbin(pi[i],1);  
  }
  
  # PRIORS ON BETAS
  b[1:Nx] ~ dmnorm(mu[1:Nx], tau[1:Nx,1:Nx]) # multivariate normal prior
  
  # Calculate ORs:
  for (l in 1:Nx) {
    OR[l] <- exp(b[l]);
  }
}

N <- length(IB_2019_Victoria) # Number of observations
Nx <- 2 # Number of variables: intercept + exposure 

# Data, parameter list and starting values
mu <- rep(0, Nx)
tau <- diag(0.001, Nx)
data.logistic <- list("N", "Nx", "IB_2016_Victoria", "IB_2019_Victoria", "mu", "tau")
parameters.logistic <-c("b","OR") # Parameters to keep track of
inits.logistic <- function() {list (b= rep(0, Nx))}

set.seed(114011)
logistic.sim<-jags(data=data.logistic, inits=inits.logistic, parameters.logistic, n.iter=100000,
                   n.burn=25000, model.file=logistic.model, n.thin=5, n.chains = 3)

print(logistic.sim, 2)
OR_IB_2016_2019_Victoria_Bayes_SA <- round(logistic.sim$BUGSoutput$summary[2,c(5,3,7)], 2)
OR_IB_2016_2019_Victoria_Bayes_SA

# sensitivity analysis

try52 <- data_orig %>%
  filter(at2016==1 & at2019==1 ) %>%
  dplyr::select(Obs, IB_2016_Victoria, IB_2019_Victoria)
try52$IB_2016_Victoria <- as.numeric(try52$IB_2016_Victoria)
try52$IB_2019_Victoria <- as.numeric(try52$IB_2019_Victoria)

table(try52$IB_2016_Victoria, try52$IB_2019_Victoria)

# Bayesian analyses
# Putting exposure and outcome variables in the environment 
IB_2016_Victoria <- try52$IB_2016_Victoria
IB_2019_Victoria <- try52$IB_2019_Victoria

### Bayesian logistic regression
logistic.model <- function() {
  # SAMPLING DISTRIBUTION
  for (i in 1:N) {
    logit(pi[i]) <- b[1] + b[2]*IB_2016_Victoria[i];
    IB_2019_Victoria[i] ~ dbin(pi[i],1);  
  }
  
  # PRIORS ON BETAS
  b[1:Nx] ~ dmnorm(mu[1:Nx], tau[1:Nx,1:Nx]) # multivariate normal prior
  
  # Calculate ORs:
  for (l in 1:Nx) {
    OR[l] <- exp(b[l]);
  }
}

N <- length(IB_2019_Victoria) # Number of observations
Nx <- 2 # Number of variables: intercept + exposure 

# Data, parameter list and starting values
mu <- rep(0, Nx)
tau <- diag(0.001, Nx)
data.logistic <- list("N", "Nx", "IB_2016_Victoria", "IB_2019_Victoria", "mu", "tau")
parameters.logistic <-c("b","OR") # Parameters to keep track of
inits.logistic <- function() {list (b= rep(0, Nx))}

set.seed(114011)
logistic.sim<-jags(data=data.logistic, inits=inits.logistic, parameters.logistic, n.iter=100000,
                   n.burn=25000, model.file=logistic.model, n.thin=5, n.chains = 3)

print(logistic.sim, 2)
OR_IB_2016_2019_Victoria_Bayes_SA <- round(logistic.sim$BUGSoutput$summary[2,c(5,3,7)], 2)
OR_IB_2016_2019_Victoria_Bayes_SA


# IB 2017 to 2018 - Victoria
try52 <- data_orig %>%
  filter(at2017==1 & at2018==1) %>%
  dplyr::select(Obs, IB_2017_Victoria, IB_2018_Victoria)
try52$IB_2017_Victoria <- as.numeric(try52$IB_2017_Victoria)
try52$IB_2018_Victoria <- as.numeric(try52$IB_2018_Victoria)

table(try52$IB_2017_Victoria, try52$IB_2018_Victoria)

# Bayesian analyses
# Putting exposure and outcome variables in the environment 
IB_2017_Victoria <- try52$IB_2017_Victoria
IB_2018_Victoria <- try52$IB_2018_Victoria

### Bayesian logistic regression
logistic.model <- function() {
  # SAMPLING DISTRIBUTION
  for (i in 1:N) {
    logit(pi[i]) <- b[1] + b[2]*IB_2017_Victoria[i];
    IB_2018_Victoria[i] ~ dbin(pi[i],1);  
  }
  
  # PRIORS ON BETAS
  b[1:Nx] ~ dmnorm(mu[1:Nx], tau[1:Nx,1:Nx]) # multivariate normal prior
  
  # Calculate ORs:
  for (l in 1:Nx) {
    OR[l] <- exp(b[l]);
  }
}

N <- length(IB_2018_Victoria) # Number of observations
Nx <- 2 # Number of variables: intercept + exposure 

# Data, parameter list and starting values
mu <- rep(0, Nx)
tau <- diag(0.001, Nx)
data.logistic <- list("N", "Nx", "IB_2017_Victoria", "IB_2018_Victoria", "mu", "tau")
parameters.logistic <-c("b","OR") # Parameters to keep track of
inits.logistic <- function() {list (b= rep(0, Nx))}

set.seed(114011)
logistic.sim<-jags(data=data.logistic, inits=inits.logistic, parameters.logistic, n.iter=100000,
                   n.burn=25000, model.file=logistic.model, n.thin=5, n.chains = 3)

print(logistic.sim, 2)
OR_IB_2017_2018_Victoria_Bayes_SA <- round(logistic.sim$BUGSoutput$summary[2,c(5,3,7)], 2)
OR_IB_2017_2018_Victoria_Bayes_SA

# Bayesian model diagnostics 
traceplot(logistic.sim, ask=TRUE)
autocorr.plot(logistic.sim) # ACF for all chains

logistic.mcmc <- as.mcmc(logistic.sim)

plot(logistic.mcmc)
autocorr.diag(logistic.mcmc)
autocorr.plot(logistic.mcmc[1]) # ACF for chain 1
autocorr.plot(logistic.mcmc[2]) # ACF for chain 2
autocorr.plot(logistic.mcmc[3]) # ACF for chain 3

geweke.diag(logistic.mcmc) # Geweke statistic
gelman.plot(logistic.mcmc,ask=TRUE) # Gelman-Rubin-Brooks plot

HPDinterval(logistic.mcmc)



# IB 2017 to 2019 - Victoria
try52 <- data_orig %>%
  filter(at2017==1 & at2019==1 & IB_2018_Victoria==0) %>%
  dplyr::select(Obs, IB_2017_Victoria, IB_2019_Victoria)
try52$IB_2017_Victoria <- as.numeric(try52$IB_2017_Victoria)
try52$IB_2019_Victoria <- as.numeric(try52$IB_2019_Victoria)

table(try52$IB_2017_Victoria, try52$IB_2019_Victoria)

# Bayesian analyses
# Putting exposure and outcome variables in the environment 
IB_2017_Victoria <- try52$IB_2017_Victoria
IB_2019_Victoria <- try52$IB_2019_Victoria

### Bayesian logistic regression
logistic.model <- function() {
  # SAMPLING DISTRIBUTION
  for (i in 1:N) {
    logit(pi[i]) <- b[1] + b[2]*IB_2017_Victoria[i];
    IB_2019_Victoria[i] ~ dbin(pi[i],1);  
  }
  
  # PRIORS ON BETAS
  b[1:Nx] ~ dmnorm(mu[1:Nx], tau[1:Nx,1:Nx]) # multivariate normal prior
  
  # Calculate ORs:
  for (l in 1:Nx) {
    OR[l] <- exp(b[l]);
  }
}

N <- length(IB_2019_Victoria) # Number of observations
Nx <- 2 # Number of variables: intercept + exposure 

# Data, parameter list and starting values
mu <- rep(0, Nx)
tau <- diag(0.001, Nx)
data.logistic <- list("N", "Nx", "IB_2017_Victoria", "IB_2019_Victoria", "mu", "tau")
parameters.logistic <-c("b","OR") # Parameters to keep track of
inits.logistic <- function() {list (b= rep(0, Nx))}

set.seed(114011)
logistic.sim<-jags(data=data.logistic, inits=inits.logistic, parameters.logistic, n.iter=100000,
                   n.burn=25000, model.file=logistic.model, n.thin=5, n.chains = 3)

print(logistic.sim, 2)
OR_IB_2017_2019_Victoria_Bayes_SA <- round(logistic.sim$BUGSoutput$summary[2,c(5,3,7)], 2)
OR_IB_2017_2019_Victoria_Bayes_SA

# sensitivity analysis

try52 <- data_orig %>%
  filter(at2017==1 & at2019==1 ) %>%
  dplyr::select(Obs, IB_2017_Victoria, IB_2019_Victoria)
try52$IB_2017_Victoria <- as.numeric(try52$IB_2017_Victoria)
try52$IB_2019_Victoria <- as.numeric(try52$IB_2019_Victoria)

table(try52$IB_2017_Victoria, try52$IB_2019_Victoria)

# Bayesian analyses
# Putting exposure and outcome variables in the environment 
IB_2017_Victoria <- try52$IB_2017_Victoria
IB_2019_Victoria <- try52$IB_2019_Victoria

### Bayesian logistic regression
logistic.model <- function() {
  # SAMPLING DISTRIBUTION
  for (i in 1:N) {
    logit(pi[i]) <- b[1] + b[2]*IB_2017_Victoria[i];
    IB_2019_Victoria[i] ~ dbin(pi[i],1);  
  }
  
  # PRIORS ON BETAS
  b[1:Nx] ~ dmnorm(mu[1:Nx], tau[1:Nx,1:Nx]) # multivariate normal prior
  
  # Calculate ORs:
  for (l in 1:Nx) {
    OR[l] <- exp(b[l]);
  }
}

N <- length(IB_2019_Victoria) # Number of observations
Nx <- 2 # Number of variables: intercept + exposure 

# Data, parameter list and starting values
mu <- rep(0, Nx)
tau <- diag(0.001, Nx)
data.logistic <- list("N", "Nx", "IB_2017_Victoria", "IB_2019_Victoria", "mu", "tau")
parameters.logistic <-c("b","OR") # Parameters to keep track of
inits.logistic <- function() {list (b= rep(0, Nx))}

set.seed(114011)
logistic.sim<-jags(data=data.logistic, inits=inits.logistic, parameters.logistic, n.iter=100000,
                   n.burn=25000, model.file=logistic.model, n.thin=5, n.chains = 3)

print(logistic.sim, 2)
OR_IB_2017_2019_Victoria_Bayes_SA <- round(logistic.sim$BUGSoutput$summary[2,c(5,3,7)], 2)
OR_IB_2017_2019_Victoria_Bayes_SA




# IB 2018 to 2019 - Victoria
try52 <- data_orig %>%
  filter(at2018==1 & at2019==1) %>%
  dplyr::select(Obs, IB_2018_Victoria, IB_2019_Victoria)
try52$IB_2018_Victoria <- as.numeric(try52$IB_2018_Victoria)
try52$IB_2019_Victoria <- as.numeric(try52$IB_2019_Victoria)

table(try52$IB_2018_Victoria, try52$IB_2019_Victoria)

# Bayesian analyses
# Putting exposure and outcome variables in the environment 
IB_2018_Victoria <- try52$IB_2018_Victoria
IB_2019_Victoria <- try52$IB_2019_Victoria

### Bayesian logistic regression
logistic.model <- function() {
  # SAMPLING DISTRIBUTION
  for (i in 1:N) {
    logit(pi[i]) <- b[1] + b[2]*IB_2018_Victoria[i];
    IB_2019_Victoria[i] ~ dbin(pi[i],1);  
  }
  
  # PRIORS ON BETAS
  b[1:Nx] ~ dmnorm(mu[1:Nx], tau[1:Nx,1:Nx]) # multivariate normal prior
  
  # Calculate ORs:
  for (l in 1:Nx) {
    OR[l] <- exp(b[l]);
  }
}

N <- length(IB_2019_Victoria) # Number of observations
Nx <- 2 # Number of variables: intercept + exposure 

# Data, parameter list and starting values
mu <- rep(0, Nx)
tau <- diag(0.001, Nx)
data.logistic <- list("N", "Nx", "IB_2018_Victoria", "IB_2019_Victoria", "mu", "tau")
parameters.logistic <-c("b","OR") # Parameters to keep track of
inits.logistic <- function() {list (b= rep(0, Nx))}

set.seed(114011)
logistic.sim<-jags(data=data.logistic, inits=inits.logistic, parameters.logistic, n.iter=100000,
                   n.burn=25000, model.file=logistic.model, n.thin=5, n.chains = 3)

print(logistic.sim, 2)
OR_IB_2018_2019_Victoria_Bayes_SA <- round(logistic.sim$BUGSoutput$summary[2,c(5,3,7)], 2)
OR_IB_2018_2019_Victoria_Bayes_SA


### AGE STRATIFIED

# Victoria 2012 to 2016 

try35 <- data_orig %>%
  filter(at2012==1 & at2016==1) %>%
  dplyr::select(Obs, IB_2012_Victoria, IB_2016_Victoria, age_2012_cat)

try35$IB_2012_Victoria <- as.numeric(try35$IB_2012_Victoria)
try35$IB_2016_Victoria <- as.numeric(try35$IB_2016_Victoria)

# 2x2 table
table(try35[try35$age_2012_cat==0,]$IB_2012_Victoria, try35[try35$age_2012_cat==0,]$IB_2016_Victoria) 
table(try35[try35$age_2012_cat==1,]$IB_2012_Victoria, try35[try35$age_2012_cat==1,]$IB_2016_Victoria) 


# Bayesian analysis of H3N2 2014 on H3N2 2016, under 5s
IB_2012_Victoria <- try35[try35$age_2012_cat==0,]$IB_2012_Victoria
IB_2016_Victoria <- try35[try35$age_2012_cat==0,]$IB_2016_Victoria

logistic.model <- function() {
  # SAMPLING DISTRIBUTION
  for (i in 1:N) {
    logit(pi[i]) <- b[1] + b[2]*IB_2012_Victoria[i];
    IB_2016_Victoria[i] ~ dbin(pi[i],1);  
  }
  
  # PRIORS ON BETAS
  b[1:Nx] ~ dmnorm(mu[1:Nx], tau[1:Nx,1:Nx]) # multivariate normal prior
  
  # Calculate ORs:
  for (l in 1:Nx) {
    OR[l] <- exp(b[l]);
  }
}

N <- length(IB_2016_Victoria) 
Nx <- 2 

mu <- rep(0, Nx)
tau <- diag(0.001, Nx)
data.logistic <- list("N", "Nx", "IB_2012_Victoria", "IB_2016_Victoria", "mu", "tau")
parameters.logistic <-c("b","OR") # Parameters to keep track of
inits.logistic <- function() {list (b= rep(0, Nx))}

set.seed(114011)
logistic.sim_IB_Victoria_2012_2016_cat0<-jags(data=data.logistic, inits=inits.logistic, parameters.logistic, n.iter=100000,
                                            n.burn=25000, model.file=logistic.model, n.thin=5, n.chains = 3)

print(logistic.sim_IB_Victoria_2012_2016_cat0, 2)
OR_IB_Victoria_2012_2016_cat0_Bayes <- round(logistic.sim_IB_Victoria_2012_2016_cat0$BUGSoutput$summary[2,c(5,3,7)], 2)
OR_IB_Victoria_2012_2016_cat0_Bayes

# Bayesian analysis of H3N2 2014 on H3N2 2016, over 5s
IB_2012_Victoria <- try35[try35$age_2012_cat==1,]$IB_2012_Victoria
IB_2016_Victoria <- try35[try35$age_2012_cat==1,]$IB_2016_Victoria

logistic.model <- function() {
  # SAMPLING DISTRIBUTION
  for (i in 1:N) {
    logit(pi[i]) <- b[1] + b[2]*IB_2012_Victoria[i];
    IB_2016_Victoria[i] ~ dbin(pi[i],1);  
  }
  
  # PRIORS ON BETAS
  b[1:Nx] ~ dmnorm(mu[1:Nx], tau[1:Nx,1:Nx]) # multivariate normal prior
  
  # Calculate ORs:
  for (l in 1:Nx) {
    OR[l] <- exp(b[l]);
  }
}

N <- length(IB_2016_Victoria) 
Nx <- 2 

mu <- rep(0, Nx)
tau <- diag(0.001, Nx)
data.logistic <- list("N", "Nx", "IB_2012_Victoria", "IB_2016_Victoria", "mu", "tau")
parameters.logistic <-c("b","OR") # Parameters to keep track of
inits.logistic <- function() {list (b= rep(0, Nx))}

set.seed(114011)
logistic.sim_IB_Victoria_2012_2016_cat1<-jags(data=data.logistic, inits=inits.logistic, parameters.logistic, n.iter=100000,
                                            n.burn=25000, model.file=logistic.model, n.thin=5, n.chains = 3)

print(logistic.sim_IB_Victoria_2012_2016_cat1, 2)
OR_IB_Victoria_2012_2016_cat1_Bayes <- round(logistic.sim_IB_Victoria_2012_2016_cat1$BUGSoutput$summary[2,c(5,3,7)], 2)
OR_IB_Victoria_2012_2016_cat1_Bayes


# Victoria 2012 to 2017 

try35 <- data_orig %>%
  filter(at2012==1 & at2017==1 & IB_2016_Victoria==0) %>%
  dplyr::select(Obs, IB_2012_Victoria, IB_2017_Victoria, age_2012_cat)

try35$IB_2012_Victoria <- as.numeric(try35$IB_2012_Victoria)
try35$IB_2017_Victoria <- as.numeric(try35$IB_2017_Victoria)

# 2x2 table
table(try35[try35$age_2012_cat==0,]$IB_2012_Victoria, try35[try35$age_2012_cat==0,]$IB_2017_Victoria) 
table(try35[try35$age_2012_cat==1,]$IB_2012_Victoria, try35[try35$age_2012_cat==1,]$IB_2017_Victoria) 


# Bayesian analysis of H3N2 2014 on H3N2 2016, under 5s
IB_2012_Victoria <- try35[try35$age_2012_cat==0,]$IB_2012_Victoria
IB_2017_Victoria <- try35[try35$age_2012_cat==0,]$IB_2017_Victoria

logistic.model <- function() {
  # SAMPLING DISTRIBUTION
  for (i in 1:N) {
    logit(pi[i]) <- b[1] + b[2]*IB_2012_Victoria[i];
    IB_2017_Victoria[i] ~ dbin(pi[i],1);  
  }
  
  # PRIORS ON BETAS
  b[1:Nx] ~ dmnorm(mu[1:Nx], tau[1:Nx,1:Nx]) # multivariate normal prior
  
  # Calculate ORs:
  for (l in 1:Nx) {
    OR[l] <- exp(b[l]);
  }
}

N <- length(IB_2017_Victoria) 
Nx <- 2 

mu <- rep(0, Nx)
tau <- diag(0.001, Nx)
data.logistic <- list("N", "Nx", "IB_2012_Victoria", "IB_2017_Victoria", "mu", "tau")
parameters.logistic <-c("b","OR") # Parameters to keep track of
inits.logistic <- function() {list (b= rep(0, Nx))}

set.seed(114011)
logistic.sim_IB_Victoria_2012_2017_cat0<-jags(data=data.logistic, inits=inits.logistic, parameters.logistic, n.iter=100000,
                                              n.burn=25000, model.file=logistic.model, n.thin=5, n.chains = 3)

print(logistic.sim_IB_Victoria_2012_2017_cat0, 2)
OR_IB_Victoria_2012_2017_cat0_Bayes <- round(logistic.sim_IB_Victoria_2012_2017_cat0$BUGSoutput$summary[2,c(5,3,7)], 2)
OR_IB_Victoria_2012_2017_cat0_Bayes

# Bayesian analysis of H3N2 2014 on H3N2 2016, over 5s
IB_2012_Victoria <- try35[try35$age_2012_cat==1,]$IB_2012_Victoria
IB_2017_Victoria <- try35[try35$age_2012_cat==1,]$IB_2017_Victoria

logistic.model <- function() {
  # SAMPLING DISTRIBUTION
  for (i in 1:N) {
    logit(pi[i]) <- b[1] + b[2]*IB_2012_Victoria[i];
    IB_2017_Victoria[i] ~ dbin(pi[i],1);  
  }
  
  # PRIORS ON BETAS
  b[1:Nx] ~ dmnorm(mu[1:Nx], tau[1:Nx,1:Nx]) # multivariate normal prior
  
  # Calculate ORs:
  for (l in 1:Nx) {
    OR[l] <- exp(b[l]);
  }
}

N <- length(IB_2017_Victoria) 
Nx <- 2 

mu <- rep(0, Nx)
tau <- diag(0.001, Nx)
data.logistic <- list("N", "Nx", "IB_2012_Victoria", "IB_2017_Victoria", "mu", "tau")
parameters.logistic <-c("b","OR") # Parameters to keep track of
inits.logistic <- function() {list (b= rep(0, Nx))}

set.seed(114011)
logistic.sim_IB_Victoria_2012_2017_cat1<-jags(data=data.logistic, inits=inits.logistic, parameters.logistic, n.iter=100000,
                                              n.burn=25000, model.file=logistic.model, n.thin=5, n.chains = 3)

print(logistic.sim_IB_Victoria_2012_2017_cat1, 2)
OR_IB_Victoria_2012_2017_cat1_Bayes <- round(logistic.sim_IB_Victoria_2012_2017_cat1$BUGSoutput$summary[2,c(5,3,7)], 2)
OR_IB_Victoria_2012_2017_cat1_Bayes

### sensitivity analysis


try35 <- data_orig %>%
  filter(at2012==1 & at2017==1) %>%
  dplyr::select(Obs, IB_2012_Victoria, IB_2017_Victoria, age_2012_cat)

try35$IB_2012_Victoria <- as.numeric(try35$IB_2012_Victoria)
try35$IB_2017_Victoria <- as.numeric(try35$IB_2017_Victoria)

# 2x2 table
table(try35[try35$age_2012_cat==0,]$IB_2012_Victoria, try35[try35$age_2012_cat==0,]$IB_2017_Victoria) 
table(try35[try35$age_2012_cat==1,]$IB_2012_Victoria, try35[try35$age_2012_cat==1,]$IB_2017_Victoria) 


# Bayesian analysis of H3N2 2014 on H3N2 2016, under 5s
IB_2012_Victoria <- try35[try35$age_2012_cat==0,]$IB_2012_Victoria
IB_2017_Victoria <- try35[try35$age_2012_cat==0,]$IB_2017_Victoria

logistic.model <- function() {
  # SAMPLING DISTRIBUTION
  for (i in 1:N) {
    logit(pi[i]) <- b[1] + b[2]*IB_2012_Victoria[i];
    IB_2017_Victoria[i] ~ dbin(pi[i],1);  
  }
  
  # PRIORS ON BETAS
  b[1:Nx] ~ dmnorm(mu[1:Nx], tau[1:Nx,1:Nx]) # multivariate normal prior
  
  # Calculate ORs:
  for (l in 1:Nx) {
    OR[l] <- exp(b[l]);
  }
}

N <- length(IB_2017_Victoria) 
Nx <- 2 

mu <- rep(0, Nx)
tau <- diag(0.001, Nx)
data.logistic <- list("N", "Nx", "IB_2012_Victoria", "IB_2017_Victoria", "mu", "tau")
parameters.logistic <-c("b","OR") # Parameters to keep track of
inits.logistic <- function() {list (b= rep(0, Nx))}

set.seed(114011)
logistic.sim_IB_Victoria_2012_2017_cat0<-jags(data=data.logistic, inits=inits.logistic, parameters.logistic, n.iter=100000,
                                              n.burn=25000, model.file=logistic.model, n.thin=5, n.chains = 3)

print(logistic.sim_IB_Victoria_2012_2017_cat0, 2)
OR_IB_Victoria_2012_2017_cat0_Bayes <- round(logistic.sim_IB_Victoria_2012_2017_cat0$BUGSoutput$summary[2,c(5,3,7)], 2)
OR_IB_Victoria_2012_2017_cat0_Bayes

# Bayesian analysis of H3N2 2014 on H3N2 2016, over 5s
IB_2012_Victoria <- try35[try35$age_2012_cat==1,]$IB_2012_Victoria
IB_2017_Victoria <- try35[try35$age_2012_cat==1,]$IB_2017_Victoria

logistic.model <- function() {
  # SAMPLING DISTRIBUTION
  for (i in 1:N) {
    logit(pi[i]) <- b[1] + b[2]*IB_2012_Victoria[i];
    IB_2017_Victoria[i] ~ dbin(pi[i],1);  
  }
  
  # PRIORS ON BETAS
  b[1:Nx] ~ dmnorm(mu[1:Nx], tau[1:Nx,1:Nx]) # multivariate normal prior
  
  # Calculate ORs:
  for (l in 1:Nx) {
    OR[l] <- exp(b[l]);
  }
}

N <- length(IB_2017_Victoria) 
Nx <- 2 

mu <- rep(0, Nx)
tau <- diag(0.001, Nx)
data.logistic <- list("N", "Nx", "IB_2012_Victoria", "IB_2017_Victoria", "mu", "tau")
parameters.logistic <-c("b","OR") # Parameters to keep track of
inits.logistic <- function() {list (b= rep(0, Nx))}

set.seed(114011)
logistic.sim_IB_Victoria_2012_2017_cat1<-jags(data=data.logistic, inits=inits.logistic, parameters.logistic, n.iter=100000,
                                              n.burn=25000, model.file=logistic.model, n.thin=5, n.chains = 3)

print(logistic.sim_IB_Victoria_2012_2017_cat1, 2)
OR_IB_Victoria_2012_2017_cat1_Bayes <- round(logistic.sim_IB_Victoria_2012_2017_cat1$BUGSoutput$summary[2,c(5,3,7)], 2)
OR_IB_Victoria_2012_2017_cat1_Bayes

# Victoria 2012 to 2018 


try35 <- data_orig %>%
  filter(at2012==1 & at2018==1 & IB_2016_Victoria==0 & IB_2017_Victoria==0) %>%
  dplyr::select(Obs, IB_2012_Victoria, IB_2018_Victoria, age_2012_cat)

try35$IB_2012_Victoria <- as.numeric(try35$IB_2012_Victoria)
try35$IB_2018_Victoria <- as.numeric(try35$IB_2018_Victoria)

# 2x2 table
table(try35[try35$age_2012_cat==0,]$IB_2012_Victoria, try35[try35$age_2012_cat==0,]$IB_2018_Victoria) 
table(try35[try35$age_2012_cat==1,]$IB_2012_Victoria, try35[try35$age_2012_cat==1,]$IB_2018_Victoria) 


# Bayesian analysis of H3N2 2014 on H3N2 2016, under 5s
IB_2012_Victoria <- try35[try35$age_2012_cat==0,]$IB_2012_Victoria
IB_2018_Victoria <- try35[try35$age_2012_cat==0,]$IB_2018_Victoria

logistic.model <- function() {
  # SAMPLING DISTRIBUTION
  for (i in 1:N) {
    logit(pi[i]) <- b[1] + b[2]*IB_2012_Victoria[i];
    IB_2018_Victoria[i] ~ dbin(pi[i],1);  
  }
  
  # PRIORS ON BETAS
  b[1:Nx] ~ dmnorm(mu[1:Nx], tau[1:Nx,1:Nx]) # multivariate normal prior
  
  # Calculate ORs:
  for (l in 1:Nx) {
    OR[l] <- exp(b[l]);
  }
}

N <- length(IB_2018_Victoria) 
Nx <- 2 

mu <- rep(0, Nx)
tau <- diag(0.001, Nx)
data.logistic <- list("N", "Nx", "IB_2012_Victoria", "IB_2018_Victoria", "mu", "tau")
parameters.logistic <-c("b","OR") # Parameters to keep track of
inits.logistic <- function() {list (b= rep(0, Nx))}

set.seed(114011)
logistic.sim_IB_Victoria_2012_2018_cat0<-jags(data=data.logistic, inits=inits.logistic, parameters.logistic, n.iter=100000,
                                              n.burn=25000, model.file=logistic.model, n.thin=5, n.chains = 3)

print(logistic.sim_IB_Victoria_2012_2018_cat0, 2)
OR_IB_Victoria_2012_2018_cat0_Bayes <- round(logistic.sim_IB_Victoria_2012_2018_cat0$BUGSoutput$summary[2,c(5,3,7)], 2)
OR_IB_Victoria_2012_2018_cat0_Bayes

# Bayesian analysis of H3N2 2014 on H3N2 2016, over 5s
IB_2012_Victoria <- try35[try35$age_2012_cat==1,]$IB_2012_Victoria
IB_2018_Victoria <- try35[try35$age_2012_cat==1,]$IB_2018_Victoria

logistic.model <- function() {
  # SAMPLING DISTRIBUTION
  for (i in 1:N) {
    logit(pi[i]) <- b[1] + b[2]*IB_2012_Victoria[i];
    IB_2018_Victoria[i] ~ dbin(pi[i],1);  
  }
  
  # PRIORS ON BETAS
  b[1:Nx] ~ dmnorm(mu[1:Nx], tau[1:Nx,1:Nx]) # multivariate normal prior
  
  # Calculate ORs:
  for (l in 1:Nx) {
    OR[l] <- exp(b[l]);
  }
}

N <- length(IB_2018_Victoria) 
Nx <- 2 

mu <- rep(0, Nx)
tau <- diag(0.001, Nx)
data.logistic <- list("N", "Nx", "IB_2012_Victoria", "IB_2018_Victoria", "mu", "tau")
parameters.logistic <-c("b","OR") # Parameters to keep track of
inits.logistic <- function() {list (b= rep(0, Nx))}

set.seed(114011)
logistic.sim_IB_Victoria_2012_2018_cat1<-jags(data=data.logistic, inits=inits.logistic, parameters.logistic, n.iter=100000,
                                              n.burn=25000, model.file=logistic.model, n.thin=5, n.chains = 3)

print(logistic.sim_IB_Victoria_2012_2018_cat1, 2)
OR_IB_Victoria_2012_2018_cat1_Bayes <- round(logistic.sim_IB_Victoria_2012_2018_cat1$BUGSoutput$summary[2,c(5,3,7)], 2)
OR_IB_Victoria_2012_2018_cat1_Bayes

### sensitivity analysis


try35 <- data_orig %>%
  filter(at2012==1 & at2018==1) %>%
  dplyr::select(Obs, IB_2012_Victoria, IB_2018_Victoria, age_2012_cat)

try35$IB_2012_Victoria <- as.numeric(try35$IB_2012_Victoria)
try35$IB_2018_Victoria <- as.numeric(try35$IB_2018_Victoria)

# 2x2 table
table(try35[try35$age_2012_cat==0,]$IB_2012_Victoria, try35[try35$age_2012_cat==0,]$IB_2018_Victoria) 
table(try35[try35$age_2012_cat==1,]$IB_2012_Victoria, try35[try35$age_2012_cat==1,]$IB_2018_Victoria) 


# Bayesian analysis of H3N2 2014 on H3N2 2016, under 5s
IB_2012_Victoria <- try35[try35$age_2012_cat==0,]$IB_2012_Victoria
IB_2018_Victoria <- try35[try35$age_2012_cat==0,]$IB_2018_Victoria

logistic.model <- function() {
  # SAMPLING DISTRIBUTION
  for (i in 1:N) {
    logit(pi[i]) <- b[1] + b[2]*IB_2012_Victoria[i];
    IB_2018_Victoria[i] ~ dbin(pi[i],1);  
  }
  
  # PRIORS ON BETAS
  b[1:Nx] ~ dmnorm(mu[1:Nx], tau[1:Nx,1:Nx]) # multivariate normal prior
  
  # Calculate ORs:
  for (l in 1:Nx) {
    OR[l] <- exp(b[l]);
  }
}

N <- length(IB_2018_Victoria) 
Nx <- 2 

mu <- rep(0, Nx)
tau <- diag(0.001, Nx)
data.logistic <- list("N", "Nx", "IB_2012_Victoria", "IB_2018_Victoria", "mu", "tau")
parameters.logistic <-c("b","OR") # Parameters to keep track of
inits.logistic <- function() {list (b= rep(0, Nx))}

set.seed(114011)
logistic.sim_IB_Victoria_2012_2018_cat0<-jags(data=data.logistic, inits=inits.logistic, parameters.logistic, n.iter=100000,
                                              n.burn=25000, model.file=logistic.model, n.thin=5, n.chains = 3)

print(logistic.sim_IB_Victoria_2012_2018_cat0, 2)
OR_IB_Victoria_2012_2018_cat0_Bayes <- round(logistic.sim_IB_Victoria_2012_2018_cat0$BUGSoutput$summary[2,c(5,3,7)], 2)
OR_IB_Victoria_2012_2018_cat0_Bayes

# Bayesian analysis of H3N2 2014 on H3N2 2016, over 5s
IB_2012_Victoria <- try35[try35$age_2012_cat==1,]$IB_2012_Victoria
IB_2018_Victoria <- try35[try35$age_2012_cat==1,]$IB_2018_Victoria

logistic.model <- function() {
  # SAMPLING DISTRIBUTION
  for (i in 1:N) {
    logit(pi[i]) <- b[1] + b[2]*IB_2012_Victoria[i];
    IB_2018_Victoria[i] ~ dbin(pi[i],1);  
  }
  
  # PRIORS ON BETAS
  b[1:Nx] ~ dmnorm(mu[1:Nx], tau[1:Nx,1:Nx]) # multivariate normal prior
  
  # Calculate ORs:
  for (l in 1:Nx) {
    OR[l] <- exp(b[l]);
  }
}

N <- length(IB_2018_Victoria) 
Nx <- 2 

mu <- rep(0, Nx)
tau <- diag(0.001, Nx)
data.logistic <- list("N", "Nx", "IB_2012_Victoria", "IB_2018_Victoria", "mu", "tau")
parameters.logistic <-c("b","OR") # Parameters to keep track of
inits.logistic <- function() {list (b= rep(0, Nx))}

set.seed(114011)
logistic.sim_IB_Victoria_2012_2018_cat1<-jags(data=data.logistic, inits=inits.logistic, parameters.logistic, n.iter=100000,
                                              n.burn=25000, model.file=logistic.model, n.thin=5, n.chains = 3)

print(logistic.sim_IB_Victoria_2012_2018_cat1, 2)
OR_IB_Victoria_2012_2018_cat1_Bayes <- round(logistic.sim_IB_Victoria_2012_2018_cat1$BUGSoutput$summary[2,c(5,3,7)], 2)
OR_IB_Victoria_2012_2018_cat1_Bayes


# Victoria 2012 to 2019 


try35 <- data_orig %>%
  filter(at2012==1 & at2019==1 & IB_2016_Victoria==0 & IB_2017_Victoria==0 & IB_2018_Victoria==0) %>%
  dplyr::select(Obs, IB_2012_Victoria, IB_2019_Victoria, age_2012_cat)

try35$IB_2012_Victoria <- as.numeric(try35$IB_2012_Victoria)
try35$IB_2019_Victoria <- as.numeric(try35$IB_2019_Victoria)

# 2x2 table
table(try35[try35$age_2012_cat==0,]$IB_2012_Victoria, try35[try35$age_2012_cat==0,]$IB_2019_Victoria) 
table(try35[try35$age_2012_cat==1,]$IB_2012_Victoria, try35[try35$age_2012_cat==1,]$IB_2019_Victoria) 


# Bayesian analysis of H3N2 2014 on H3N2 2016, under 5s
IB_2012_Victoria <- try35[try35$age_2012_cat==0,]$IB_2012_Victoria
IB_2019_Victoria <- try35[try35$age_2012_cat==0,]$IB_2019_Victoria

logistic.model <- function() {
  # SAMPLING DISTRIBUTION
  for (i in 1:N) {
    logit(pi[i]) <- b[1] + b[2]*IB_2012_Victoria[i];
    IB_2019_Victoria[i] ~ dbin(pi[i],1);  
  }
  
  # PRIORS ON BETAS
  b[1:Nx] ~ dmnorm(mu[1:Nx], tau[1:Nx,1:Nx]) # multivariate normal prior
  
  # Calculate ORs:
  for (l in 1:Nx) {
    OR[l] <- exp(b[l]);
  }
}

N <- length(IB_2019_Victoria) 
Nx <- 2 

mu <- rep(0, Nx)
tau <- diag(0.001, Nx)
data.logistic <- list("N", "Nx", "IB_2012_Victoria", "IB_2019_Victoria", "mu", "tau")
parameters.logistic <-c("b","OR") # Parameters to keep track of
inits.logistic <- function() {list (b= rep(0, Nx))}

set.seed(114011)
logistic.sim_IB_Victoria_2012_2019_cat0<-jags(data=data.logistic, inits=inits.logistic, parameters.logistic, n.iter=100000,
                                              n.burn=25000, model.file=logistic.model, n.thin=5, n.chains = 3)

print(logistic.sim_IB_Victoria_2012_2019_cat0, 2)
OR_IB_Victoria_2012_2019_cat0_Bayes <- round(logistic.sim_IB_Victoria_2012_2019_cat0$BUGSoutput$summary[2,c(5,3,7)], 2)
OR_IB_Victoria_2012_2019_cat0_Bayes

# Bayesian analysis of H3N2 2014 on H3N2 2016, over 5s
IB_2012_Victoria <- try35[try35$age_2012_cat==1,]$IB_2012_Victoria
IB_2019_Victoria <- try35[try35$age_2012_cat==1,]$IB_2019_Victoria

logistic.model <- function() {
  # SAMPLING DISTRIBUTION
  for (i in 1:N) {
    logit(pi[i]) <- b[1] + b[2]*IB_2012_Victoria[i];
    IB_2019_Victoria[i] ~ dbin(pi[i],1);  
  }
  
  # PRIORS ON BETAS
  b[1:Nx] ~ dmnorm(mu[1:Nx], tau[1:Nx,1:Nx]) # multivariate normal prior
  
  # Calculate ORs:
  for (l in 1:Nx) {
    OR[l] <- exp(b[l]);
  }
}

N <- length(IB_2019_Victoria) 
Nx <- 2 

mu <- rep(0, Nx)
tau <- diag(0.001, Nx)
data.logistic <- list("N", "Nx", "IB_2012_Victoria", "IB_2019_Victoria", "mu", "tau")
parameters.logistic <-c("b","OR") # Parameters to keep track of
inits.logistic <- function() {list (b= rep(0, Nx))}

set.seed(114011)
logistic.sim_IB_Victoria_2012_2019_cat1<-jags(data=data.logistic, inits=inits.logistic, parameters.logistic, n.iter=100000,
                                              n.burn=25000, model.file=logistic.model, n.thin=5, n.chains = 3)

print(logistic.sim_IB_Victoria_2012_2019_cat1, 2)
OR_IB_Victoria_2012_2019_cat1_Bayes <- round(logistic.sim_IB_Victoria_2012_2019_cat1$BUGSoutput$summary[2,c(5,3,7)], 2)
OR_IB_Victoria_2012_2019_cat1_Bayes

### sensitivity analysis


try35 <- data_orig %>%
  filter(at2012==1 & at2019==1 ) %>%
  dplyr::select(Obs, IB_2012_Victoria, IB_2019_Victoria, age_2012_cat)

try35$IB_2012_Victoria <- as.numeric(try35$IB_2012_Victoria)
try35$IB_2019_Victoria <- as.numeric(try35$IB_2019_Victoria)

# 2x2 table
table(try35[try35$age_2012_cat==0,]$IB_2012_Victoria, try35[try35$age_2012_cat==0,]$IB_2019_Victoria) 
table(try35[try35$age_2012_cat==1,]$IB_2012_Victoria, try35[try35$age_2012_cat==1,]$IB_2019_Victoria) 


# Bayesian analysis of H3N2 2014 on H3N2 2016, under 5s
IB_2012_Victoria <- try35[try35$age_2012_cat==0,]$IB_2012_Victoria
IB_2019_Victoria <- try35[try35$age_2012_cat==0,]$IB_2019_Victoria

logistic.model <- function() {
  # SAMPLING DISTRIBUTION
  for (i in 1:N) {
    logit(pi[i]) <- b[1] + b[2]*IB_2012_Victoria[i];
    IB_2019_Victoria[i] ~ dbin(pi[i],1);  
  }
  
  # PRIORS ON BETAS
  b[1:Nx] ~ dmnorm(mu[1:Nx], tau[1:Nx,1:Nx]) # multivariate normal prior
  
  # Calculate ORs:
  for (l in 1:Nx) {
    OR[l] <- exp(b[l]);
  }
}

N <- length(IB_2019_Victoria) 
Nx <- 2 

mu <- rep(0, Nx)
tau <- diag(0.001, Nx)
data.logistic <- list("N", "Nx", "IB_2012_Victoria", "IB_2019_Victoria", "mu", "tau")
parameters.logistic <-c("b","OR") # Parameters to keep track of
inits.logistic <- function() {list (b= rep(0, Nx))}

set.seed(114011)
logistic.sim_IB_Victoria_2012_2019_cat0<-jags(data=data.logistic, inits=inits.logistic, parameters.logistic, n.iter=100000,
                                              n.burn=25000, model.file=logistic.model, n.thin=5, n.chains = 3)

print(logistic.sim_IB_Victoria_2012_2019_cat0, 2)
OR_IB_Victoria_2012_2019_cat0_Bayes <- round(logistic.sim_IB_Victoria_2012_2019_cat0$BUGSoutput$summary[2,c(5,3,7)], 2)
OR_IB_Victoria_2012_2019_cat0_Bayes

# Bayesian analysis of H3N2 2014 on H3N2 2016, over 5s
IB_2012_Victoria <- try35[try35$age_2012_cat==1,]$IB_2012_Victoria
IB_2019_Victoria <- try35[try35$age_2012_cat==1,]$IB_2019_Victoria

logistic.model <- function() {
  # SAMPLING DISTRIBUTION
  for (i in 1:N) {
    logit(pi[i]) <- b[1] + b[2]*IB_2012_Victoria[i];
    IB_2019_Victoria[i] ~ dbin(pi[i],1);  
  }
  
  # PRIORS ON BETAS
  b[1:Nx] ~ dmnorm(mu[1:Nx], tau[1:Nx,1:Nx]) # multivariate normal prior
  
  # Calculate ORs:
  for (l in 1:Nx) {
    OR[l] <- exp(b[l]);
  }
}

N <- length(IB_2019_Victoria) 
Nx <- 2 

mu <- rep(0, Nx)
tau <- diag(0.001, Nx)
data.logistic <- list("N", "Nx", "IB_2012_Victoria", "IB_2019_Victoria", "mu", "tau")
parameters.logistic <-c("b","OR") # Parameters to keep track of
inits.logistic <- function() {list (b= rep(0, Nx))}

set.seed(114011)
logistic.sim_IB_Victoria_2012_2019_cat1<-jags(data=data.logistic, inits=inits.logistic, parameters.logistic, n.iter=100000,
                                              n.burn=25000, model.file=logistic.model, n.thin=5, n.chains = 3)

print(logistic.sim_IB_Victoria_2012_2019_cat1, 2)
OR_IB_Victoria_2012_2019_cat1_Bayes <- round(logistic.sim_IB_Victoria_2012_2019_cat1$BUGSoutput$summary[2,c(5,3,7)], 2)
OR_IB_Victoria_2012_2019_cat1_Bayes


# Victoria 2016 to 2017 

try35 <- data_orig %>%
  filter(at2016==1 & at2017==1 ) %>%
  dplyr::select(Obs, IB_2016_Victoria, IB_2017_Victoria, age_2016_cat)

try35$IB_2016_Victoria <- as.numeric(try35$IB_2016_Victoria)
try35$IB_2017_Victoria <- as.numeric(try35$IB_2017_Victoria)

# 2x2 table
table(try35[try35$age_2016_cat==0,]$IB_2016_Victoria, try35[try35$age_2016_cat==0,]$IB_2017_Victoria) 
table(try35[try35$age_2016_cat==1,]$IB_2016_Victoria, try35[try35$age_2016_cat==1,]$IB_2017_Victoria) 


# Bayesian analysis of H3N2 2014 on H3N2 2016, under 5s
IB_2016_Victoria <- try35[try35$age_2016_cat==0,]$IB_2016_Victoria
IB_2017_Victoria <- try35[try35$age_2016_cat==0,]$IB_2017_Victoria

logistic.model <- function() {
  # SAMPLING DISTRIBUTION
  for (i in 1:N) {
    logit(pi[i]) <- b[1] + b[2]*IB_2016_Victoria[i];
    IB_2017_Victoria[i] ~ dbin(pi[i],1);  
  }
  
  # PRIORS ON BETAS
  b[1:Nx] ~ dmnorm(mu[1:Nx], tau[1:Nx,1:Nx]) # multivariate normal prior
  
  # Calculate ORs:
  for (l in 1:Nx) {
    OR[l] <- exp(b[l]);
  }
}

N <- length(IB_2017_Victoria) 
Nx <- 2 

mu <- rep(0, Nx)
tau <- diag(0.001, Nx)
data.logistic <- list("N", "Nx", "IB_2016_Victoria", "IB_2017_Victoria", "mu", "tau")
parameters.logistic <-c("b","OR") # Parameters to keep track of
inits.logistic <- function() {list (b= rep(0, Nx))}

set.seed(114011)
logistic.sim_IB_Victoria_2016_2017_cat0<-jags(data=data.logistic, inits=inits.logistic, parameters.logistic, n.iter=100000,
                                              n.burn=25000, model.file=logistic.model, n.thin=5, n.chains = 3)

print(logistic.sim_IB_Victoria_2016_2017_cat0, 2)
OR_IB_Victoria_2016_2017_cat0_Bayes <- round(logistic.sim_IB_Victoria_2016_2017_cat0$BUGSoutput$summary[2,c(5,3,7)], 2)
OR_IB_Victoria_2016_2017_cat0_Bayes

# Bayesian analysis of H3N2 2014 on H3N2 2016, over 5s
IB_2016_Victoria <- try35[try35$age_2016_cat==1,]$IB_2016_Victoria
IB_2017_Victoria <- try35[try35$age_2016_cat==1,]$IB_2017_Victoria

logistic.model <- function() {
  # SAMPLING DISTRIBUTION
  for (i in 1:N) {
    logit(pi[i]) <- b[1] + b[2]*IB_2016_Victoria[i];
    IB_2017_Victoria[i] ~ dbin(pi[i],1);  
  }
  
  # PRIORS ON BETAS
  b[1:Nx] ~ dmnorm(mu[1:Nx], tau[1:Nx,1:Nx]) # multivariate normal prior
  
  # Calculate ORs:
  for (l in 1:Nx) {
    OR[l] <- exp(b[l]);
  }
}

N <- length(IB_2017_Victoria) 
Nx <- 2 

mu <- rep(0, Nx)
tau <- diag(0.001, Nx)
data.logistic <- list("N", "Nx", "IB_2016_Victoria", "IB_2017_Victoria", "mu", "tau")
parameters.logistic <-c("b","OR") # Parameters to keep track of
inits.logistic <- function() {list (b= rep(0, Nx))}

set.seed(114011)
logistic.sim_IB_Victoria_2016_2017_cat1<-jags(data=data.logistic, inits=inits.logistic, parameters.logistic, n.iter=100000,
                                              n.burn=25000, model.file=logistic.model, n.thin=5, n.chains = 3)

print(logistic.sim_IB_Victoria_2016_2017_cat1, 2)
OR_IB_Victoria_2016_2017_cat1_Bayes <- round(logistic.sim_IB_Victoria_2016_2017_cat1$BUGSoutput$summary[2,c(5,3,7)], 2)
OR_IB_Victoria_2016_2017_cat1_Bayes

# Victoria 2016 to 2018

try35 <- data_orig %>%
  filter(at2016==1 & at2018==1 & IB_2017_Victoria==0) %>%
  dplyr::select(Obs, IB_2016_Victoria, IB_2018_Victoria, age_2016_cat)

try35$IB_2016_Victoria <- as.numeric(try35$IB_2016_Victoria)
try35$IB_2018_Victoria <- as.numeric(try35$IB_2018_Victoria)

# 2x2 table
table(try35[try35$age_2016_cat==0,]$IB_2016_Victoria, try35[try35$age_2016_cat==0,]$IB_2018_Victoria) 
table(try35[try35$age_2016_cat==1,]$IB_2016_Victoria, try35[try35$age_2016_cat==1,]$IB_2018_Victoria) 


# Bayesian analysis of H3N2 2014 on H3N2 2016, under 5s
IB_2016_Victoria <- try35[try35$age_2016_cat==0,]$IB_2016_Victoria
IB_2018_Victoria <- try35[try35$age_2016_cat==0,]$IB_2018_Victoria

logistic.model <- function() {
  # SAMPLING DISTRIBUTION
  for (i in 1:N) {
    logit(pi[i]) <- b[1] + b[2]*IB_2016_Victoria[i];
    IB_2018_Victoria[i] ~ dbin(pi[i],1);  
  }
  
  # PRIORS ON BETAS
  b[1:Nx] ~ dmnorm(mu[1:Nx], tau[1:Nx,1:Nx]) # multivariate normal prior
  
  # Calculate ORs:
  for (l in 1:Nx) {
    OR[l] <- exp(b[l]);
  }
}

N <- length(IB_2018_Victoria) 
Nx <- 2 

mu <- rep(0, Nx)
tau <- diag(0.001, Nx)
data.logistic <- list("N", "Nx", "IB_2016_Victoria", "IB_2018_Victoria", "mu", "tau")
parameters.logistic <-c("b","OR") # Parameters to keep track of
inits.logistic <- function() {list (b= rep(0, Nx))}

set.seed(114011)
logistic.sim_IB_Victoria_2016_2018_cat0<-jags(data=data.logistic, inits=inits.logistic, parameters.logistic, n.iter=100000,
                                              n.burn=25000, model.file=logistic.model, n.thin=5, n.chains = 3)

print(logistic.sim_IB_Victoria_2016_2018_cat0, 2)
OR_IB_Victoria_2016_2018_cat0_Bayes <- round(logistic.sim_IB_Victoria_2016_2018_cat0$BUGSoutput$summary[2,c(5,3,7)], 2)
OR_IB_Victoria_2016_2018_cat0_Bayes

# Bayesian analysis of H3N2 2014 on H3N2 2016, over 5s
IB_2016_Victoria <- try35[try35$age_2016_cat==1,]$IB_2016_Victoria
IB_2018_Victoria <- try35[try35$age_2016_cat==1,]$IB_2018_Victoria

logistic.model <- function() {
  # SAMPLING DISTRIBUTION
  for (i in 1:N) {
    logit(pi[i]) <- b[1] + b[2]*IB_2016_Victoria[i];
    IB_2018_Victoria[i] ~ dbin(pi[i],1);  
  }
  
  # PRIORS ON BETAS
  b[1:Nx] ~ dmnorm(mu[1:Nx], tau[1:Nx,1:Nx]) # multivariate normal prior
  
  # Calculate ORs:
  for (l in 1:Nx) {
    OR[l] <- exp(b[l]);
  }
}

N <- length(IB_2018_Victoria) 
Nx <- 2 

mu <- rep(0, Nx)
tau <- diag(0.001, Nx)
data.logistic <- list("N", "Nx", "IB_2016_Victoria", "IB_2018_Victoria", "mu", "tau")
parameters.logistic <-c("b","OR") # Parameters to keep track of
inits.logistic <- function() {list (b= rep(0, Nx))}

set.seed(114011)
logistic.sim_IB_Victoria_2016_2018_cat1<-jags(data=data.logistic, inits=inits.logistic, parameters.logistic, n.iter=100000,
                                              n.burn=25000, model.file=logistic.model, n.thin=5, n.chains = 3)

print(logistic.sim_IB_Victoria_2016_2018_cat1, 2)
OR_IB_Victoria_2016_2018_cat1_Bayes <- round(logistic.sim_IB_Victoria_2016_2018_cat1$BUGSoutput$summary[2,c(5,3,7)], 2)
OR_IB_Victoria_2016_2018_cat1_Bayes

### sensitivity analysis

try35 <- data_orig %>%
  filter(at2016==1 & at2018==1 ) %>%
  dplyr::select(Obs, IB_2016_Victoria, IB_2018_Victoria, age_2016_cat)

try35$IB_2016_Victoria <- as.numeric(try35$IB_2016_Victoria)
try35$IB_2018_Victoria <- as.numeric(try35$IB_2018_Victoria)

# 2x2 table
table(try35[try35$age_2016_cat==0,]$IB_2016_Victoria, try35[try35$age_2016_cat==0,]$IB_2018_Victoria) 
table(try35[try35$age_2016_cat==1,]$IB_2016_Victoria, try35[try35$age_2016_cat==1,]$IB_2018_Victoria) 


# Bayesian analysis of H3N2 2014 on H3N2 2016, under 5s
IB_2016_Victoria <- try35[try35$age_2016_cat==0,]$IB_2016_Victoria
IB_2018_Victoria <- try35[try35$age_2016_cat==0,]$IB_2018_Victoria

logistic.model <- function() {
  # SAMPLING DISTRIBUTION
  for (i in 1:N) {
    logit(pi[i]) <- b[1] + b[2]*IB_2016_Victoria[i];
    IB_2018_Victoria[i] ~ dbin(pi[i],1);  
  }
  
  # PRIORS ON BETAS
  b[1:Nx] ~ dmnorm(mu[1:Nx], tau[1:Nx,1:Nx]) # multivariate normal prior
  
  # Calculate ORs:
  for (l in 1:Nx) {
    OR[l] <- exp(b[l]);
  }
}

N <- length(IB_2018_Victoria) 
Nx <- 2 

mu <- rep(0, Nx)
tau <- diag(0.001, Nx)
data.logistic <- list("N", "Nx", "IB_2016_Victoria", "IB_2018_Victoria", "mu", "tau")
parameters.logistic <-c("b","OR") # Parameters to keep track of
inits.logistic <- function() {list (b= rep(0, Nx))}

set.seed(114011)
logistic.sim_IB_Victoria_2016_2018_cat0<-jags(data=data.logistic, inits=inits.logistic, parameters.logistic, n.iter=100000,
                                              n.burn=25000, model.file=logistic.model, n.thin=5, n.chains = 3)

print(logistic.sim_IB_Victoria_2016_2018_cat0, 2)
OR_IB_Victoria_2016_2018_cat0_Bayes <- round(logistic.sim_IB_Victoria_2016_2018_cat0$BUGSoutput$summary[2,c(5,3,7)], 2)
OR_IB_Victoria_2016_2018_cat0_Bayes

# Bayesian analysis of H3N2 2014 on H3N2 2016, over 5s
IB_2016_Victoria <- try35[try35$age_2016_cat==1,]$IB_2016_Victoria
IB_2018_Victoria <- try35[try35$age_2016_cat==1,]$IB_2018_Victoria

logistic.model <- function() {
  # SAMPLING DISTRIBUTION
  for (i in 1:N) {
    logit(pi[i]) <- b[1] + b[2]*IB_2016_Victoria[i];
    IB_2018_Victoria[i] ~ dbin(pi[i],1);  
  }
  
  # PRIORS ON BETAS
  b[1:Nx] ~ dmnorm(mu[1:Nx], tau[1:Nx,1:Nx]) # multivariate normal prior
  
  # Calculate ORs:
  for (l in 1:Nx) {
    OR[l] <- exp(b[l]);
  }
}

N <- length(IB_2018_Victoria) 
Nx <- 2 

mu <- rep(0, Nx)
tau <- diag(0.001, Nx)
data.logistic <- list("N", "Nx", "IB_2016_Victoria", "IB_2018_Victoria", "mu", "tau")
parameters.logistic <-c("b","OR") # Parameters to keep track of
inits.logistic <- function() {list (b= rep(0, Nx))}

set.seed(114011)
logistic.sim_IB_Victoria_2016_2018_cat1<-jags(data=data.logistic, inits=inits.logistic, parameters.logistic, n.iter=100000,
                                              n.burn=25000, model.file=logistic.model, n.thin=5, n.chains = 3)

print(logistic.sim_IB_Victoria_2016_2018_cat1, 2)
OR_IB_Victoria_2016_2018_cat1_Bayes <- round(logistic.sim_IB_Victoria_2016_2018_cat1$BUGSoutput$summary[2,c(5,3,7)], 2)
OR_IB_Victoria_2016_2018_cat1_Bayes


# Victoria 2016 to 2019

try35 <- data_orig %>%
  filter(at2016==1 & at2019==1 & IB_2017_Victoria==0 & IB_2018_Victoria==0) %>%
  dplyr::select(Obs, IB_2016_Victoria, IB_2019_Victoria, age_2016_cat)

try35$IB_2016_Victoria <- as.numeric(try35$IB_2016_Victoria)
try35$IB_2019_Victoria <- as.numeric(try35$IB_2019_Victoria)

# 2x2 table
table(try35[try35$age_2016_cat==0,]$IB_2016_Victoria, try35[try35$age_2016_cat==0,]$IB_2019_Victoria) 
table(try35[try35$age_2016_cat==1,]$IB_2016_Victoria, try35[try35$age_2016_cat==1,]$IB_2019_Victoria) 


# Bayesian analysis of H3N2 2014 on H3N2 2016, under 5s
IB_2016_Victoria <- try35[try35$age_2016_cat==0,]$IB_2016_Victoria
IB_2019_Victoria <- try35[try35$age_2016_cat==0,]$IB_2019_Victoria

logistic.model <- function() {
  # SAMPLING DISTRIBUTION
  for (i in 1:N) {
    logit(pi[i]) <- b[1] + b[2]*IB_2016_Victoria[i];
    IB_2019_Victoria[i] ~ dbin(pi[i],1);  
  }
  
  # PRIORS ON BETAS
  b[1:Nx] ~ dmnorm(mu[1:Nx], tau[1:Nx,1:Nx]) # multivariate normal prior
  
  # Calculate ORs:
  for (l in 1:Nx) {
    OR[l] <- exp(b[l]);
  }
}

N <- length(IB_2019_Victoria) 
Nx <- 2 

mu <- rep(0, Nx)
tau <- diag(0.001, Nx)
data.logistic <- list("N", "Nx", "IB_2016_Victoria", "IB_2019_Victoria", "mu", "tau")
parameters.logistic <-c("b","OR") # Parameters to keep track of
inits.logistic <- function() {list (b= rep(0, Nx))}

set.seed(114011)
logistic.sim_IB_Victoria_2016_2019_cat0<-jags(data=data.logistic, inits=inits.logistic, parameters.logistic, n.iter=100000,
                                              n.burn=25000, model.file=logistic.model, n.thin=5, n.chains = 3)

print(logistic.sim_IB_Victoria_2016_2019_cat0, 2)
OR_IB_Victoria_2016_2019_cat0_Bayes <- round(logistic.sim_IB_Victoria_2016_2019_cat0$BUGSoutput$summary[2,c(5,3,7)], 2)
OR_IB_Victoria_2016_2019_cat0_Bayes

# Bayesian analysis of H3N2 2014 on H3N2 2016, over 5s
IB_2016_Victoria <- try35[try35$age_2016_cat==1,]$IB_2016_Victoria
IB_2019_Victoria <- try35[try35$age_2016_cat==1,]$IB_2019_Victoria

logistic.model <- function() {
  # SAMPLING DISTRIBUTION
  for (i in 1:N) {
    logit(pi[i]) <- b[1] + b[2]*IB_2016_Victoria[i];
    IB_2019_Victoria[i] ~ dbin(pi[i],1);  
  }
  
  # PRIORS ON BETAS
  b[1:Nx] ~ dmnorm(mu[1:Nx], tau[1:Nx,1:Nx]) # multivariate normal prior
  
  # Calculate ORs:
  for (l in 1:Nx) {
    OR[l] <- exp(b[l]);
  }
}

N <- length(IB_2019_Victoria) 
Nx <- 2 

mu <- rep(0, Nx)
tau <- diag(0.001, Nx)
data.logistic <- list("N", "Nx", "IB_2016_Victoria", "IB_2019_Victoria", "mu", "tau")
parameters.logistic <-c("b","OR") # Parameters to keep track of
inits.logistic <- function() {list (b= rep(0, Nx))}

set.seed(114011)
logistic.sim_IB_Victoria_2016_2019_cat1<-jags(data=data.logistic, inits=inits.logistic, parameters.logistic, n.iter=100000,
                                              n.burn=25000, model.file=logistic.model, n.thin=5, n.chains = 3)

print(logistic.sim_IB_Victoria_2016_2019_cat1, 2)
OR_IB_Victoria_2016_2019_cat1_Bayes <- round(logistic.sim_IB_Victoria_2016_2019_cat1$BUGSoutput$summary[2,c(5,3,7)], 2)
OR_IB_Victoria_2016_2019_cat1_Bayes

### sensitivity analysis

try35 <- data_orig %>%
  filter(at2016==1 & at2019==1) %>%
  dplyr::select(Obs, IB_2016_Victoria, IB_2019_Victoria, age_2016_cat)

try35$IB_2016_Victoria <- as.numeric(try35$IB_2016_Victoria)
try35$IB_2019_Victoria <- as.numeric(try35$IB_2019_Victoria)

# 2x2 table
table(try35[try35$age_2016_cat==0,]$IB_2016_Victoria, try35[try35$age_2016_cat==0,]$IB_2019_Victoria) 
table(try35[try35$age_2016_cat==1,]$IB_2016_Victoria, try35[try35$age_2016_cat==1,]$IB_2019_Victoria) 


# Bayesian analysis of H3N2 2014 on H3N2 2016, under 5s
IB_2016_Victoria <- try35[try35$age_2016_cat==0,]$IB_2016_Victoria
IB_2019_Victoria <- try35[try35$age_2016_cat==0,]$IB_2019_Victoria

logistic.model <- function() {
  # SAMPLING DISTRIBUTION
  for (i in 1:N) {
    logit(pi[i]) <- b[1] + b[2]*IB_2016_Victoria[i];
    IB_2019_Victoria[i] ~ dbin(pi[i],1);  
  }
  
  # PRIORS ON BETAS
  b[1:Nx] ~ dmnorm(mu[1:Nx], tau[1:Nx,1:Nx]) # multivariate normal prior
  
  # Calculate ORs:
  for (l in 1:Nx) {
    OR[l] <- exp(b[l]);
  }
}

N <- length(IB_2019_Victoria) 
Nx <- 2 

mu <- rep(0, Nx)
tau <- diag(0.001, Nx)
data.logistic <- list("N", "Nx", "IB_2016_Victoria", "IB_2019_Victoria", "mu", "tau")
parameters.logistic <-c("b","OR") # Parameters to keep track of
inits.logistic <- function() {list (b= rep(0, Nx))}

set.seed(114011)
logistic.sim_IB_Victoria_2016_2019_cat0<-jags(data=data.logistic, inits=inits.logistic, parameters.logistic, n.iter=100000,
                                              n.burn=25000, model.file=logistic.model, n.thin=5, n.chains = 3)

print(logistic.sim_IB_Victoria_2016_2019_cat0, 2)
OR_IB_Victoria_2016_2019_cat0_Bayes <- round(logistic.sim_IB_Victoria_2016_2019_cat0$BUGSoutput$summary[2,c(5,3,7)], 2)
OR_IB_Victoria_2016_2019_cat0_Bayes

# Bayesian analysis of H3N2 2014 on H3N2 2016, over 5s
IB_2016_Victoria <- try35[try35$age_2016_cat==1,]$IB_2016_Victoria
IB_2019_Victoria <- try35[try35$age_2016_cat==1,]$IB_2019_Victoria

logistic.model <- function() {
  # SAMPLING DISTRIBUTION
  for (i in 1:N) {
    logit(pi[i]) <- b[1] + b[2]*IB_2016_Victoria[i];
    IB_2019_Victoria[i] ~ dbin(pi[i],1);  
  }
  
  # PRIORS ON BETAS
  b[1:Nx] ~ dmnorm(mu[1:Nx], tau[1:Nx,1:Nx]) # multivariate normal prior
  
  # Calculate ORs:
  for (l in 1:Nx) {
    OR[l] <- exp(b[l]);
  }
}

N <- length(IB_2019_Victoria) 
Nx <- 2 

mu <- rep(0, Nx)
tau <- diag(0.001, Nx)
data.logistic <- list("N", "Nx", "IB_2016_Victoria", "IB_2019_Victoria", "mu", "tau")
parameters.logistic <-c("b","OR") # Parameters to keep track of
inits.logistic <- function() {list (b= rep(0, Nx))}

set.seed(114011)
logistic.sim_IB_Victoria_2016_2019_cat1<-jags(data=data.logistic, inits=inits.logistic, parameters.logistic, n.iter=100000,
                                              n.burn=25000, model.file=logistic.model, n.thin=5, n.chains = 3)

print(logistic.sim_IB_Victoria_2016_2019_cat1, 2)
OR_IB_Victoria_2016_2019_cat1_Bayes <- round(logistic.sim_IB_Victoria_2016_2019_cat1$BUGSoutput$summary[2,c(5,3,7)], 2)
OR_IB_Victoria_2016_2019_cat1_Bayes


# Victoria 2017 to 2018

try35 <- data_orig %>%
  filter(at2017==1 & at2018==1 ) %>%
  dplyr::select(Obs, IB_2017_Victoria, IB_2018_Victoria, age_2017_cat)

try35$IB_2017_Victoria <- as.numeric(try35$IB_2017_Victoria)
try35$IB_2018_Victoria <- as.numeric(try35$IB_2018_Victoria)

# 2x2 table
table(try35[try35$age_2017_cat==0,]$IB_2017_Victoria, try35[try35$age_2017_cat==0,]$IB_2018_Victoria) 
table(try35[try35$age_2017_cat==1,]$IB_2017_Victoria, try35[try35$age_2017_cat==1,]$IB_2018_Victoria) 


# Bayesian analysis of H3N2 2014 on H3N2 2016, under 5s
IB_2017_Victoria <- try35[try35$age_2017_cat==0,]$IB_2017_Victoria
IB_2018_Victoria <- try35[try35$age_2017_cat==0,]$IB_2018_Victoria

logistic.model <- function() {
  # SAMPLING DISTRIBUTION
  for (i in 1:N) {
    logit(pi[i]) <- b[1] + b[2]*IB_2017_Victoria[i];
    IB_2018_Victoria[i] ~ dbin(pi[i],1);  
  }
  
  # PRIORS ON BETAS
  b[1:Nx] ~ dmnorm(mu[1:Nx], tau[1:Nx,1:Nx]) # multivariate normal prior
  
  # Calculate ORs:
  for (l in 1:Nx) {
    OR[l] <- exp(b[l]);
  }
}

N <- length(IB_2018_Victoria) 
Nx <- 2 

mu <- rep(0, Nx)
tau <- diag(0.001, Nx)
data.logistic <- list("N", "Nx", "IB_2017_Victoria", "IB_2018_Victoria", "mu", "tau")
parameters.logistic <-c("b","OR") # Parameters to keep track of
inits.logistic <- function() {list (b= rep(0, Nx))}

set.seed(114011)
logistic.sim_IB_Victoria_2017_2018_cat0<-jags(data=data.logistic, inits=inits.logistic, parameters.logistic, n.iter=100000,
                                              n.burn=25000, model.file=logistic.model, n.thin=5, n.chains = 3)

print(logistic.sim_IB_Victoria_2017_2018_cat0, 2)
OR_IB_Victoria_2017_2018_cat0_Bayes <- round(logistic.sim_IB_Victoria_2017_2018_cat0$BUGSoutput$summary[2,c(5,3,7)], 2)
OR_IB_Victoria_2017_2018_cat0_Bayes

# Bayesian analysis of H3N2 2014 on H3N2 2016, over 5s
IB_2017_Victoria <- try35[try35$age_2017_cat==1,]$IB_2017_Victoria
IB_2018_Victoria <- try35[try35$age_2017_cat==1,]$IB_2018_Victoria

logistic.model <- function() {
  # SAMPLING DISTRIBUTION
  for (i in 1:N) {
    logit(pi[i]) <- b[1] + b[2]*IB_2017_Victoria[i];
    IB_2018_Victoria[i] ~ dbin(pi[i],1);  
  }
  
  # PRIORS ON BETAS
  b[1:Nx] ~ dmnorm(mu[1:Nx], tau[1:Nx,1:Nx]) # multivariate normal prior
  
  # Calculate ORs:
  for (l in 1:Nx) {
    OR[l] <- exp(b[l]);
  }
}

N <- length(IB_2018_Victoria) 
Nx <- 2 

mu <- rep(0, Nx)
tau <- diag(0.001, Nx)
data.logistic <- list("N", "Nx", "IB_2017_Victoria", "IB_2018_Victoria", "mu", "tau")
parameters.logistic <-c("b","OR") # Parameters to keep track of
inits.logistic <- function() {list (b= rep(0, Nx))}

set.seed(114011)
logistic.sim_IB_Victoria_2017_2018_cat1<-jags(data=data.logistic, inits=inits.logistic, parameters.logistic, n.iter=100000,
                                              n.burn=25000, model.file=logistic.model, n.thin=5, n.chains = 3)

print(logistic.sim_IB_Victoria_2017_2018_cat1, 2)
OR_IB_Victoria_2017_2018_cat1_Bayes <- round(logistic.sim_IB_Victoria_2017_2018_cat1$BUGSoutput$summary[2,c(5,3,7)], 2)
OR_IB_Victoria_2017_2018_cat1_Bayes

# Victoria 2017 to 2019

try35 <- data_orig %>%
  filter(at2017==1 & at2019==1 & IB_2018_Victoria==0) %>%
  dplyr::select(Obs, IB_2017_Victoria, IB_2019_Victoria, age_2017_cat)

try35$IB_2017_Victoria <- as.numeric(try35$IB_2017_Victoria)
try35$IB_2019_Victoria <- as.numeric(try35$IB_2019_Victoria)

# 2x2 table
table(try35[try35$age_2017_cat==0,]$IB_2017_Victoria, try35[try35$age_2017_cat==0,]$IB_2019_Victoria) 
table(try35[try35$age_2017_cat==1,]$IB_2017_Victoria, try35[try35$age_2017_cat==1,]$IB_2019_Victoria) 


# Bayesian analysis of H3N2 2014 on H3N2 2016, under 5s
IB_2017_Victoria <- try35[try35$age_2017_cat==0,]$IB_2017_Victoria
IB_2019_Victoria <- try35[try35$age_2017_cat==0,]$IB_2019_Victoria

logistic.model <- function() {
  # SAMPLING DISTRIBUTION
  for (i in 1:N) {
    logit(pi[i]) <- b[1] + b[2]*IB_2017_Victoria[i];
    IB_2019_Victoria[i] ~ dbin(pi[i],1);  
  }
  
  # PRIORS ON BETAS
  b[1:Nx] ~ dmnorm(mu[1:Nx], tau[1:Nx,1:Nx]) # multivariate normal prior
  
  # Calculate ORs:
  for (l in 1:Nx) {
    OR[l] <- exp(b[l]);
  }
}

N <- length(IB_2019_Victoria) 
Nx <- 2 

mu <- rep(0, Nx)
tau <- diag(0.001, Nx)
data.logistic <- list("N", "Nx", "IB_2017_Victoria", "IB_2019_Victoria", "mu", "tau")
parameters.logistic <-c("b","OR") # Parameters to keep track of
inits.logistic <- function() {list (b= rep(0, Nx))}

set.seed(114011)
logistic.sim_IB_Victoria_2017_2019_cat0<-jags(data=data.logistic, inits=inits.logistic, parameters.logistic, n.iter=100000,
                                              n.burn=25000, model.file=logistic.model, n.thin=5, n.chains = 3)

print(logistic.sim_IB_Victoria_2017_2019_cat0, 2)
OR_IB_Victoria_2017_2019_cat0_Bayes <- round(logistic.sim_IB_Victoria_2017_2019_cat0$BUGSoutput$summary[2,c(5,3,7)], 2)
OR_IB_Victoria_2017_2019_cat0_Bayes

# Bayesian analysis of H3N2 2014 on H3N2 2016, over 5s
IB_2017_Victoria <- try35[try35$age_2017_cat==1,]$IB_2017_Victoria
IB_2019_Victoria <- try35[try35$age_2017_cat==1,]$IB_2019_Victoria

logistic.model <- function() {
  # SAMPLING DISTRIBUTION
  for (i in 1:N) {
    logit(pi[i]) <- b[1] + b[2]*IB_2017_Victoria[i];
    IB_2019_Victoria[i] ~ dbin(pi[i],1);  
  }
  
  # PRIORS ON BETAS
  b[1:Nx] ~ dmnorm(mu[1:Nx], tau[1:Nx,1:Nx]) # multivariate normal prior
  
  # Calculate ORs:
  for (l in 1:Nx) {
    OR[l] <- exp(b[l]);
  }
}

N <- length(IB_2019_Victoria) 
Nx <- 2 

mu <- rep(0, Nx)
tau <- diag(0.001, Nx)
data.logistic <- list("N", "Nx", "IB_2017_Victoria", "IB_2019_Victoria", "mu", "tau")
parameters.logistic <-c("b","OR") # Parameters to keep track of
inits.logistic <- function() {list (b= rep(0, Nx))}

set.seed(114011)
logistic.sim_IB_Victoria_2017_2019_cat1<-jags(data=data.logistic, inits=inits.logistic, parameters.logistic, n.iter=100000,
                                              n.burn=25000, model.file=logistic.model, n.thin=5, n.chains = 3)

print(logistic.sim_IB_Victoria_2017_2019_cat1, 2)
OR_IB_Victoria_2017_2019_cat1_Bayes <- round(logistic.sim_IB_Victoria_2017_2019_cat1$BUGSoutput$summary[2,c(5,3,7)], 2)
OR_IB_Victoria_2017_2019_cat1_Bayes

### sensitivity analysis


try35 <- data_orig %>%
  filter(at2017==1 & at2019==1) %>%
  dplyr::select(Obs, IB_2017_Victoria, IB_2019_Victoria, age_2017_cat)

try35$IB_2017_Victoria <- as.numeric(try35$IB_2017_Victoria)
try35$IB_2019_Victoria <- as.numeric(try35$IB_2019_Victoria)

# 2x2 table
table(try35[try35$age_2017_cat==0,]$IB_2017_Victoria, try35[try35$age_2017_cat==0,]$IB_2019_Victoria) 
table(try35[try35$age_2017_cat==1,]$IB_2017_Victoria, try35[try35$age_2017_cat==1,]$IB_2019_Victoria) 


# Bayesian analysis of H3N2 2014 on H3N2 2016, under 5s
IB_2017_Victoria <- try35[try35$age_2017_cat==0,]$IB_2017_Victoria
IB_2019_Victoria <- try35[try35$age_2017_cat==0,]$IB_2019_Victoria

logistic.model <- function() {
  # SAMPLING DISTRIBUTION
  for (i in 1:N) {
    logit(pi[i]) <- b[1] + b[2]*IB_2017_Victoria[i];
    IB_2019_Victoria[i] ~ dbin(pi[i],1);  
  }
  
  # PRIORS ON BETAS
  b[1:Nx] ~ dmnorm(mu[1:Nx], tau[1:Nx,1:Nx]) # multivariate normal prior
  
  # Calculate ORs:
  for (l in 1:Nx) {
    OR[l] <- exp(b[l]);
  }
}

N <- length(IB_2019_Victoria) 
Nx <- 2 

mu <- rep(0, Nx)
tau <- diag(0.001, Nx)
data.logistic <- list("N", "Nx", "IB_2017_Victoria", "IB_2019_Victoria", "mu", "tau")
parameters.logistic <-c("b","OR") # Parameters to keep track of
inits.logistic <- function() {list (b= rep(0, Nx))}

set.seed(114011)
logistic.sim_IB_Victoria_2017_2019_cat0<-jags(data=data.logistic, inits=inits.logistic, parameters.logistic, n.iter=100000,
                                              n.burn=25000, model.file=logistic.model, n.thin=5, n.chains = 3)

print(logistic.sim_IB_Victoria_2017_2019_cat0, 2)
OR_IB_Victoria_2017_2019_cat0_Bayes <- round(logistic.sim_IB_Victoria_2017_2019_cat0$BUGSoutput$summary[2,c(5,3,7)], 2)
OR_IB_Victoria_2017_2019_cat0_Bayes

# Bayesian analysis of H3N2 2014 on H3N2 2016, over 5s
IB_2017_Victoria <- try35[try35$age_2017_cat==1,]$IB_2017_Victoria
IB_2019_Victoria <- try35[try35$age_2017_cat==1,]$IB_2019_Victoria

logistic.model <- function() {
  # SAMPLING DISTRIBUTION
  for (i in 1:N) {
    logit(pi[i]) <- b[1] + b[2]*IB_2017_Victoria[i];
    IB_2019_Victoria[i] ~ dbin(pi[i],1);  
  }
  
  # PRIORS ON BETAS
  b[1:Nx] ~ dmnorm(mu[1:Nx], tau[1:Nx,1:Nx]) # multivariate normal prior
  
  # Calculate ORs:
  for (l in 1:Nx) {
    OR[l] <- exp(b[l]);
  }
}

N <- length(IB_2019_Victoria) 
Nx <- 2 

mu <- rep(0, Nx)
tau <- diag(0.001, Nx)
data.logistic <- list("N", "Nx", "IB_2017_Victoria", "IB_2019_Victoria", "mu", "tau")
parameters.logistic <-c("b","OR") # Parameters to keep track of
inits.logistic <- function() {list (b= rep(0, Nx))}

set.seed(114011)
logistic.sim_IB_Victoria_2017_2019_cat1<-jags(data=data.logistic, inits=inits.logistic, parameters.logistic, n.iter=100000,
                                              n.burn=25000, model.file=logistic.model, n.thin=5, n.chains = 3)

print(logistic.sim_IB_Victoria_2017_2019_cat1, 2)
OR_IB_Victoria_2017_2019_cat1_Bayes <- round(logistic.sim_IB_Victoria_2017_2019_cat1$BUGSoutput$summary[2,c(5,3,7)], 2)
OR_IB_Victoria_2017_2019_cat1_Bayes


# 2018 to 2019


try35 <- data_orig %>%
  filter(at2018==1 & at2019==1) %>%
  dplyr::select(Obs, IB_2018_Victoria, IB_2019_Victoria, age_2018_cat)

try35$IB_2018_Victoria <- as.numeric(try35$IB_2018_Victoria)
try35$IB_2019_Victoria <- as.numeric(try35$IB_2019_Victoria)

# 2x2 table
table(try35[try35$age_2018_cat==0,]$IB_2018_Victoria, try35[try35$age_2018_cat==0,]$IB_2019_Victoria) 
table(try35[try35$age_2018_cat==1,]$IB_2018_Victoria, try35[try35$age_2018_cat==1,]$IB_2019_Victoria) 


# Bayesian analysis of H3N2 2014 on H3N2 2016, under 5s
IB_2018_Victoria <- try35[try35$age_2018_cat==0,]$IB_2018_Victoria
IB_2019_Victoria <- try35[try35$age_2018_cat==0,]$IB_2019_Victoria

logistic.model <- function() {
  # SAMPLING DISTRIBUTION
  for (i in 1:N) {
    logit(pi[i]) <- b[1] + b[2]*IB_2018_Victoria[i];
    IB_2019_Victoria[i] ~ dbin(pi[i],1);  
  }
  
  # PRIORS ON BETAS
  b[1:Nx] ~ dmnorm(mu[1:Nx], tau[1:Nx,1:Nx]) # multivariate normal prior
  
  # Calculate ORs:
  for (l in 1:Nx) {
    OR[l] <- exp(b[l]);
  }
}

N <- length(IB_2019_Victoria) 
Nx <- 2 

mu <- rep(0, Nx)
tau <- diag(0.001, Nx)
data.logistic <- list("N", "Nx", "IB_2018_Victoria", "IB_2019_Victoria", "mu", "tau")
parameters.logistic <-c("b","OR") # Parameters to keep track of
inits.logistic <- function() {list (b= rep(0, Nx))}

set.seed(114011)
logistic.sim_IB_Victoria_2018_2019_cat0<-jags(data=data.logistic, inits=inits.logistic, parameters.logistic, n.iter=100000,
                                              n.burn=25000, model.file=logistic.model, n.thin=5, n.chains = 3)

print(logistic.sim_IB_Victoria_2018_2019_cat0, 2)
OR_IB_Victoria_2018_2019_cat0_Bayes <- round(logistic.sim_IB_Victoria_2018_2019_cat0$BUGSoutput$summary[2,c(5,3,7)], 2)
OR_IB_Victoria_2018_2019_cat0_Bayes

# Bayesian analysis of H3N2 2014 on H3N2 2016, over 5s
IB_2018_Victoria <- try35[try35$age_2018_cat==1,]$IB_2018_Victoria
IB_2019_Victoria <- try35[try35$age_2018_cat==1,]$IB_2019_Victoria

logistic.model <- function() {
  # SAMPLING DISTRIBUTION
  for (i in 1:N) {
    logit(pi[i]) <- b[1] + b[2]*IB_2018_Victoria[i];
    IB_2019_Victoria[i] ~ dbin(pi[i],1);  
  }
  
  # PRIORS ON BETAS
  b[1:Nx] ~ dmnorm(mu[1:Nx], tau[1:Nx,1:Nx]) # multivariate normal prior
  
  # Calculate ORs:
  for (l in 1:Nx) {
    OR[l] <- exp(b[l]);
  }
}

N <- length(IB_2019_Victoria) 
Nx <- 2 

mu <- rep(0, Nx)
tau <- diag(0.001, Nx)
data.logistic <- list("N", "Nx", "IB_2018_Victoria", "IB_2019_Victoria", "mu", "tau")
parameters.logistic <-c("b","OR") # Parameters to keep track of
inits.logistic <- function() {list (b= rep(0, Nx))}

set.seed(114011)
logistic.sim_IB_Victoria_2018_2019_cat1<-jags(data=data.logistic, inits=inits.logistic, parameters.logistic, n.iter=100000,
                                              n.burn=25000, model.file=logistic.model, n.thin=5, n.chains = 3)

print(logistic.sim_IB_Victoria_2018_2019_cat1, 2)
OR_IB_Victoria_2018_2019_cat1_Bayes <- round(logistic.sim_IB_Victoria_2017_2019_cat1$BUGSoutput$summary[2,c(5,3,7)], 2)
OR_IB_Victoria_2018_2019_cat1_Bayes




