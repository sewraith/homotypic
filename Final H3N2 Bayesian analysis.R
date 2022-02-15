## Analyses for H3N2

library(R2jags)
library(coda)
library(lmtest)
library(foreign)
library(plyr)
library(dplyr)
library(varhandle)

setwd("/Users/sewraith/Box/Research/Homotypic/Data/tempwork/Repeat")


### Analyses of interest: H3N2 ----

# Bring in the data
data_orig <- read.csv2("faustodata_updated.csv", header = T, sep = ",") 

# H3N2 2012 to 2013
try1 <- data_orig %>%
  filter(at2012==1 & at2013==1) %>%
  dplyr::select(Obs, H3N2_2012, H3N2_2013)
try1$H3N2_2012 <- as.numeric(try1$H3N2_2012)
try1$H3N2_2013 <- as.numeric(try1$H3N2_2013)

table(try1$H3N2_2012, try1$H3N2_2013)

# Bayesian analyses
# Putting exposure and outcome variables in the environment 
H3N2_2012 <- try1$H3N2_2012
H3N2_2013 <- try1$H3N2_2013

### Bayesian logistic regression
logistic.model <- function() {
  # SAMPLING DISTRIBUTION
  for (i in 1:N) {
    logit(pi[i]) <- b[1] + b[2]*H3N2_2012[i];
    H3N2_2013[i] ~ dbin(pi[i],1);  
  }
  
  # PRIORS ON BETAS
  b[1:Nx] ~ dmnorm(mu[1:Nx], tau[1:Nx,1:Nx]) # multivariate normal prior
  
  # Calculate ORs:
  for (l in 1:Nx) {
    OR[l] <- exp(b[l]);
  }
}

N <- length(H3N2_2013) # Number of observations
Nx <- 2 # Number of variables: intercept + exposure 

# Data, parameter list and starting values
mu <- rep(0, Nx)
tau <- diag(0.001, Nx)
data.logistic <- list("N", "Nx", "H3N2_2012", "H3N2_2013", "mu", "tau")
parameters.logistic <-c("b","OR") # Parameters to keep track of
inits.logistic <- function() {list (b= rep(0, Nx))}

set.seed(114011)
logistic.sim<-jags(data=data.logistic, inits=inits.logistic, parameters.logistic, n.iter=100000,
                   n.burn=25000, model.file=logistic.model, n.thin=5, n.chains = 3)

print(logistic.sim, 2)
OR_H3N2_2012_2013_Bayes <- round(logistic.sim$BUGSoutput$summary[2,c(5,3,7)], 2)
OR_H3N2_2012_2013_Bayes 

# H3N2 2013 to 2014
try2 <- data_orig %>%
  filter(at2013==1 & at2014==1) %>%
  dplyr::select(Obs, H3N2_2013, H3N2_2014)
try2$H3N2_2013 <- as.numeric(try2$H3N2_2013)
try2$H3N2_2014 <- as.numeric(try2$H3N2_2014)

table(try2$H3N2_2013, try2$H3N2_2014)

# Bayesian analyses
# Putting exposure and outcome variables in the environment 
H3N2_2013 <- try2$H3N2_2013
H3N2_2014 <- try2$H3N2_2014

### Bayesian logistic regression
logistic.model <- function() {
  # SAMPLING DISTRIBUTION
  for (i in 1:N) {
    logit(pi[i]) <- b[1] + b[2]*H3N2_2013[i];
    H3N2_2014[i] ~ dbin(pi[i],1);  
  }
  
  # PRIORS ON BETAS
  b[1:Nx] ~ dmnorm(mu[1:Nx], tau[1:Nx,1:Nx]) # multivariate normal prior
  
  # Calculate ORs:
  for (l in 1:Nx) {
    OR[l] <- exp(b[l]);
  }
}

N <- length(H3N2_2014) # Number of observations
Nx <- 2 # Number of variables: intercept + exposure 

# Data, parameter list and starting values
mu <- rep(0, Nx)
tau <- diag(0.001, Nx)
data.logistic <- list("N", "Nx", "H3N2_2013", "H3N2_2014", "mu", "tau")
parameters.logistic <-c("b","OR") # Parameters to keep track of
inits.logistic <- function() {list (b= rep(0, Nx))}

set.seed(114011)
logistic.sim<-jags(data=data.logistic, inits=inits.logistic, parameters.logistic, n.iter=100000,
                   n.burn=25000, model.file=logistic.model, n.thin=5, n.chains = 3)

print(logistic.sim, 2)
OR_H3N2_2013_2014_Bayes <- round(logistic.sim$BUGSoutput$summary[2,c(5,3,7)], 2)



# H3N2 2016 to 2017
try4 <- data_orig %>%
  filter(at2016==1 & at2017==1) %>%
  dplyr::select(Obs, H3N2_2016, H3N2_2017)
try4$H3N2_2016 <- as.numeric(try4$H3N2_2016)
try4$H3N2_2017 <- as.numeric(try4$H3N2_2017)

table(try4$H3N2_2016, try4$H3N2_2017)


# Bayesian analyses
# Putting exposure and outcome variables in the environment 
H3N2_2016 <- try4$H3N2_2016
H3N2_2017 <- try4$H3N2_2017

### Bayesian logistic regression
logistic.model <- function() {
  # SAMPLING DISTRIBUTION
  for (i in 1:N) {
    logit(pi[i]) <- b[1] + b[2]*H3N2_2016[i];
    H3N2_2017[i] ~ dbin(pi[i],1);  
  }
  
  # PRIORS ON BETAS
  b[1:Nx] ~ dmnorm(mu[1:Nx], tau[1:Nx,1:Nx]) # multivariate normal prior
  
  # Calculate ORs:
  for (l in 1:Nx) {
    OR[l] <- exp(b[l]);
  }
}

N <- length(H3N2_2017) # Number of observations
Nx <- 2 # Number of variables: intercept + exposure 

# Data, parameter list and starting values
mu <- rep(0, Nx)
tau <- diag(0.001, Nx)
data.logistic <- list("N", "Nx", "H3N2_2016", "H3N2_2017", "mu", "tau")
parameters.logistic <-c("b","OR") # Parameters to keep track of
inits.logistic <- function() {list (b= rep(0, Nx))}

set.seed(114011)
logistic.sim<-jags(data=data.logistic, inits=inits.logistic, parameters.logistic, n.iter=100000,
                   n.burn=25000, model.file=logistic.model, n.thin=5, n.chains = 3)

print(logistic.sim, 2)
OR_H3N2_2016_2017_Bayes <- round(logistic.sim$BUGSoutput$summary[2,c(5,3,7)], 2)


# H3N2 2012 to 2014
try51 <- data_orig %>%
  filter(at2012==1 & at2014==1 & H3N2_2013==0) %>%
  dplyr::select(Obs, H3N2_2012, H3N2_2014)
try51$H3N2_2012 <- as.numeric(try51$H3N2_2012)
try51$H3N2_2014 <- as.numeric(try51$H3N2_2014)

table(try51$H3N2_2012, try51$H3N2_2014)

# Bayesian analyses
# Putting exposure and outcome variables in the environment 
H3N2_2012 <- try51$H3N2_2012
H3N2_2014 <- try51$H3N2_2014

### Bayesian logistic regression
logistic.model <- function() {
  # SAMPLING DISTRIBUTION
  for (i in 1:N) {
    logit(pi[i]) <- b[1] + b[2]*H3N2_2012[i];
    H3N2_2014[i] ~ dbin(pi[i],1);  
  }
  
  # PRIORS ON BETAS
  b[1:Nx] ~ dmnorm(mu[1:Nx], tau[1:Nx,1:Nx]) # multivariate normal prior
  
  # Calculate ORs:
  for (l in 1:Nx) {
    OR[l] <- exp(b[l]);
  }
}

N <- length(H3N2_2014) # Number of observations
Nx <- 2 # Number of variables: intercept + exposure 

# Data, parameter list and starting values
mu <- rep(0, Nx)
tau <- diag(0.001, Nx)
data.logistic <- list("N", "Nx", "H3N2_2012", "H3N2_2014", "mu", "tau")
parameters.logistic <-c("b","OR") # Parameters to keep track of
inits.logistic <- function() {list (b= rep(0, Nx))}

set.seed(114011)
logistic.sim<-jags(data=data.logistic, inits=inits.logistic, parameters.logistic, n.iter=100000,
                   n.burn=25000, model.file=logistic.model, n.thin=5, n.chains = 3)

print(logistic.sim, 2)
OR_H3N2_2012_2014_Bayes <- round(logistic.sim$BUGSoutput$summary[2,c(5,3,7)], 2)
OR_H3N2_2012_2014_Bayes

### Sensitivity analysis

try51 <- data_orig %>%
  filter(at2012==1 & at2014==1) %>%
  dplyr::select(Obs, H3N2_2012, H3N2_2014)
try51$H3N2_2012 <- as.numeric(try51$H3N2_2012)
try51$H3N2_2014 <- as.numeric(try51$H3N2_2014)

table(try51$H3N2_2012, try51$H3N2_2014)

# Bayesian analyses
# Putting exposure and outcome variables in the environment 
H3N2_2012 <- try51$H3N2_2012
H3N2_2014 <- try51$H3N2_2014

### Bayesian logistic regression
logistic.model <- function() {
  # SAMPLING DISTRIBUTION
  for (i in 1:N) {
    logit(pi[i]) <- b[1] + b[2]*H3N2_2012[i];
    H3N2_2014[i] ~ dbin(pi[i],1);  
  }
  
  # PRIORS ON BETAS
  b[1:Nx] ~ dmnorm(mu[1:Nx], tau[1:Nx,1:Nx]) # multivariate normal prior
  
  # Calculate ORs:
  for (l in 1:Nx) {
    OR[l] <- exp(b[l]);
  }
}

N <- length(H3N2_2014) # Number of observations
Nx <- 2 # Number of variables: intercept + exposure 

# Data, parameter list and starting values
mu <- rep(0, Nx)
tau <- diag(0.001, Nx)
data.logistic <- list("N", "Nx", "H3N2_2012", "H3N2_2014", "mu", "tau")
parameters.logistic <-c("b","OR") # Parameters to keep track of
inits.logistic <- function() {list (b= rep(0, Nx))}

set.seed(114011)
logistic.sim<-jags(data=data.logistic, inits=inits.logistic, parameters.logistic, n.iter=100000,
                   n.burn=25000, model.file=logistic.model, n.thin=5, n.chains = 3)

print(logistic.sim, 2)
OR_H3N2_2012_2014_Bayes_SA <- round(logistic.sim$BUGSoutput$summary[2,c(5,3,7)], 2)
OR_H3N2_2012_2014_Bayes_SA


# H3N2 2014 to 2016
try51 <- data_orig %>%
  filter(at2014==1 & at2016==1 & H3N2_2015==0) %>%
  dplyr::select(Obs, H3N2_2014, H3N2_2016)
try51$H3N2_2014 <- as.numeric(try51$H3N2_2014)
try51$H3N2_2016 <- as.numeric(try51$H3N2_2016)

table(try51$H3N2_2014, try51$H3N2_2016)

# Bayesian analyses
# Putting exposure and outcome variables in the environment 
H3N2_2014 <- try51$H3N2_2014
H3N2_2016 <- try51$H3N2_2016

### Bayesian logistic regression
logistic.model <- function() {
  # SAMPLING DISTRIBUTION
  for (i in 1:N) {
    logit(pi[i]) <- b[1] + b[2]*H3N2_2014[i];
    H3N2_2016[i] ~ dbin(pi[i],1);  
  }
  
  # PRIORS ON BETAS
  b[1:Nx] ~ dmnorm(mu[1:Nx], tau[1:Nx,1:Nx]) # multivariate normal prior
  
  # Calculate ORs:
  for (l in 1:Nx) {
    OR[l] <- exp(b[l]);
  }
}

N <- length(H3N2_2016) # Number of observations
Nx <- 2 # Number of variables: intercept + exposure 

# Data, parameter list and starting values
mu <- rep(0, Nx)
tau <- diag(0.001, Nx)
data.logistic <- list("N", "Nx", "H3N2_2014", "H3N2_2016", "mu", "tau")
parameters.logistic <-c("b","OR") # Parameters to keep track of
inits.logistic <- function() {list (b= rep(0, Nx))}

set.seed(114011)
logistic.sim<-jags(data=data.logistic, inits=inits.logistic, parameters.logistic, n.iter=100000,
                   n.burn=25000, model.file=logistic.model, n.thin=5, n.chains = 3)

print(logistic.sim, 2)
OR_H3N2_2014_2016_Bayes <- round(logistic.sim$BUGSoutput$summary[2,c(5,3,7)], 2)
OR_H3N2_2014_2016_Bayes

## Sensitivity analysis

try51 <- data_orig %>%
  filter(at2014==1 & at2016==1) %>%
  dplyr::select(Obs, H3N2_2014, H3N2_2016)
try51$H3N2_2014 <- as.numeric(try51$H3N2_2014)
try51$H3N2_2016 <- as.numeric(try51$H3N2_2016)

table(try51$H3N2_2014, try51$H3N2_2016)

# Bayesian analyses
# Putting exposure and outcome variables in the environment 
H3N2_2014 <- try51$H3N2_2014
H3N2_2016 <- try51$H3N2_2016

### Bayesian logistic regression
logistic.model <- function() {
  # SAMPLING DISTRIBUTION
  for (i in 1:N) {
    logit(pi[i]) <- b[1] + b[2]*H3N2_2014[i];
    H3N2_2016[i] ~ dbin(pi[i],1);  
  }
  
  # PRIORS ON BETAS
  b[1:Nx] ~ dmnorm(mu[1:Nx], tau[1:Nx,1:Nx]) # multivariate normal prior
  
  # Calculate ORs:
  for (l in 1:Nx) {
    OR[l] <- exp(b[l]);
  }
}

N <- length(H3N2_2016) # Number of observations
Nx <- 2 # Number of variables: intercept + exposure 

# Data, parameter list and starting values
mu <- rep(0, Nx)
tau <- diag(0.001, Nx)
data.logistic <- list("N", "Nx", "H3N2_2014", "H3N2_2016", "mu", "tau")
parameters.logistic <-c("b","OR") # Parameters to keep track of
inits.logistic <- function() {list (b= rep(0, Nx))}

set.seed(114011)
logistic.sim<-jags(data=data.logistic, inits=inits.logistic, parameters.logistic, n.iter=100000,
                   n.burn=25000, model.file=logistic.model, n.thin=5, n.chains = 3)

print(logistic.sim, 2)
OR_H3N2_2014_2016_Bayes_SA <- round(logistic.sim$BUGSoutput$summary[2,c(5,3,7)], 2)
OR_H3N2_2014_2016_Bayes_SA


# H3N2 2017 to 2019
try5 <- data_orig %>%
  filter(at2017==1 & at2019==1) %>%
  dplyr::select(Obs, H3N2_2017, H3N2_2019)
try5$H3N2_2017 <- as.numeric(try5$H3N2_2017)
try5$H3N2_2019 <- as.numeric(try5$H3N2_2019)

table(try5$H3N2_2017, try5$H3N2_2019)


# Bayesian analyses
# Putting exposure and outcome variables in the environment 
H3N2_2017 <- try5$H3N2_2017
H3N2_2019 <- try5$H3N2_2019

### Bayesian logistic regression
logistic.model <- function() {
  # SAMPLING DISTRIBUTION
  for (i in 1:N) {
    logit(pi[i]) <- b[1] + b[2]*H3N2_2017[i];
    H3N2_2019[i] ~ dbin(pi[i],1);  
  }
  
  # PRIORS ON BETAS
  b[1:Nx] ~ dmnorm(mu[1:Nx], tau[1:Nx,1:Nx]) # multivariate normal prior
  
  # Calculate ORs:
  for (l in 1:Nx) {
    OR[l] <- exp(b[l]);
  }
}

N <- length(H3N2_2019) # Number of observations
Nx <- 2 # Number of variables: intercept + exposure 

# Data, parameter list and starting values
mu <- rep(0, Nx)
tau <- diag(0.001, Nx)
data.logistic <- list("N", "Nx", "H3N2_2017", "H3N2_2019", "mu", "tau")
parameters.logistic <-c("b","OR") # Parameters to keep track of
inits.logistic <- function() {list (b= rep(0, Nx))}

set.seed(114011)
logistic.sim<-jags(data=data.logistic, inits=inits.logistic, parameters.logistic, n.iter=100000,
                   n.burn=25000, model.file=logistic.model, n.thin=5, n.chains = 3)

print(logistic.sim, 2)
OR_H3N2_2017_2019_Bayes <- round(logistic.sim$BUGSoutput$summary[2,c(5,3,7)], 2)



# H3N2 2012 to 2016
try51 <- data_orig %>%
  filter(at2012==1 & at2016==1 & H3N2_2013==0 & H3N2_2014==0 & H3N2_2015==0) %>%
  dplyr::select(Obs, H3N2_2012, H3N2_2016)
try51$H3N2_2012 <- as.numeric(try51$H3N2_2012)
try51$H3N2_2016 <- as.numeric(try51$H3N2_2016)

table(try51$H3N2_2012, try51$H3N2_2016)

# Bayesian analyses
# Putting exposure and outcome variables in the environment 
H3N2_2012 <- try51$H3N2_2012
H3N2_2016 <- try51$H3N2_2016

### Bayesian logistic regression
logistic.model <- function() {
  # SAMPLING DISTRIBUTION
  for (i in 1:N) {
    logit(pi[i]) <- b[1] + b[2]*H3N2_2012[i];
    H3N2_2016[i] ~ dbin(pi[i],1);  
  }
  
  # PRIORS ON BETAS
  b[1:Nx] ~ dmnorm(mu[1:Nx], tau[1:Nx,1:Nx]) # multivariate normal prior
  
  # Calculate ORs:
  for (l in 1:Nx) {
    OR[l] <- exp(b[l]);
  }
}

N <- length(H3N2_2016) # Number of observations
Nx <- 2 # Number of variables: intercept + exposure 

# Data, parameter list and starting values
mu <- rep(0, Nx)
tau <- diag(0.001, Nx)
data.logistic <- list("N", "Nx", "H3N2_2012", "H3N2_2016", "mu", "tau")
parameters.logistic <-c("b","OR") # Parameters to keep track of
inits.logistic <- function() {list (b= rep(0, Nx))}

set.seed(114011)
logistic.sim<-jags(data=data.logistic, inits=inits.logistic, parameters.logistic, n.iter=100000,
                   n.burn=25000, model.file=logistic.model, n.thin=5, n.chains = 3)

print(logistic.sim, 2)
OR_H3N2_2012_2016_Bayes <- round(logistic.sim$BUGSoutput$summary[2,c(5,3,7)], 2)
OR_H3N2_2012_2016_Bayes

### Sensitivity analysis

try51 <- data_orig %>%
  filter(at2012==1 & at2016==1 ) %>%
  dplyr::select(Obs, H3N2_2012, H3N2_2016)
try51$H3N2_2012 <- as.numeric(try51$H3N2_2012)
try51$H3N2_2016 <- as.numeric(try51$H3N2_2016)

table(try51$H3N2_2012, try51$H3N2_2016)

# Bayesian analyses
# Putting exposure and outcome variables in the environment 
H3N2_2012 <- try51$H3N2_2012
H3N2_2016 <- try51$H3N2_2016

### Bayesian logistic regression
logistic.model <- function() {
  # SAMPLING DISTRIBUTION
  for (i in 1:N) {
    logit(pi[i]) <- b[1] + b[2]*H3N2_2012[i];
    H3N2_2016[i] ~ dbin(pi[i],1);  
  }
  
  # PRIORS ON BETAS
  b[1:Nx] ~ dmnorm(mu[1:Nx], tau[1:Nx,1:Nx]) # multivariate normal prior
  
  # Calculate ORs:
  for (l in 1:Nx) {
    OR[l] <- exp(b[l]);
  }
}

N <- length(H3N2_2016) # Number of observations
Nx <- 2 # Number of variables: intercept + exposure 

# Data, parameter list and starting values
mu <- rep(0, Nx)
tau <- diag(0.001, Nx)
data.logistic <- list("N", "Nx", "H3N2_2012", "H3N2_2016", "mu", "tau")
parameters.logistic <-c("b","OR") # Parameters to keep track of
inits.logistic <- function() {list (b= rep(0, Nx))}

set.seed(114011)
logistic.sim<-jags(data=data.logistic, inits=inits.logistic, parameters.logistic, n.iter=100000,
                   n.burn=25000, model.file=logistic.model, n.thin=5, n.chains = 3)

print(logistic.sim, 2)
OR_H3N2_2012_2016_Bayes_SA <- round(logistic.sim$BUGSoutput$summary[2,c(5,3,7)], 2)
OR_H3N2_2012_2016_Bayes_SA


# H3N2 2012 to 2017
try51 <- data_orig %>%
  filter(at2012==1 & at2017==1 & H3N2_2013==0 & H3N2_2014==0 & H3N2_2015==0 & H3N2_2016==0) %>%
  dplyr::select(Obs, H3N2_2012, H3N2_2017)
try51$H3N2_2012 <- as.numeric(try51$H3N2_2012)
try51$H3N2_2017 <- as.numeric(try51$H3N2_2017)

table(try51$H3N2_2012, try51$H3N2_2017)

# Bayesian analyses
# Putting exposure and outcome variables in the environment 
H3N2_2012 <- try51$H3N2_2012
H3N2_2017 <- try51$H3N2_2017

### Bayesian logistic regression
logistic.model <- function() {
  # SAMPLING DISTRIBUTION
  for (i in 1:N) {
    logit(pi[i]) <- b[1] + b[2]*H3N2_2012[i];
    H3N2_2017[i] ~ dbin(pi[i],1);  
  }
  
  # PRIORS ON BETAS
  b[1:Nx] ~ dmnorm(mu[1:Nx], tau[1:Nx,1:Nx]) # multivariate normal prior
  
  # Calculate ORs:
  for (l in 1:Nx) {
    OR[l] <- exp(b[l]);
  }
}

N <- length(H3N2_2017) # Number of observations
Nx <- 2 # Number of variables: intercept + exposure 

# Data, parameter list and starting values
mu <- rep(0, Nx)
tau <- diag(0.001, Nx)
data.logistic <- list("N", "Nx", "H3N2_2012", "H3N2_2017", "mu", "tau")
parameters.logistic <-c("b","OR") # Parameters to keep track of
inits.logistic <- function() {list (b= rep(0, Nx))}

set.seed(114011)
logistic.sim<-jags(data=data.logistic, inits=inits.logistic, parameters.logistic, n.iter=100000,
                   n.burn=25000, model.file=logistic.model, n.thin=5, n.chains = 3)

print(logistic.sim, 2)
OR_H3N2_2012_2017_Bayes <- round(logistic.sim$BUGSoutput$summary[2,c(5,3,7)], 2)
OR_H3N2_2012_2017_Bayes

### Sensitivity analysis

try51 <- data_orig %>%
  filter(at2012==1 & at2017==1 ) %>%
  dplyr::select(Obs, H3N2_2012, H3N2_2017)
try51$H3N2_2012 <- as.numeric(try51$H3N2_2012)
try51$H3N2_2017 <- as.numeric(try51$H3N2_2017)

table(try51$H3N2_2012, try51$H3N2_2017)

# Bayesian analyses
# Putting exposure and outcome variables in the environment 
H3N2_2012 <- try51$H3N2_2012
H3N2_2017 <- try51$H3N2_2017

### Bayesian logistic regression
logistic.model <- function() {
  # SAMPLING DISTRIBUTION
  for (i in 1:N) {
    logit(pi[i]) <- b[1] + b[2]*H3N2_2012[i];
    H3N2_2017[i] ~ dbin(pi[i],1);  
  }
  
  # PRIORS ON BETAS
  b[1:Nx] ~ dmnorm(mu[1:Nx], tau[1:Nx,1:Nx]) # multivariate normal prior
  
  # Calculate ORs:
  for (l in 1:Nx) {
    OR[l] <- exp(b[l]);
  }
}

N <- length(H3N2_2017) # Number of observations
Nx <- 2 # Number of variables: intercept + exposure 

# Data, parameter list and starting values
mu <- rep(0, Nx)
tau <- diag(0.001, Nx)
data.logistic <- list("N", "Nx", "H3N2_2012", "H3N2_2017", "mu", "tau")
parameters.logistic <-c("b","OR") # Parameters to keep track of
inits.logistic <- function() {list (b= rep(0, Nx))}

set.seed(114011)
logistic.sim<-jags(data=data.logistic, inits=inits.logistic, parameters.logistic, n.iter=100000,
                   n.burn=25000, model.file=logistic.model, n.thin=5, n.chains = 3)

print(logistic.sim, 2)
OR_H3N2_2012_2017_Bayes_SA <- round(logistic.sim$BUGSoutput$summary[2,c(5,3,7)], 2)
OR_H3N2_2012_2017_Bayes_SA

# H3N2 2012 to 2019
try51 <- data_orig %>%
  filter(at2012==1 & at2019==1 & H3N2_2013==0 & H3N2_2014==0 & H3N2_2015==0 & H3N2_2016==0 & H3N2_2017==0) %>%
  dplyr::select(Obs, H3N2_2012, H3N2_2019)
try51$H3N2_2012 <- as.numeric(try51$H3N2_2012)
try51$H3N2_2019 <- as.numeric(try51$H3N2_2019)

table(try51$H3N2_2012, try51$H3N2_2019)

# Bayesian analyses
# Putting exposure and outcome variables in the environment 
H3N2_2012 <- try51$H3N2_2012
H3N2_2019 <- try51$H3N2_2019

### Bayesian logistic regression
logistic.model <- function() {
  # SAMPLING DISTRIBUTION
  for (i in 1:N) {
    logit(pi[i]) <- b[1] + b[2]*H3N2_2012[i];
    H3N2_2019[i] ~ dbin(pi[i],1);  
  }
  
  # PRIORS ON BETAS
  b[1:Nx] ~ dmnorm(mu[1:Nx], tau[1:Nx,1:Nx]) # multivariate normal prior
  
  # Calculate ORs:
  for (l in 1:Nx) {
    OR[l] <- exp(b[l]);
  }
}

N <- length(H3N2_2019) # Number of observations
Nx <- 2 # Number of variables: intercept + exposure 

# Data, parameter list and starting values
mu <- rep(0, Nx)
tau <- diag(0.001, Nx)
data.logistic <- list("N", "Nx", "H3N2_2012", "H3N2_2019", "mu", "tau")
parameters.logistic <-c("b","OR") # Parameters to keep track of
inits.logistic <- function() {list (b= rep(0, Nx))}

set.seed(114011)
logistic.sim<-jags(data=data.logistic, inits=inits.logistic, parameters.logistic, n.iter=100000,
                   n.burn=25000, model.file=logistic.model, n.thin=5, n.chains = 3)

print(logistic.sim, 2)
OR_H3N2_2012_2019_Bayes <- round(logistic.sim$BUGSoutput$summary[2,c(5,3,7)], 2)
OR_H3N2_2012_2019_Bayes

### Sensitivity analysis

try51 <- data_orig %>%
  filter(at2012==1 & at2019==1 ) %>%
  dplyr::select(Obs, H3N2_2012, H3N2_2019)
try51$H3N2_2012 <- as.numeric(try51$H3N2_2012)
try51$H3N2_2019 <- as.numeric(try51$H3N2_2019)

table(try51$H3N2_2012, try51$H3N2_2019)

# Bayesian analyses
# Putting exposure and outcome variables in the environment 
H3N2_2012 <- try51$H3N2_2012
H3N2_2019 <- try51$H3N2_2019

### Bayesian logistic regression
logistic.model <- function() {
  # SAMPLING DISTRIBUTION
  for (i in 1:N) {
    logit(pi[i]) <- b[1] + b[2]*H3N2_2012[i];
    H3N2_2019[i] ~ dbin(pi[i],1);  
  }
  
  # PRIORS ON BETAS
  b[1:Nx] ~ dmnorm(mu[1:Nx], tau[1:Nx,1:Nx]) # multivariate normal prior
  
  # Calculate ORs:
  for (l in 1:Nx) {
    OR[l] <- exp(b[l]);
  }
}

N <- length(H3N2_2019) # Number of observations
Nx <- 2 # Number of variables: intercept + exposure 

# Data, parameter list and starting values
mu <- rep(0, Nx)
tau <- diag(0.001, Nx)
data.logistic <- list("N", "Nx", "H3N2_2012", "H3N2_2019", "mu", "tau")
parameters.logistic <-c("b","OR") # Parameters to keep track of
inits.logistic <- function() {list (b= rep(0, Nx))}

set.seed(114011)
logistic.sim<-jags(data=data.logistic, inits=inits.logistic, parameters.logistic, n.iter=100000,
                   n.burn=25000, model.file=logistic.model, n.thin=5, n.chains = 3)

print(logistic.sim, 2)
OR_H3N2_2012_2019_Bayes_SA <- round(logistic.sim$BUGSoutput$summary[2,c(5,3,7)], 2)
OR_H3N2_2012_2019_Bayes_SA


# H3N2 2013 to 2016
try51 <- data_orig %>%
  filter(at2013==1 & at2016==1 & H3N2_2014==0 & H3N2_2015==0 ) %>%
  dplyr::select(Obs, H3N2_2013, H3N2_2016)
try51$H3N2_2013 <- as.numeric(try51$H3N2_2013)
try51$H3N2_2016 <- as.numeric(try51$H3N2_2016)

table(try51$H3N2_2013, try51$H3N2_2016)

# Bayesian analyses
# Putting exposure and outcome variables in the environment 
H3N2_2013 <- try51$H3N2_2013
H3N2_2016 <- try51$H3N2_2016

### Bayesian logistic regression
logistic.model <- function() {
  # SAMPLING DISTRIBUTION
  for (i in 1:N) {
    logit(pi[i]) <- b[1] + b[2]*H3N2_2013[i];
    H3N2_2016[i] ~ dbin(pi[i],1);  
  }
  
  # PRIORS ON BETAS
  b[1:Nx] ~ dmnorm(mu[1:Nx], tau[1:Nx,1:Nx]) # multivariate normal prior
  
  # Calculate ORs:
  for (l in 1:Nx) {
    OR[l] <- exp(b[l]);
  }
}

N <- length(H3N2_2016) # Number of observations
Nx <- 2 # Number of variables: intercept + exposure 

# Data, parameter list and starting values
mu <- rep(0, Nx)
tau <- diag(0.001, Nx)
data.logistic <- list("N", "Nx", "H3N2_2013", "H3N2_2016", "mu", "tau")
parameters.logistic <-c("b","OR") # Parameters to keep track of
inits.logistic <- function() {list (b= rep(0, Nx))}

set.seed(114011)
logistic.sim<-jags(data=data.logistic, inits=inits.logistic, parameters.logistic, n.iter=100000,
                   n.burn=25000, model.file=logistic.model, n.thin=5, n.chains = 3)

print(logistic.sim, 2)
OR_H3N2_2013_2016_Bayes <- round(logistic.sim$BUGSoutput$summary[2,c(5,3,7)], 2)
OR_H3N2_2013_2016_Bayes

### Sensitivity Analysis

try51 <- data_orig %>%
  filter(at2013==1 & at2016==1) %>%
  dplyr::select(Obs, H3N2_2013, H3N2_2016)
try51$H3N2_2013 <- as.numeric(try51$H3N2_2013)
try51$H3N2_2016 <- as.numeric(try51$H3N2_2016)

table(try51$H3N2_2013, try51$H3N2_2016)

# Bayesian analyses
# Putting exposure and outcome variables in the environment 
H3N2_2013 <- try51$H3N2_2013
H3N2_2016 <- try51$H3N2_2016

### Bayesian logistic regression
logistic.model <- function() {
  # SAMPLING DISTRIBUTION
  for (i in 1:N) {
    logit(pi[i]) <- b[1] + b[2]*H3N2_2013[i];
    H3N2_2016[i] ~ dbin(pi[i],1);  
  }
  
  # PRIORS ON BETAS
  b[1:Nx] ~ dmnorm(mu[1:Nx], tau[1:Nx,1:Nx]) # multivariate normal prior
  
  # Calculate ORs:
  for (l in 1:Nx) {
    OR[l] <- exp(b[l]);
  }
}

N <- length(H3N2_2016) # Number of observations
Nx <- 2 # Number of variables: intercept + exposure 

# Data, parameter list and starting values
mu <- rep(0, Nx)
tau <- diag(0.001, Nx)
data.logistic <- list("N", "Nx", "H3N2_2013", "H3N2_2016", "mu", "tau")
parameters.logistic <-c("b","OR") # Parameters to keep track of
inits.logistic <- function() {list (b= rep(0, Nx))}

set.seed(114011)
logistic.sim<-jags(data=data.logistic, inits=inits.logistic, parameters.logistic, n.iter=100000,
                   n.burn=25000, model.file=logistic.model, n.thin=5, n.chains = 3)

print(logistic.sim, 2)
OR_H3N2_2013_2016_Bayes_SA <- round(logistic.sim$BUGSoutput$summary[2,c(5,3,7)], 2)
OR_H3N2_2013_2016_Bayes_SA




# H3N2 2013 to 2017
try51 <- data_orig %>%
  filter(at2013==1 & at2017==1 & H3N2_2014==0 & H3N2_2015==0 & H3N2_2016==0) %>%
  dplyr::select(Obs, H3N2_2013, H3N2_2017)
try51$H3N2_2013 <- as.numeric(try51$H3N2_2013)
try51$H3N2_2017 <- as.numeric(try51$H3N2_2017)

table(try51$H3N2_2013, try51$H3N2_2017)

# Bayesian analyses
# Putting exposure and outcome variables in the environment 
H3N2_2013 <- try51$H3N2_2013
H3N2_2017 <- try51$H3N2_2017

### Bayesian logistic regression
logistic.model <- function() {
  # SAMPLING DISTRIBUTION
  for (i in 1:N) {
    logit(pi[i]) <- b[1] + b[2]*H3N2_2013[i];
    H3N2_2017[i] ~ dbin(pi[i],1);  
  }
  
  # PRIORS ON BETAS
  b[1:Nx] ~ dmnorm(mu[1:Nx], tau[1:Nx,1:Nx]) # multivariate normal prior
  
  # Calculate ORs:
  for (l in 1:Nx) {
    OR[l] <- exp(b[l]);
  }
}

N <- length(H3N2_2017) # Number of observations
Nx <- 2 # Number of variables: intercept + exposure 

# Data, parameter list and starting values
mu <- rep(0, Nx)
tau <- diag(0.001, Nx)
data.logistic <- list("N", "Nx", "H3N2_2013", "H3N2_2017", "mu", "tau")
parameters.logistic <-c("b","OR") # Parameters to keep track of
inits.logistic <- function() {list (b= rep(0, Nx))}

set.seed(114011)
logistic.sim<-jags(data=data.logistic, inits=inits.logistic, parameters.logistic, n.iter=100000,
                   n.burn=25000, model.file=logistic.model, n.thin=5, n.chains = 3)

print(logistic.sim, 2)
OR_H3N2_2013_2017_Bayes <- round(logistic.sim$BUGSoutput$summary[2,c(5,3,7)], 2)
OR_H3N2_2013_2017_Bayes

### Sensitivity analysis

try51 <- data_orig %>%
  filter(at2013==1 & at2017==1) %>%
  dplyr::select(Obs, H3N2_2013, H3N2_2017)
try51$H3N2_2013 <- as.numeric(try51$H3N2_2013)
try51$H3N2_2017 <- as.numeric(try51$H3N2_2017)

table(try51$H3N2_2013, try51$H3N2_2017)

# Bayesian analyses
# Putting exposure and outcome variables in the environment 
H3N2_2013 <- try51$H3N2_2013
H3N2_2017 <- try51$H3N2_2017

### Bayesian logistic regression
logistic.model <- function() {
  # SAMPLING DISTRIBUTION
  for (i in 1:N) {
    logit(pi[i]) <- b[1] + b[2]*H3N2_2013[i];
    H3N2_2017[i] ~ dbin(pi[i],1);  
  }
  
  # PRIORS ON BETAS
  b[1:Nx] ~ dmnorm(mu[1:Nx], tau[1:Nx,1:Nx]) # multivariate normal prior
  
  # Calculate ORs:
  for (l in 1:Nx) {
    OR[l] <- exp(b[l]);
  }
}

N <- length(H3N2_2017) # Number of observations
Nx <- 2 # Number of variables: intercept + exposure 

# Data, parameter list and starting values
mu <- rep(0, Nx)
tau <- diag(0.001, Nx)
data.logistic <- list("N", "Nx", "H3N2_2013", "H3N2_2017", "mu", "tau")
parameters.logistic <-c("b","OR") # Parameters to keep track of
inits.logistic <- function() {list (b= rep(0, Nx))}

set.seed(114011)
logistic.sim<-jags(data=data.logistic, inits=inits.logistic, parameters.logistic, n.iter=100000,
                   n.burn=25000, model.file=logistic.model, n.thin=5, n.chains = 3)

print(logistic.sim, 2)
OR_H3N2_2013_2017_Bayes_SA <- round(logistic.sim$BUGSoutput$summary[2,c(5,3,7)], 2)
OR_H3N2_2013_2017_Bayes_SA



# H3N2 2013 to 2019
try51 <- data_orig %>%
  filter(at2013==1 & at2019==1 & H3N2_2014==0 & H3N2_2015==0 & H3N2_2016==0 & H3N2_2017==0) %>%
  dplyr::select(Obs, H3N2_2013, H3N2_2019)
try51$H3N2_2013 <- as.numeric(try51$H3N2_2013)
try51$H3N2_2019 <- as.numeric(try51$H3N2_2019)

table(try51$H3N2_2013, try51$H3N2_2019)

# Bayesian analyses
# Putting exposure and outcome variables in the environment 
H3N2_2013 <- try51$H3N2_2013
H3N2_2019 <- try51$H3N2_2019

### Bayesian logistic regression
logistic.model <- function() {
  # SAMPLING DISTRIBUTION
  for (i in 1:N) {
    logit(pi[i]) <- b[1] + b[2]*H3N2_2013[i];
    H3N2_2019[i] ~ dbin(pi[i],1);  
  }
  
  # PRIORS ON BETAS
  b[1:Nx] ~ dmnorm(mu[1:Nx], tau[1:Nx,1:Nx]) # multivariate normal prior
  
  # Calculate ORs:
  for (l in 1:Nx) {
    OR[l] <- exp(b[l]);
  }
}

N <- length(H3N2_2019) # Number of observations
Nx <- 2 # Number of variables: intercept + exposure 

# Data, parameter list and starting values
mu <- rep(0, Nx)
tau <- diag(0.001, Nx)
data.logistic <- list("N", "Nx", "H3N2_2013", "H3N2_2019", "mu", "tau")
parameters.logistic <-c("b","OR") # Parameters to keep track of
inits.logistic <- function() {list (b= rep(0, Nx))}

set.seed(114011)
logistic.sim<-jags(data=data.logistic, inits=inits.logistic, parameters.logistic, n.iter=100000,
                   n.burn=25000, model.file=logistic.model, n.thin=5, n.chains = 3)

print(logistic.sim, 2)
OR_H3N2_2013_2019_Bayes <- round(logistic.sim$BUGSoutput$summary[2,c(5,3,7)], 2)
OR_H3N2_2013_2019_Bayes

### Sensitivity analysis

try51 <- data_orig %>%
  filter(at2013==1 & at2019==1 ) %>%
  dplyr::select(Obs, H3N2_2013, H3N2_2019)
try51$H3N2_2013 <- as.numeric(try51$H3N2_2013)
try51$H3N2_2019 <- as.numeric(try51$H3N2_2019)

table(try51$H3N2_2013, try51$H3N2_2019)

# Bayesian analyses
# Putting exposure and outcome variables in the environment 
H3N2_2013 <- try51$H3N2_2013
H3N2_2019 <- try51$H3N2_2019

### Bayesian logistic regression
logistic.model <- function() {
  # SAMPLING DISTRIBUTION
  for (i in 1:N) {
    logit(pi[i]) <- b[1] + b[2]*H3N2_2013[i];
    H3N2_2019[i] ~ dbin(pi[i],1);  
  }
  
  # PRIORS ON BETAS
  b[1:Nx] ~ dmnorm(mu[1:Nx], tau[1:Nx,1:Nx]) # multivariate normal prior
  
  # Calculate ORs:
  for (l in 1:Nx) {
    OR[l] <- exp(b[l]);
  }
}

N <- length(H3N2_2019) # Number of observations
Nx <- 2 # Number of variables: intercept + exposure 

# Data, parameter list and starting values
mu <- rep(0, Nx)
tau <- diag(0.001, Nx)
data.logistic <- list("N", "Nx", "H3N2_2013", "H3N2_2019", "mu", "tau")
parameters.logistic <-c("b","OR") # Parameters to keep track of
inits.logistic <- function() {list (b= rep(0, Nx))}

set.seed(114011)
logistic.sim<-jags(data=data.logistic, inits=inits.logistic, parameters.logistic, n.iter=100000,
                   n.burn=25000, model.file=logistic.model, n.thin=5, n.chains = 3)

print(logistic.sim, 2)
OR_H3N2_2013_2019_Bayes_SA <- round(logistic.sim$BUGSoutput$summary[2,c(5,3,7)], 2)
OR_H3N2_2013_2019_Bayes_SA


# H3N2 2014 to 2017
try51 <- data_orig %>%
  filter(at2014==1 & at2017==1 & H3N2_2015==0 & H3N2_2016==0 ) %>%
  dplyr::select(Obs, H3N2_2014, H3N2_2017)
try51$H3N2_2014 <- as.numeric(try51$H3N2_2014)
try51$H3N2_2017 <- as.numeric(try51$H3N2_2017)

table(try51$H3N2_2014, try51$H3N2_2017)

# Bayesian analyses
# Putting exposure and outcome variables in the environment 
H3N2_2014 <- try51$H3N2_2014
H3N2_2017 <- try51$H3N2_2017

### Bayesian logistic regression
logistic.model <- function() {
  # SAMPLING DISTRIBUTION
  for (i in 1:N) {
    logit(pi[i]) <- b[1] + b[2]*H3N2_2014[i];
    H3N2_2017[i] ~ dbin(pi[i],1);  
  }
  
  # PRIORS ON BETAS
  b[1:Nx] ~ dmnorm(mu[1:Nx], tau[1:Nx,1:Nx]) # multivariate normal prior
  
  # Calculate ORs:
  for (l in 1:Nx) {
    OR[l] <- exp(b[l]);
  }
}

N <- length(H3N2_2017) # Number of observations
Nx <- 2 # Number of variables: intercept + exposure 

# Data, parameter list and starting values
mu <- rep(0, Nx)
tau <- diag(0.001, Nx)
data.logistic <- list("N", "Nx", "H3N2_2014", "H3N2_2017", "mu", "tau")
parameters.logistic <-c("b","OR") # Parameters to keep track of
inits.logistic <- function() {list (b= rep(0, Nx))}

set.seed(114011)
logistic.sim<-jags(data=data.logistic, inits=inits.logistic, parameters.logistic, n.iter=100000,
                   n.burn=25000, model.file=logistic.model, n.thin=5, n.chains = 3)

print(logistic.sim, 2)
OR_H3N2_2014_2017_Bayes <- round(logistic.sim$BUGSoutput$summary[2,c(5,3,7)], 2)
OR_H3N2_2014_2017_Bayes

### Sensitivity Analysis

try51 <- data_orig %>%
  filter(at2014==1 & at2017==1  ) %>%
  dplyr::select(Obs, H3N2_2014, H3N2_2017)
try51$H3N2_2014 <- as.numeric(try51$H3N2_2014)
try51$H3N2_2017 <- as.numeric(try51$H3N2_2017)

table(try51$H3N2_2014, try51$H3N2_2017)

# Bayesian analyses
# Putting exposure and outcome variables in the environment 
H3N2_2014 <- try51$H3N2_2014
H3N2_2017 <- try51$H3N2_2017

### Bayesian logistic regression
logistic.model <- function() {
  # SAMPLING DISTRIBUTION
  for (i in 1:N) {
    logit(pi[i]) <- b[1] + b[2]*H3N2_2014[i];
    H3N2_2017[i] ~ dbin(pi[i],1);  
  }
  
  # PRIORS ON BETAS
  b[1:Nx] ~ dmnorm(mu[1:Nx], tau[1:Nx,1:Nx]) # multivariate normal prior
  
  # Calculate ORs:
  for (l in 1:Nx) {
    OR[l] <- exp(b[l]);
  }
}

N <- length(H3N2_2017) # Number of observations
Nx <- 2 # Number of variables: intercept + exposure 

# Data, parameter list and starting values
mu <- rep(0, Nx)
tau <- diag(0.001, Nx)
data.logistic <- list("N", "Nx", "H3N2_2014", "H3N2_2017", "mu", "tau")
parameters.logistic <-c("b","OR") # Parameters to keep track of
inits.logistic <- function() {list (b= rep(0, Nx))}

set.seed(114011)
logistic.sim<-jags(data=data.logistic, inits=inits.logistic, parameters.logistic, n.iter=100000,
                   n.burn=25000, model.file=logistic.model, n.thin=5, n.chains = 3)

print(logistic.sim, 2)
OR_H3N2_2014_2017_Bayes_SA <- round(logistic.sim$BUGSoutput$summary[2,c(5,3,7)], 2)
OR_H3N2_2014_2017_Bayes_SA


# H3N2 2014 to 2019
try51 <- data_orig %>%
  filter(at2014==1 & at2019==1 & H3N2_2015==0 & H3N2_2016==0 & H3N2_2017==0) %>%
  dplyr::select(Obs, H3N2_2014, H3N2_2019)
try51$H3N2_2014 <- as.numeric(try51$H3N2_2014)
try51$H3N2_2019 <- as.numeric(try51$H3N2_2019)

table(try51$H3N2_2014, try51$H3N2_2019)

# Bayesian analyses
# Putting exposure and outcome variables in the environment 
H3N2_2014 <- try51$H3N2_2014
H3N2_2019 <- try51$H3N2_2019

### Bayesian logistic regression
logistic.model <- function() {
  # SAMPLING DISTRIBUTION
  for (i in 1:N) {
    logit(pi[i]) <- b[1] + b[2]*H3N2_2014[i];
    H3N2_2019[i] ~ dbin(pi[i],1);  
  }
  
  # PRIORS ON BETAS
  b[1:Nx] ~ dmnorm(mu[1:Nx], tau[1:Nx,1:Nx]) # multivariate normal prior
  
  # Calculate ORs:
  for (l in 1:Nx) {
    OR[l] <- exp(b[l]);
  }
}

N <- length(H3N2_2019) # Number of observations
Nx <- 2 # Number of variables: intercept + exposure 

# Data, parameter list and starting values
mu <- rep(0, Nx)
tau <- diag(0.001, Nx)
data.logistic <- list("N", "Nx", "H3N2_2014", "H3N2_2019", "mu", "tau")
parameters.logistic <-c("b","OR") # Parameters to keep track of
inits.logistic <- function() {list (b= rep(0, Nx))}

set.seed(114011)
logistic.sim<-jags(data=data.logistic, inits=inits.logistic, parameters.logistic, n.iter=100000,
                   n.burn=25000, model.file=logistic.model, n.thin=5, n.chains = 3)

print(logistic.sim, 2)
OR_H3N2_2014_2019_Bayes <- round(logistic.sim$BUGSoutput$summary[2,c(5,3,7)], 2)
OR_H3N2_2014_2019_Bayes

### Sensitivity analysis

try51 <- data_orig %>%
  filter(at2014==1 & at2019==1 ) %>%
  dplyr::select(Obs, H3N2_2014, H3N2_2019)
try51$H3N2_2014 <- as.numeric(try51$H3N2_2014)
try51$H3N2_2019 <- as.numeric(try51$H3N2_2019)

table(try51$H3N2_2014, try51$H3N2_2019)

# Bayesian analyses
# Putting exposure and outcome variables in the environment 
H3N2_2014 <- try51$H3N2_2014
H3N2_2019 <- try51$H3N2_2019

### Bayesian logistic regression
logistic.model <- function() {
  # SAMPLING DISTRIBUTION
  for (i in 1:N) {
    logit(pi[i]) <- b[1] + b[2]*H3N2_2014[i];
    H3N2_2019[i] ~ dbin(pi[i],1);  
  }
  
  # PRIORS ON BETAS
  b[1:Nx] ~ dmnorm(mu[1:Nx], tau[1:Nx,1:Nx]) # multivariate normal prior
  
  # Calculate ORs:
  for (l in 1:Nx) {
    OR[l] <- exp(b[l]);
  }
}

N <- length(H3N2_2019) # Number of observations
Nx <- 2 # Number of variables: intercept + exposure 

# Data, parameter list and starting values
mu <- rep(0, Nx)
tau <- diag(0.001, Nx)
data.logistic <- list("N", "Nx", "H3N2_2014", "H3N2_2019", "mu", "tau")
parameters.logistic <-c("b","OR") # Parameters to keep track of
inits.logistic <- function() {list (b= rep(0, Nx))}

set.seed(114011)
logistic.sim<-jags(data=data.logistic, inits=inits.logistic, parameters.logistic, n.iter=100000,
                   n.burn=25000, model.file=logistic.model, n.thin=5, n.chains = 3)

print(logistic.sim, 2)
OR_H3N2_2014_2019_Bayes_SA <- round(logistic.sim$BUGSoutput$summary[2,c(5,3,7)], 2)
OR_H3N2_2014_2019_Bayes_SA


# H3N2 2016 to 2019
try51 <- data_orig %>%
  filter(at2016==1 & at2019==1 & H3N2_2017==0) %>%
  dplyr::select(Obs, H3N2_2016, H3N2_2019)
try51$H3N2_2016 <- as.numeric(try51$H3N2_2016)
try51$H3N2_2019 <- as.numeric(try51$H3N2_2019)

table(try51$H3N2_2016, try51$H3N2_2019)

# Bayesian analyses
# Putting exposure and outcome variables in the environment 
H3N2_2016 <- try51$H3N2_2016
H3N2_2019 <- try51$H3N2_2019

### Bayesian logistic regression
logistic.model <- function() {
  # SAMPLING DISTRIBUTION
  for (i in 1:N) {
    logit(pi[i]) <- b[1] + b[2]*H3N2_2016[i];
    H3N2_2019[i] ~ dbin(pi[i],1);  
  }
  
  # PRIORS ON BETAS
  b[1:Nx] ~ dmnorm(mu[1:Nx], tau[1:Nx,1:Nx]) # multivariate normal prior
  
  # Calculate ORs:
  for (l in 1:Nx) {
    OR[l] <- exp(b[l]);
  }
}

N <- length(H3N2_2019) # Number of observations
Nx <- 2 # Number of variables: intercept + exposure 

# Data, parameter list and starting values
mu <- rep(0, Nx)
tau <- diag(0.001, Nx)
data.logistic <- list("N", "Nx", "H3N2_2016", "H3N2_2019", "mu", "tau")
parameters.logistic <-c("b","OR") # Parameters to keep track of
inits.logistic <- function() {list (b= rep(0, Nx))}

set.seed(114011)
logistic.sim<-jags(data=data.logistic, inits=inits.logistic, parameters.logistic, n.iter=100000,
                   n.burn=25000, model.file=logistic.model, n.thin=5, n.chains = 3)

print(logistic.sim, 2)
OR_H3N2_2016_2019_Bayes <- round(logistic.sim$BUGSoutput$summary[2,c(5,3,7)], 2)
OR_H3N2_2016_2019_Bayes


### Sensitivity Analysis ###

try51 <- data_orig %>%
  filter(at2016==1 & at2019==1 ) %>%
  dplyr::select(Obs, H3N2_2016, H3N2_2019)
try51$H3N2_2016 <- as.numeric(try51$H3N2_2016)
try51$H3N2_2019 <- as.numeric(try51$H3N2_2019)

table(try51$H3N2_2016, try51$H3N2_2019)

# Bayesian analyses
# Putting exposure and outcome variables in the environment 
H3N2_2016 <- try51$H3N2_2016
H3N2_2019 <- try51$H3N2_2019

### Bayesian logistic regression
logistic.model <- function() {
  # SAMPLING DISTRIBUTION
  for (i in 1:N) {
    logit(pi[i]) <- b[1] + b[2]*H3N2_2016[i];
    H3N2_2019[i] ~ dbin(pi[i],1);  
  }
  
  # PRIORS ON BETAS
  b[1:Nx] ~ dmnorm(mu[1:Nx], tau[1:Nx,1:Nx]) # multivariate normal prior
  
  # Calculate ORs:
  for (l in 1:Nx) {
    OR[l] <- exp(b[l]);
  }
}

N <- length(H3N2_2019) # Number of observations
Nx <- 2 # Number of variables: intercept + exposure 

# Data, parameter list and starting values
mu <- rep(0, Nx)
tau <- diag(0.001, Nx)
data.logistic <- list("N", "Nx", "H3N2_2016", "H3N2_2019", "mu", "tau")
parameters.logistic <-c("b","OR") # Parameters to keep track of
inits.logistic <- function() {list (b= rep(0, Nx))}

set.seed(114011)
logistic.sim<-jags(data=data.logistic, inits=inits.logistic, parameters.logistic, n.iter=100000,
                   n.burn=25000, model.file=logistic.model, n.thin=5, n.chains = 3)

print(logistic.sim, 2)
OR_H3N2_2016_2019_Bayes_SA <- round(logistic.sim$BUGSoutput$summary[2,c(5,3,7)], 2)
OR_H3N2_2016_2019_Bayes_SA


### Age-stratified analyses
try5 <- data_orig %>%
  filter(at2012==1 & at2013==1) %>%
  dplyr::select(Obs, H3N2_2012, H3N2_2013, age_2012_cat)

try5$H3N2_2012 <- as.numeric(try5$H3N2_2012)
try5$H3N2_2013 <- as.numeric(try5$H3N2_2013)

# 2x2 table
table(try5[try5$age_2012_cat==0,]$H3N2_2012, try5[try5$age_2012_cat==0,]$H3N2_2013) 
table(try5[try5$age_2012_cat==1,]$H3N2_2012, try5[try5$age_2012_cat==1,]$H3N2_2013) 


# Bayesian analysis of H3N2 2012 on H3N2 2013, under 5s
H3N2_2012 <- try5[try5$age_2012_cat==0,]$H3N2_2012
H3N2_2013 <- try5[try5$age_2012_cat==0,]$H3N2_2013

logistic.model <- function() {
  # SAMPLING DISTRIBUTION
  for (i in 1:N) {
    logit(pi[i]) <- b[1] + b[2]*H3N2_2012[i];
    H3N2_2013[i] ~ dbin(pi[i],1);  
  }
  
  # PRIORS ON BETAS
  b[1:Nx] ~ dmnorm(mu[1:Nx], tau[1:Nx,1:Nx]) # multivariate normal prior
  
  # Calculate ORs:
  for (l in 1:Nx) {
    OR[l] <- exp(b[l]);
  }
}

N <- length(H3N2_2013) 
Nx <- 2 

mu <- rep(0, Nx)
tau <- diag(0.001, Nx)
data.logistic <- list("N", "Nx", "H3N2_2012", "H3N2_2013", "mu", "tau")
parameters.logistic <-c("b","OR") # Parameters to keep track of
inits.logistic <- function() {list (b= rep(0, Nx))}

set.seed(114011)
logistic.sim_2012_2013_cat0<-jags(data=data.logistic, inits=inits.logistic, parameters.logistic, n.iter=100000,
                                  n.burn=25000, model.file=logistic.model, n.thin=5, n.chains = 3)

print(logistic.sim_2012_2013_cat0, 2)
OR_H3N2_2012_2013_cat0_Bayes <- round(logistic.sim_2012_2013_cat0$BUGSoutput$summary[2,c(5,3,7)], 2)
OR_H3N2_2012_2013_cat0_Bayes


# Bayesian analysis of H3N2 2012 on H3N2 2013, over 5s
H3N2_2012 <- try5[try5$age_2012_cat==1,]$H3N2_2012
H3N2_2013 <- try5[try5$age_2012_cat==1,]$H3N2_2013

logistic.model <- function() {
  # SAMPLING DISTRIBUTION
  for (i in 1:N) {
    logit(pi[i]) <- b[1] + b[2]*H3N2_2012[i];
    H3N2_2013[i] ~ dbin(pi[i],1);  
  }
  
  # PRIORS ON BETAS
  b[1:Nx] ~ dmnorm(mu[1:Nx], tau[1:Nx,1:Nx]) # multivariate normal prior
  
  # Calculate ORs:
  for (l in 1:Nx) {
    OR[l] <- exp(b[l]);
  }
}

N <- length(H3N2_2013) 
Nx <- 2 

mu <- rep(0, Nx)
tau <- diag(0.001, Nx)
data.logistic <- list("N", "Nx", "H3N2_2012", "H3N2_2013", "mu", "tau")
parameters.logistic <-c("b","OR") # Parameters to keep track of
inits.logistic <- function() {list (b= rep(0, Nx))}

set.seed(114011)
logistic.sim_2012_2013_cat1<-jags(data=data.logistic, inits=inits.logistic, parameters.logistic, n.iter=100000,
                                  n.burn=25000, model.file=logistic.model, n.thin=5, n.chains = 3)

print(logistic.sim_2012_2013_cat1, 2)
OR_H3N2_2012_2013_cat1 <- round(logistic.sim_2012_2013_cat1$BUGSoutput$summary[2,c(5,3,7)], 2)
OR_H3N2_2012_2013_cat1

# no sensitivity analysis needed

try6 <- data_orig %>%
  filter(at2013==1 & at2014==1) %>%
  dplyr::select(Obs, H3N2_2013, H3N2_2014, age_2013_cat)

try6$H3N2_2013 <- as.numeric(try6$H3N2_2013)
try6$H3N2_2014 <- as.numeric(try6$H3N2_2014)

# 2x2 table
table(try6[try6$age_2013_cat==0,]$H3N2_2013, try6[try6$age_2013_cat==0,]$H3N2_2014) # Numbers match Steph's
table(try6[try6$age_2013_cat==1,]$H3N2_2013, try6[try6$age_2013_cat==1,]$H3N2_2014) # Numbers match Steph's

# Freq analysis of H3N2 2013 on H3N2 2013, under 5s
freq.fit <- glm(H3N2_2014 ~ H3N2_2013, family=binomial, data=try6[try6$age_2013_cat==0,])
summary(freq.fit) 
OR.freq <- exp(cbind(coef(freq.fit), confint(freq.fit)))
colnames(OR.freq) <- c("OR","95% LL","95% UL")
round(OR.freq, digits=2) 
OR_H3N2_2013_2014_cat0_freq <- round(OR.freq, digits=2)[2,]

# Freq analysis of H3N2 2013 on H3N2 2014, over 5s
freq.fit <- glm(H3N2_2014 ~ H3N2_2013, family=binomial, data=try6[try6$age_2013_cat==1,])
summary(freq.fit) 
OR.freq <- exp(cbind(coef(freq.fit), confint(freq.fit)))
colnames(OR.freq) <- c("OR","95% LL","95% UL")
round(OR.freq, digits=2) 
OR_H3N2_2013_2014_cat1_freq <- round(OR.freq, digits=2)[2,]

# Bayesian analysis of H3N2 2013 on H3N2 2014, under 5s
H3N2_2013 <- try6[try6$age_2013_cat==0,]$H3N2_2013
H3N2_2014 <- try6[try6$age_2013_cat==0,]$H3N2_2014

logistic.model <- function() {
  # SAMPLING DISTRIBUTION
  for (i in 1:N) {
    logit(pi[i]) <- b[1] + b[2]*H3N2_2013[i];
    H3N2_2014[i] ~ dbin(pi[i],1);  
  }
  
  # PRIORS ON BETAS
  b[1:Nx] ~ dmnorm(mu[1:Nx], tau[1:Nx,1:Nx]) # multivariate normal prior
  
  # Calculate ORs:
  for (l in 1:Nx) {
    OR[l] <- exp(b[l]);
  }
}

N <- length(H3N2_2013) 
Nx <- 2 

mu <- rep(0, Nx)
tau <- diag(0.001, Nx)
data.logistic <- list("N", "Nx", "H3N2_2013", "H3N2_2014", "mu", "tau")
parameters.logistic <-c("b","OR") # Parameters to keep track of
inits.logistic <- function() {list (b= rep(0, Nx))}

set.seed(114011)
logistic.sim_2013_2014_cat0<-jags(data=data.logistic, inits=inits.logistic, parameters.logistic, n.iter=100000,
                                  n.burn=25000, model.file=logistic.model, n.thin=5, n.chains = 3)

print(logistic.sim_2013_2014_cat0, 2)
OR_H3N2_2013_2014_cat0_Bayes <- round(logistic.sim_2013_2014_cat0$BUGSoutput$summary[2,c(5,3,7)], 2)

OR_H3N2_2013_2014_cat0_freq
OR_H3N2_2013_2014_cat0_Bayes

# Bayesian analysis of H3N2 2013 on H3N2 2014, over 5s
H3N2_2013 <- try6[try6$age_2013_cat==1,]$H3N2_2013
H3N2_2014 <- try6[try6$age_2013_cat==1,]$H3N2_2014

logistic.model <- function() {
  # SAMPLING DISTRIBUTION
  for (i in 1:N) {
    logit(pi[i]) <- b[1] + b[2]*H3N2_2013[i];
    H3N2_2014[i] ~ dbin(pi[i],1);  
  }
  
  # PRIORS ON BETAS
  b[1:Nx] ~ dmnorm(mu[1:Nx], tau[1:Nx,1:Nx]) # multivariate normal prior
  
  # Calculate ORs:
  for (l in 1:Nx) {
    OR[l] <- exp(b[l]);
  }
}

N <- length(H3N2_2013) 
Nx <- 2 

mu <- rep(0, Nx)
tau <- diag(0.001, Nx)
data.logistic <- list("N", "Nx", "H3N2_2013", "H3N2_2014", "mu", "tau")
parameters.logistic <-c("b","OR") # Parameters to keep track of
inits.logistic <- function() {list (b= rep(0, Nx))}

set.seed(114011)
logistic.sim_2013_2014_cat1<-jags(data=data.logistic, inits=inits.logistic, parameters.logistic, n.iter=100000,
                                  n.burn=25000, model.file=logistic.model, n.thin=5, n.chains = 3)

print(logistic.sim_2013_2014_cat1, 2)
OR_H3N2_2013_2014_cat1_Bayes <- round(logistic.sim_2013_2014_cat1$BUGSoutput$summary[2,c(5,3,7)], 2)

OR_H3N2_2013_2014_cat1_freq
OR_H3N2_2013_2014_cat1_Bayes

# no sensitivity analysis needed

# 2016 to 2017

### Analyses of interest: Age-stratified H3N2 2016 to 2017 ----
try20 <- data_orig %>%
  filter(at2016==1 & at2017==1) %>%
  dplyr::select(Obs, H3N2_2016, H3N2_2017, age_2016_cat)

try20$H3N2_2016 <- as.numeric(try20$H3N2_2016)
try20$H3N2_2017 <- as.numeric(try20$H3N2_2017)

# 2x2 table
table(try20[try20$age_2016_cat==0,]$H3N2_2016, try20[try20$age_2016_cat==0,]$H3N2_2017) 
table(try20[try20$age_2016_cat==1,]$H3N2_2016, try20[try20$age_2016_cat==1,]$H3N2_2017) 


# Freq analysis of H3N2 2016 on H3N2 2017, under 5s
freq.fit <- glm(H3N2_2017 ~ H3N2_2016, family=binomial, data=try20[try20$age_2016_cat==0,])
summary(freq.fit) 
OR.freq <- exp(cbind(coef(freq.fit), confint(freq.fit)))
colnames(OR.freq) <- c("OR","95% LL","95% UL")
round(OR.freq, digits=2) 
OR_H3N2_2016_2017_cat0_freq <- round(OR.freq, digits=2)[2,]

# Freq analysis of H3N2 2016 on H3N2 2017, over 5s
freq.fit <- glm(H3N2_2017 ~ H3N2_2016, family=binomial, data=try20[try20$age_2016_cat==1,])
summary(freq.fit) 
OR.freq <- exp(cbind(coef(freq.fit), confint(freq.fit)))
colnames(OR.freq) <- c("OR","95% LL","95% UL")
round(OR.freq, digits=2) 
OR_H3N2_2016_2017_cat1_freq <- round(OR.freq, digits=2)[2,]

# Bayesian analysis of H3N2 2014 on H3N2 2016, under 5s
H3N2_2016 <- try20[try20$age_2016_cat==0,]$H3N2_2016
H3N2_2017 <- try20[try20$age_2016_cat==0,]$H3N2_2017

logistic.model <- function() {
  # SAMPLING DISTRIBUTION
  for (i in 1:N) {
    logit(pi[i]) <- b[1] + b[2]*H3N2_2016[i];
    H3N2_2017[i] ~ dbin(pi[i],1);  
  }
  
  # PRIORS ON BETAS
  b[1:Nx] ~ dmnorm(mu[1:Nx], tau[1:Nx,1:Nx]) # multivariate normal prior
  
  # Calculate ORs:
  for (l in 1:Nx) {
    OR[l] <- exp(b[l]);
  }
}

N <- length(H3N2_2017) 
Nx <- 2 

mu <- rep(0, Nx)
tau <- diag(0.001, Nx)
data.logistic <- list("N", "Nx", "H3N2_2016", "H3N2_2017", "mu", "tau")
parameters.logistic <-c("b","OR") # Parameters to keep track of
inits.logistic <- function() {list (b= rep(0, Nx))}

set.seed(114011)
logistic.sim_2016_2017_cat0<-jags(data=data.logistic, inits=inits.logistic, parameters.logistic, n.iter=100000,
                                  n.burn=25000, model.file=logistic.model, n.thin=5, n.chains = 3)

print(logistic.sim_2016_2017_cat0, 2)
OR_H3N2_2016_2017_cat0_Bayes <- round(logistic.sim_2016_2017_cat0$BUGSoutput$summary[2,c(5,3,7)], 2)

OR_H3N2_2016_2017_cat0_freq
OR_H3N2_2016_2017_cat0_Bayes

# Bayesian analysis of H3N2 2014 on H3N2 2016, over 5s
H3N2_2016 <- try20[try20$age_2016_cat==1,]$H3N2_2016
H3N2_2017 <- try20[try20$age_2016_cat==1,]$H3N2_2017

logistic.model <- function() {
  # SAMPLING DISTRIBUTION
  for (i in 1:N) {
    logit(pi[i]) <- b[1] + b[2]*H3N2_2016[i];
    H3N2_2017[i] ~ dbin(pi[i],1);  
  }
  
  # PRIORS ON BETAS
  b[1:Nx] ~ dmnorm(mu[1:Nx], tau[1:Nx,1:Nx]) # multivariate normal prior
  
  # Calculate ORs:
  for (l in 1:Nx) {
    OR[l] <- exp(b[l]);
  }
}

N <- length(H3N2_2017) 
Nx <- 2 

mu <- rep(0, Nx)
tau <- diag(0.001, Nx)
data.logistic <- list("N", "Nx", "H3N2_2016", "H3N2_2017", "mu", "tau")
parameters.logistic <-c("b","OR") # Parameters to keep track of
inits.logistic <- function() {list (b= rep(0, Nx))}

set.seed(114011)
logistic.sim_2016_2017_cat1<-jags(data=data.logistic, inits=inits.logistic, parameters.logistic, n.iter=100000,
                                  n.burn=25000, model.file=logistic.model, n.thin=5, n.chains = 3)

print(logistic.sim_2016_2017_cat1, 2)
OR_H3N2_2016_2017_cat1_Bayes <- round(logistic.sim_2016_2017_cat1$BUGSoutput$summary[2,c(5,3,7)], 2)

OR_H3N2_2016_2017_cat1_freq
OR_H3N2_2016_2017_cat1_Bayes
 # no sensitivity analysis needed


# 2012 to 2014

try7 <- data_orig %>%
  filter(at2012==1 & at2014==1 & H3N2_2013==0) %>%
  dplyr::select(Obs, H3N2_2012, H3N2_2014, age_2012_cat)

try7$H3N2_2012 <- as.numeric(try7$H3N2_2012)
try7$H3N2_2014 <- as.numeric(try7$H3N2_2014)

# 2x2 table
table(try7[try7$age_2012_cat==0,]$H3N2_2012, try7[try7$age_2012_cat==0,]$H3N2_2014) 
table(try7[try7$age_2012_cat==1,]$H3N2_2012, try7[try7$age_2012_cat==1,]$H3N2_2014) 

# Bayesian analysis of H3N2 2014 on H3N2 2016, under 5s
H3N2_2012 <- try7[try7$age_2012_cat==0,]$H3N2_2012
H3N2_2014 <- try7[try7$age_2012_cat==0,]$H3N2_2014

logistic.model <- function() {
  # SAMPLING DISTRIBUTION
  for (i in 1:N) {
    logit(pi[i]) <- b[1] + b[2]*H3N2_2012[i];
    H3N2_2014[i] ~ dbin(pi[i],1);  
  }
  
  # PRIORS ON BETAS
  b[1:Nx] ~ dmnorm(mu[1:Nx], tau[1:Nx,1:Nx]) # multivariate normal prior
  
  # Calculate ORs:
  for (l in 1:Nx) {
    OR[l] <- exp(b[l]);
  }
}

N <- length(H3N2_2014) 
Nx <- 2 

mu <- rep(0, Nx)
tau <- diag(0.001, Nx)
data.logistic <- list("N", "Nx", "H3N2_2012", "H3N2_2014", "mu", "tau")
parameters.logistic <-c("b","OR") # Parameters to keep track of
inits.logistic <- function() {list (b= rep(0, Nx))}

set.seed(114011)
logistic.sim_2012_2014_cat0<-jags(data=data.logistic, inits=inits.logistic, parameters.logistic, n.iter=100000,
                                  n.burn=25000, model.file=logistic.model, n.thin=5, n.chains = 3)

print(logistic.sim_2012_2014_cat0, 2)
OR_H3N2_2012_2014_cat0_Bayes <- round(logistic.sim_2012_2014_cat0$BUGSoutput$summary[2,c(5,3,7)], 2)
OR_H3N2_2012_2014_cat0_Bayes

# Bayesian analysis of H3N2 2014 on H3N2 2016, over 5s
H3N2_2012 <- try7[try7$age_2012_cat==1,]$H3N2_2012
H3N2_2014 <- try7[try7$age_2012_cat==1,]$H3N2_2014

logistic.model <- function() {
  # SAMPLING DISTRIBUTION
  for (i in 1:N) {
    logit(pi[i]) <- b[1] + b[2]*H3N2_2012[i];
    H3N2_2014[i] ~ dbin(pi[i],1);  
  }
  
  # PRIORS ON BETAS
  b[1:Nx] ~ dmnorm(mu[1:Nx], tau[1:Nx,1:Nx]) # multivariate normal prior
  
  # Calculate ORs:
  for (l in 1:Nx) {
    OR[l] <- exp(b[l]);
  }
}

N <- length(H3N2_2014) 
Nx <- 2 

mu <- rep(0, Nx)
tau <- diag(0.001, Nx)
data.logistic <- list("N", "Nx", "H3N2_2012", "H3N2_2014", "mu", "tau")
parameters.logistic <-c("b","OR") # Parameters to keep track of
inits.logistic <- function() {list (b= rep(0, Nx))}

set.seed(114011)
logistic.sim_2012_2014_cat1<-jags(data=data.logistic, inits=inits.logistic, parameters.logistic, n.iter=100000,
                                  n.burn=25000, model.file=logistic.model, n.thin=5, n.chains = 3)

print(logistic.sim_2012_2014_cat1, 2)
OR_H3N2_2012_2014_cat1_Bayes <- round(logistic.sim_2012_2014_cat1$BUGSoutput$summary[2,c(5,3,7)], 2)
OR_H3N2_2012_2014_cat1_Bayes

### sensitivity analysis

try7 <- data_orig %>%
  filter(at2012==1 & at2014==1) %>%
  dplyr::select(Obs, H3N2_2012, H3N2_2014, age_2012_cat)

try7$H3N2_2012 <- as.numeric(try7$H3N2_2012)
try7$H3N2_2014 <- as.numeric(try7$H3N2_2014)

# 2x2 table
table(try7[try7$age_2012_cat==0,]$H3N2_2012, try7[try7$age_2012_cat==0,]$H3N2_2014) 
table(try7[try7$age_2012_cat==1,]$H3N2_2012, try7[try7$age_2012_cat==1,]$H3N2_2014) 

# Bayesian analysis of H3N2 2014 on H3N2 2016, under 5s
H3N2_2012 <- try7[try7$age_2012_cat==0,]$H3N2_2012
H3N2_2014 <- try7[try7$age_2012_cat==0,]$H3N2_2014

logistic.model <- function() {
  # SAMPLING DISTRIBUTION
  for (i in 1:N) {
    logit(pi[i]) <- b[1] + b[2]*H3N2_2012[i];
    H3N2_2014[i] ~ dbin(pi[i],1);  
  }
  
  # PRIORS ON BETAS
  b[1:Nx] ~ dmnorm(mu[1:Nx], tau[1:Nx,1:Nx]) # multivariate normal prior
  
  # Calculate ORs:
  for (l in 1:Nx) {
    OR[l] <- exp(b[l]);
  }
}

N <- length(H3N2_2014) 
Nx <- 2 

mu <- rep(0, Nx)
tau <- diag(0.001, Nx)
data.logistic <- list("N", "Nx", "H3N2_2012", "H3N2_2014", "mu", "tau")
parameters.logistic <-c("b","OR") # Parameters to keep track of
inits.logistic <- function() {list (b= rep(0, Nx))}

set.seed(114011)
logistic.sim_2012_2014_cat0<-jags(data=data.logistic, inits=inits.logistic, parameters.logistic, n.iter=100000,
                                  n.burn=25000, model.file=logistic.model, n.thin=5, n.chains = 3)

print(logistic.sim_2012_2014_cat0, 2)
OR_H3N2_2012_2014_cat0_Bayes <- round(logistic.sim_2012_2014_cat0$BUGSoutput$summary[2,c(5,3,7)], 2)
OR_H3N2_2012_2014_cat0_Bayes

# Bayesian analysis of H3N2 2014 on H3N2 2016, over 5s
H3N2_2012 <- try7[try7$age_2012_cat==1,]$H3N2_2012
H3N2_2014 <- try7[try7$age_2012_cat==1,]$H3N2_2014

logistic.model <- function() {
  # SAMPLING DISTRIBUTION
  for (i in 1:N) {
    logit(pi[i]) <- b[1] + b[2]*H3N2_2012[i];
    H3N2_2014[i] ~ dbin(pi[i],1);  
  }
  
  # PRIORS ON BETAS
  b[1:Nx] ~ dmnorm(mu[1:Nx], tau[1:Nx,1:Nx]) # multivariate normal prior
  
  # Calculate ORs:
  for (l in 1:Nx) {
    OR[l] <- exp(b[l]);
  }
}

N <- length(H3N2_2014) 
Nx <- 2 

mu <- rep(0, Nx)
tau <- diag(0.001, Nx)
data.logistic <- list("N", "Nx", "H3N2_2012", "H3N2_2014", "mu", "tau")
parameters.logistic <-c("b","OR") # Parameters to keep track of
inits.logistic <- function() {list (b= rep(0, Nx))}

set.seed(114011)
logistic.sim_2012_2014_cat1<-jags(data=data.logistic, inits=inits.logistic, parameters.logistic, n.iter=100000,
                                  n.burn=25000, model.file=logistic.model, n.thin=5, n.chains = 3)

print(logistic.sim_2012_2014_cat1, 2)
OR_H3N2_2012_2014_cat1_Bayes <- round(logistic.sim_2012_2014_cat1$BUGSoutput$summary[2,c(5,3,7)], 2)
OR_H3N2_2012_2014_cat1_Bayes

# 2014 to 2016
try7 <- data_orig %>%
  filter(at2014==1 & at2016==1 & H3N2_2015==0) %>%
  dplyr::select(Obs, H3N2_2014, H3N2_2016, age_2014_cat)

try7$H3N2_2014 <- as.numeric(try7$H3N2_2014)
try7$H3N2_2016 <- as.numeric(try7$H3N2_2016)

# 2x2 table
table(try7[try7$age_2014_cat==0,]$H3N2_2014, try7[try7$age_2014_cat==0,]$H3N2_2016) 
table(try7[try7$age_2014_cat==1,]$H3N2_2014, try7[try7$age_2014_cat==1,]$H3N2_2016) 

# Bayesian analysis of H3N2 2014 on H3N2 2016, under 5s
H3N2_2014 <- try7[try7$age_2014_cat==0,]$H3N2_2014
H3N2_2016 <- try7[try7$age_2014_cat==0,]$H3N2_2016

logistic.model <- function() {
  # SAMPLING DISTRIBUTION
  for (i in 1:N) {
    logit(pi[i]) <- b[1] + b[2]*H3N2_2014[i];
    H3N2_2016[i] ~ dbin(pi[i],1);  
  }
  
  # PRIORS ON BETAS
  b[1:Nx] ~ dmnorm(mu[1:Nx], tau[1:Nx,1:Nx]) # multivariate normal prior
  
  # Calculate ORs:
  for (l in 1:Nx) {
    OR[l] <- exp(b[l]);
  }
}

N <- length(H3N2_2016) 
Nx <- 2 

mu <- rep(0, Nx)
tau <- diag(0.001, Nx)
data.logistic <- list("N", "Nx", "H3N2_2014", "H3N2_2016", "mu", "tau")
parameters.logistic <-c("b","OR") # Parameters to keep track of
inits.logistic <- function() {list (b= rep(0, Nx))}

set.seed(114011)
logistic.sim_2014_2016_cat0<-jags(data=data.logistic, inits=inits.logistic, parameters.logistic, n.iter=100000,
                                  n.burn=25000, model.file=logistic.model, n.thin=5, n.chains = 3)

print(logistic.sim_2014_2016_cat0, 2)
OR_H3N2_2014_2016_cat0_Bayes <- round(logistic.sim_2014_2016_cat0$BUGSoutput$summary[2,c(5,3,7)], 2)
OR_H3N2_2014_2016_cat0_Bayes

# Bayesian analysis of H3N2 2014 on H3N2 2016, over 5s
H3N2_2014 <- try7[try7$age_2014_cat==1,]$H3N2_2014
H3N2_2016 <- try7[try7$age_2014_cat==1,]$H3N2_2016

logistic.model <- function() {
  # SAMPLING DISTRIBUTION
  for (i in 1:N) {
    logit(pi[i]) <- b[1] + b[2]*H3N2_2014[i];
    H3N2_2016[i] ~ dbin(pi[i],1);  
  }
  
  # PRIORS ON BETAS
  b[1:Nx] ~ dmnorm(mu[1:Nx], tau[1:Nx,1:Nx]) # multivariate normal prior
  
  # Calculate ORs:
  for (l in 1:Nx) {
    OR[l] <- exp(b[l]);
  }
}

N <- length(H3N2_2016) 
Nx <- 2 

mu <- rep(0, Nx)
tau <- diag(0.001, Nx)
data.logistic <- list("N", "Nx", "H3N2_2014", "H3N2_2016", "mu", "tau")
parameters.logistic <-c("b","OR") # Parameters to keep track of
inits.logistic <- function() {list (b= rep(0, Nx))}

set.seed(114011)
logistic.sim_2014_2016_cat1<-jags(data=data.logistic, inits=inits.logistic, parameters.logistic, n.iter=100000,
                                  n.burn=25000, model.file=logistic.model, n.thin=5, n.chains = 3)

print(logistic.sim_2014_2016_cat1, 2)
OR_H3N2_2014_2016_cat1_Bayes <- round(logistic.sim_2014_2016_cat1$BUGSoutput$summary[2,c(5,3,7)], 2)
OR_H3N2_2014_2016_cat1_Bayes

### sensitivity analysis

try7 <- data_orig %>%
  filter(at2014==1 & at2016==1) %>%
  dplyr::select(Obs, H3N2_2014, H3N2_2016, age_2014_cat)

try7$H3N2_2014 <- as.numeric(try7$H3N2_2014)
try7$H3N2_2016 <- as.numeric(try7$H3N2_2016)

# 2x2 table
table(try7[try7$age_2014_cat==0,]$H3N2_2014, try7[try7$age_2014_cat==0,]$H3N2_2016) 
table(try7[try7$age_2014_cat==1,]$H3N2_2014, try7[try7$age_2014_cat==1,]$H3N2_2016) 

# Bayesian analysis of H3N2 2014 on H3N2 2016, under 5s
H3N2_2014 <- try7[try7$age_2014_cat==0,]$H3N2_2014
H3N2_2016 <- try7[try7$age_2014_cat==0,]$H3N2_2016

logistic.model <- function() {
  # SAMPLING DISTRIBUTION
  for (i in 1:N) {
    logit(pi[i]) <- b[1] + b[2]*H3N2_2014[i];
    H3N2_2016[i] ~ dbin(pi[i],1);  
  }
  
  # PRIORS ON BETAS
  b[1:Nx] ~ dmnorm(mu[1:Nx], tau[1:Nx,1:Nx]) # multivariate normal prior
  
  # Calculate ORs:
  for (l in 1:Nx) {
    OR[l] <- exp(b[l]);
  }
}

N <- length(H3N2_2016) 
Nx <- 2 

mu <- rep(0, Nx)
tau <- diag(0.001, Nx)
data.logistic <- list("N", "Nx", "H3N2_2014", "H3N2_2016", "mu", "tau")
parameters.logistic <-c("b","OR") # Parameters to keep track of
inits.logistic <- function() {list (b= rep(0, Nx))}

set.seed(114011)
logistic.sim_2014_2016_cat0<-jags(data=data.logistic, inits=inits.logistic, parameters.logistic, n.iter=100000,
                                  n.burn=25000, model.file=logistic.model, n.thin=5, n.chains = 3)

print(logistic.sim_2014_2016_cat0, 2)
OR_H3N2_2014_2016_cat0_Bayes <- round(logistic.sim_2014_2016_cat0$BUGSoutput$summary[2,c(5,3,7)], 2)
OR_H3N2_2014_2016_cat0_Bayes

# Bayesian analysis of H3N2 2014 on H3N2 2016, over 5s
H3N2_2014 <- try7[try7$age_2014_cat==1,]$H3N2_2014
H3N2_2016 <- try7[try7$age_2014_cat==1,]$H3N2_2016

logistic.model <- function() {
  # SAMPLING DISTRIBUTION
  for (i in 1:N) {
    logit(pi[i]) <- b[1] + b[2]*H3N2_2014[i];
    H3N2_2016[i] ~ dbin(pi[i],1);  
  }
  
  # PRIORS ON BETAS
  b[1:Nx] ~ dmnorm(mu[1:Nx], tau[1:Nx,1:Nx]) # multivariate normal prior
  
  # Calculate ORs:
  for (l in 1:Nx) {
    OR[l] <- exp(b[l]);
  }
}

N <- length(H3N2_2016) 
Nx <- 2 

mu <- rep(0, Nx)
tau <- diag(0.001, Nx)
data.logistic <- list("N", "Nx", "H3N2_2014", "H3N2_2016", "mu", "tau")
parameters.logistic <-c("b","OR") # Parameters to keep track of
inits.logistic <- function() {list (b= rep(0, Nx))}

set.seed(114011)
logistic.sim_2014_2016_cat1<-jags(data=data.logistic, inits=inits.logistic, parameters.logistic, n.iter=100000,
                                  n.burn=25000, model.file=logistic.model, n.thin=5, n.chains = 3)

print(logistic.sim_2014_2016_cat1, 2)
OR_H3N2_2014_2016_cat1_Bayes <- round(logistic.sim_2014_2016_cat1$BUGSoutput$summary[2,c(5,3,7)], 2)
OR_H3N2_2014_2016_cat1_Bayes



### Analyses of interest: Age-stratified H3N2 2017 to 2019 ----
try21 <- data_orig %>%
  filter(at2017==1 & at2019==1) %>%
  dplyr::select(Obs, H3N2_2017, H3N2_2019, age_2017_cat)

try21$H3N2_2017 <- as.numeric(try21$H3N2_2017)
try21$H3N2_2019 <- as.numeric(try21$H3N2_2019)

# 2x2 table
table(try21[try21$age_2017_cat==0,]$H3N2_2017, try21[try21$age_2017_cat==0,]$H3N2_2019) 
table(try21[try21$age_2017_cat==1,]$H3N2_2017, try21[try21$age_2017_cat==1,]$H3N2_2019) 


# Freq analysis of H3N2 2016 on H3N2 2017, under 5s
freq.fit <- glm(H3N2_2019 ~ H3N2_2017, family=binomial, data=try21[try21$age_2017_cat==0,])
summary(freq.fit) 
OR.freq <- exp(cbind(coef(freq.fit), confint(freq.fit)))
colnames(OR.freq) <- c("OR","95% LL","95% UL")
round(OR.freq, digits=2) 
OR_H3N2_2017_2019_cat0_freq <- round(OR.freq, digits=2)[2,]

# Freq analysis of H3N2 2016 on H3N2 2017, over 5s
freq.fit <- glm(H3N2_2019 ~ H3N2_2017, family=binomial, data=try21[try21$age_2017_cat==1,])
summary(freq.fit) 
OR.freq <- exp(cbind(coef(freq.fit), confint(freq.fit)))
colnames(OR.freq) <- c("OR","95% LL","95% UL")
round(OR.freq, digits=2) 
OR_H3N2_2017_2019_cat1_freq <- round(OR.freq, digits=2)[2,]

# Bayesian analysis of H3N2 2014 on H3N2 2016, under 5s
H3N2_2017 <- try21[try21$age_2017_cat==0,]$H3N2_2017
H3N2_2019 <- try21[try21$age_2017_cat==0,]$H3N2_2019

logistic.model <- function() {
  # SAMPLING DISTRIBUTION
  for (i in 1:N) {
    logit(pi[i]) <- b[1] + b[2]*H3N2_2017[i];
    H3N2_2019[i] ~ dbin(pi[i],1);  
  }
  
  # PRIORS ON BETAS
  b[1:Nx] ~ dmnorm(mu[1:Nx], tau[1:Nx,1:Nx]) # multivariate normal prior
  
  # Calculate ORs:
  for (l in 1:Nx) {
    OR[l] <- exp(b[l]);
  }
}

N <- length(H3N2_2019) 
Nx <- 2 

mu <- rep(0, Nx)
tau <- diag(0.001, Nx)
data.logistic <- list("N", "Nx", "H3N2_2017", "H3N2_2019", "mu", "tau")
parameters.logistic <-c("b","OR") # Parameters to keep track of
inits.logistic <- function() {list (b= rep(0, Nx))}

set.seed(114011)
logistic.sim_2017_2019_cat0<-jags(data=data.logistic, inits=inits.logistic, parameters.logistic, n.iter=100000,
                                  n.burn=25000, model.file=logistic.model, n.thin=5, n.chains = 3)

print(logistic.sim_2017_2019_cat0, 2)
OR_H3N2_2017_2019_cat0_Bayes <- round(logistic.sim_2017_2019_cat0$BUGSoutput$summary[2,c(5,3,7)], 2)

OR_H3N2_2017_2019_cat0_freq
OR_H3N2_2017_2019_cat0_Bayes

# Bayesian analysis of H3N2 2014 on H3N2 2016, over 5s
H3N2_2017 <- try21[try21$age_2017_cat==1,]$H3N2_2017
H3N2_2019 <- try21[try21$age_2017_cat==1,]$H3N2_2019

logistic.model <- function() {
  # SAMPLING DISTRIBUTION
  for (i in 1:N) {
    logit(pi[i]) <- b[1] + b[2]*H3N2_2017[i];
    H3N2_2019[i] ~ dbin(pi[i],1);  
  }
  
  # PRIORS ON BETAS
  b[1:Nx] ~ dmnorm(mu[1:Nx], tau[1:Nx,1:Nx]) # multivariate normal prior
  
  # Calculate ORs:
  for (l in 1:Nx) {
    OR[l] <- exp(b[l]);
  }
}

N <- length(H3N2_2019) 
Nx <- 2 

mu <- rep(0, Nx)
tau <- diag(0.001, Nx)
data.logistic <- list("N", "Nx", "H3N2_2017", "H3N2_2019", "mu", "tau")
parameters.logistic <-c("b","OR") # Parameters to keep track of
inits.logistic <- function() {list (b= rep(0, Nx))}

set.seed(114011)
logistic.sim_2017_2019_cat1<-jags(data=data.logistic, inits=inits.logistic, parameters.logistic, n.iter=100000,
                                  n.burn=25000, model.file=logistic.model, n.thin=5, n.chains = 3)

print(logistic.sim_2017_2019_cat1, 2)
OR_H3N2_2017_2019_cat1_Bayes <- round(logistic.sim_2017_2019_cat1$BUGSoutput$summary[2,c(5,3,7)], 2)

OR_H3N2_2017_2019_cat1_freq
OR_H3N2_2017_2019_cat1_Bayes

# no sensitivity analysis needed


# 2012 to 2016

try7 <- data_orig %>%
  filter(at2012==1 & at2016==1 & H3N2_2013==0 & H3N2_2014==0 & H3N2_2015==0) %>%
  dplyr::select(Obs, H3N2_2012, H3N2_2016, age_2012_cat)

try7$H3N2_2012 <- as.numeric(try7$H3N2_2012)
try7$H3N2_2016 <- as.numeric(try7$H3N2_2016)

# 2x2 table
table(try7[try7$age_2012_cat==0,]$H3N2_2012, try7[try7$age_2012_cat==0,]$H3N2_2016) 
table(try7[try7$age_2012_cat==1,]$H3N2_2012, try7[try7$age_2012_cat==1,]$H3N2_2016) 

# Bayesian analysis of H3N2 2014 on H3N2 2016, under 5s
H3N2_2012 <- try7[try7$age_2012_cat==0,]$H3N2_2012
H3N2_2016 <- try7[try7$age_2012_cat==0,]$H3N2_2016

logistic.model <- function() {
  # SAMPLING DISTRIBUTION
  for (i in 1:N) {
    logit(pi[i]) <- b[1] + b[2]*H3N2_2012[i];
    H3N2_2016[i] ~ dbin(pi[i],1);  
  }
  
  # PRIORS ON BETAS
  b[1:Nx] ~ dmnorm(mu[1:Nx], tau[1:Nx,1:Nx]) # multivariate normal prior
  
  # Calculate ORs:
  for (l in 1:Nx) {
    OR[l] <- exp(b[l]);
  }
}

N <- length(H3N2_2016) 
Nx <- 2 

mu <- rep(0, Nx)
tau <- diag(0.001, Nx)
data.logistic <- list("N", "Nx", "H3N2_2012", "H3N2_2016", "mu", "tau")
parameters.logistic <-c("b","OR") # Parameters to keep track of
inits.logistic <- function() {list (b= rep(0, Nx))}

set.seed(114011)
logistic.sim_2012_2016_cat0<-jags(data=data.logistic, inits=inits.logistic, parameters.logistic, n.iter=100000,
                                  n.burn=25000, model.file=logistic.model, n.thin=5, n.chains = 3)

print(logistic.sim_2012_2016_cat0, 2)
OR_H3N2_2012_2016_cat0_Bayes <- round(logistic.sim_2012_2016_cat0$BUGSoutput$summary[2,c(5,3,7)], 2)
OR_H3N2_2012_2016_cat0_Bayes

# Bayesian analysis of H3N2 2014 on H3N2 2016, over 5s
H3N2_2012 <- try7[try7$age_2012_cat==1,]$H3N2_2012
H3N2_2016 <- try7[try7$age_2012_cat==1,]$H3N2_2016

logistic.model <- function() {
  # SAMPLING DISTRIBUTION
  for (i in 1:N) {
    logit(pi[i]) <- b[1] + b[2]*H3N2_2012[i];
    H3N2_2016[i] ~ dbin(pi[i],1);  
  }
  
  # PRIORS ON BETAS
  b[1:Nx] ~ dmnorm(mu[1:Nx], tau[1:Nx,1:Nx]) # multivariate normal prior
  
  # Calculate ORs:
  for (l in 1:Nx) {
    OR[l] <- exp(b[l]);
  }
}

N <- length(H3N2_2016) 
Nx <- 2 

mu <- rep(0, Nx)
tau <- diag(0.001, Nx)
data.logistic <- list("N", "Nx", "H3N2_2012", "H3N2_2016", "mu", "tau")
parameters.logistic <-c("b","OR") # Parameters to keep track of
inits.logistic <- function() {list (b= rep(0, Nx))}

set.seed(114011)
logistic.sim_2012_2016_cat1<-jags(data=data.logistic, inits=inits.logistic, parameters.logistic, n.iter=100000,
                                  n.burn=25000, model.file=logistic.model, n.thin=5, n.chains = 3)

print(logistic.sim_2012_2016_cat1, 2)
OR_H3N2_2012_2016_cat1_Bayes <- round(logistic.sim_2012_2016_cat1$BUGSoutput$summary[2,c(5,3,7)], 2)
OR_H3N2_2012_2016_cat1_Bayes

### sensitivity analysis

try7 <- data_orig %>%
  filter(at2012==1 & at2016==1) %>%
  dplyr::select(Obs, H3N2_2012, H3N2_2016, age_2012_cat)

try7$H3N2_2012 <- as.numeric(try7$H3N2_2012)
try7$H3N2_2016 <- as.numeric(try7$H3N2_2016)

# 2x2 table
table(try7[try7$age_2012_cat==0,]$H3N2_2012, try7[try7$age_2012_cat==0,]$H3N2_2016) 
table(try7[try7$age_2012_cat==1,]$H3N2_2012, try7[try7$age_2012_cat==1,]$H3N2_2016) 

# Bayesian analysis of H3N2 2014 on H3N2 2016, under 5s
H3N2_2012 <- try7[try7$age_2012_cat==0,]$H3N2_2012
H3N2_2016 <- try7[try7$age_2012_cat==0,]$H3N2_2016

logistic.model <- function() {
  # SAMPLING DISTRIBUTION
  for (i in 1:N) {
    logit(pi[i]) <- b[1] + b[2]*H3N2_2012[i];
    H3N2_2016[i] ~ dbin(pi[i],1);  
  }
  
  # PRIORS ON BETAS
  b[1:Nx] ~ dmnorm(mu[1:Nx], tau[1:Nx,1:Nx]) # multivariate normal prior
  
  # Calculate ORs:
  for (l in 1:Nx) {
    OR[l] <- exp(b[l]);
  }
}

N <- length(H3N2_2016) 
Nx <- 2 

mu <- rep(0, Nx)
tau <- diag(0.001, Nx)
data.logistic <- list("N", "Nx", "H3N2_2012", "H3N2_2016", "mu", "tau")
parameters.logistic <-c("b","OR") # Parameters to keep track of
inits.logistic <- function() {list (b= rep(0, Nx))}

set.seed(114011)
logistic.sim_2012_2016_cat0<-jags(data=data.logistic, inits=inits.logistic, parameters.logistic, n.iter=100000,
                                  n.burn=25000, model.file=logistic.model, n.thin=5, n.chains = 3)

print(logistic.sim_2012_2016_cat0, 2)
OR_H3N2_2012_2016_cat0_Bayes <- round(logistic.sim_2012_2016_cat0$BUGSoutput$summary[2,c(5,3,7)], 2)
OR_H3N2_2012_2016_cat0_Bayes

# Bayesian analysis of H3N2 2014 on H3N2 2016, over 5s
H3N2_2012 <- try7[try7$age_2012_cat==1,]$H3N2_2012
H3N2_2016 <- try7[try7$age_2012_cat==1,]$H3N2_2016

logistic.model <- function() {
  # SAMPLING DISTRIBUTION
  for (i in 1:N) {
    logit(pi[i]) <- b[1] + b[2]*H3N2_2012[i];
    H3N2_2016[i] ~ dbin(pi[i],1);  
  }
  
  # PRIORS ON BETAS
  b[1:Nx] ~ dmnorm(mu[1:Nx], tau[1:Nx,1:Nx]) # multivariate normal prior
  
  # Calculate ORs:
  for (l in 1:Nx) {
    OR[l] <- exp(b[l]);
  }
}

N <- length(H3N2_2016) 
Nx <- 2 

mu <- rep(0, Nx)
tau <- diag(0.001, Nx)
data.logistic <- list("N", "Nx", "H3N2_2012", "H3N2_2016", "mu", "tau")
parameters.logistic <-c("b","OR") # Parameters to keep track of
inits.logistic <- function() {list (b= rep(0, Nx))}

set.seed(114011)
logistic.sim_2012_2016_cat1<-jags(data=data.logistic, inits=inits.logistic, parameters.logistic, n.iter=100000,
                                  n.burn=25000, model.file=logistic.model, n.thin=5, n.chains = 3)

print(logistic.sim_2012_2016_cat1, 2)
OR_H3N2_2012_2016_cat1_Bayes <- round(logistic.sim_2012_2016_cat1$BUGSoutput$summary[2,c(5,3,7)], 2)
OR_H3N2_2012_2016_cat1_Bayes


# 2012 to 2017

try7 <- data_orig %>%
  filter(at2012==1 & at2017==1 & H3N2_2013==0 & H3N2_2014==0 & H3N2_2015==0 & H3N2_2016==0) %>%
  dplyr::select(Obs, H3N2_2012, H3N2_2017, age_2012_cat)

try7$H3N2_2012 <- as.numeric(try7$H3N2_2012)
try7$H3N2_2017 <- as.numeric(try7$H3N2_2017)

# 2x2 table
table(try7[try7$age_2012_cat==0,]$H3N2_2012, try7[try7$age_2012_cat==0,]$H3N2_2017) 
table(try7[try7$age_2012_cat==1,]$H3N2_2012, try7[try7$age_2012_cat==1,]$H3N2_2017) 

# Bayesian analysis of H3N2 2014 on H3N2 2016, under 5s
H3N2_2012 <- try7[try7$age_2012_cat==0,]$H3N2_2012
H3N2_2017 <- try7[try7$age_2012_cat==0,]$H3N2_2017

logistic.model <- function() {
  # SAMPLING DISTRIBUTION
  for (i in 1:N) {
    logit(pi[i]) <- b[1] + b[2]*H3N2_2012[i];
    H3N2_2017[i] ~ dbin(pi[i],1);  
  }
  
  # PRIORS ON BETAS
  b[1:Nx] ~ dmnorm(mu[1:Nx], tau[1:Nx,1:Nx]) # multivariate normal prior
  
  # Calculate ORs:
  for (l in 1:Nx) {
    OR[l] <- exp(b[l]);
  }
}

N <- length(H3N2_2017) 
Nx <- 2 

mu <- rep(0, Nx)
tau <- diag(0.001, Nx)
data.logistic <- list("N", "Nx", "H3N2_2012", "H3N2_2017", "mu", "tau")
parameters.logistic <-c("b","OR") # Parameters to keep track of
inits.logistic <- function() {list (b= rep(0, Nx))}

set.seed(114011)
logistic.sim_2012_2017_cat0<-jags(data=data.logistic, inits=inits.logistic, parameters.logistic, n.iter=100000,
                                  n.burn=25000, model.file=logistic.model, n.thin=5, n.chains = 3)

print(logistic.sim_2012_2017_cat0, 2)
OR_H3N2_2012_2017_cat0_Bayes <- round(logistic.sim_2012_2017_cat0$BUGSoutput$summary[2,c(5,3,7)], 2)
OR_H3N2_2012_2017_cat0_Bayes

# Bayesian analysis of H3N2 2014 on H3N2 2016, over 5s
H3N2_2012 <- try7[try7$age_2012_cat==1,]$H3N2_2012
H3N2_2017 <- try7[try7$age_2012_cat==1,]$H3N2_2017

logistic.model <- function() {
  # SAMPLING DISTRIBUTION
  for (i in 1:N) {
    logit(pi[i]) <- b[1] + b[2]*H3N2_2012[i];
    H3N2_2017[i] ~ dbin(pi[i],1);  
  }
  
  # PRIORS ON BETAS
  b[1:Nx] ~ dmnorm(mu[1:Nx], tau[1:Nx,1:Nx]) # multivariate normal prior
  
  # Calculate ORs:
  for (l in 1:Nx) {
    OR[l] <- exp(b[l]);
  }
}

N <- length(H3N2_2017) 
Nx <- 2 

mu <- rep(0, Nx)
tau <- diag(0.001, Nx)
data.logistic <- list("N", "Nx", "H3N2_2012", "H3N2_2017", "mu", "tau")
parameters.logistic <-c("b","OR") # Parameters to keep track of
inits.logistic <- function() {list (b= rep(0, Nx))}

set.seed(114011)
logistic.sim_2012_2017_cat1<-jags(data=data.logistic, inits=inits.logistic, parameters.logistic, n.iter=100000,
                                  n.burn=25000, model.file=logistic.model, n.thin=5, n.chains = 3)

print(logistic.sim_2012_2017_cat1, 2)
OR_H3N2_2012_2017_cat1_Bayes <- round(logistic.sim_2012_2017_cat1$BUGSoutput$summary[2,c(5,3,7)], 2)
OR_H3N2_2012_2017_cat1_Bayes

### sensitivity analysis 

try7 <- data_orig %>%
  filter(at2012==1 & at2017==1) %>%
  dplyr::select(Obs, H3N2_2012, H3N2_2017, age_2012_cat)

try7$H3N2_2012 <- as.numeric(try7$H3N2_2012)
try7$H3N2_2017 <- as.numeric(try7$H3N2_2017)

# 2x2 table
table(try7[try7$age_2012_cat==0,]$H3N2_2012, try7[try7$age_2012_cat==0,]$H3N2_2017) 
table(try7[try7$age_2012_cat==1,]$H3N2_2012, try7[try7$age_2012_cat==1,]$H3N2_2017) 

# Bayesian analysis of H3N2 2014 on H3N2 2016, under 5s
H3N2_2012 <- try7[try7$age_2012_cat==0,]$H3N2_2012
H3N2_2017 <- try7[try7$age_2012_cat==0,]$H3N2_2017

logistic.model <- function() {
  # SAMPLING DISTRIBUTION
  for (i in 1:N) {
    logit(pi[i]) <- b[1] + b[2]*H3N2_2012[i];
    H3N2_2017[i] ~ dbin(pi[i],1);  
  }
  
  # PRIORS ON BETAS
  b[1:Nx] ~ dmnorm(mu[1:Nx], tau[1:Nx,1:Nx]) # multivariate normal prior
  
  # Calculate ORs:
  for (l in 1:Nx) {
    OR[l] <- exp(b[l]);
  }
}

N <- length(H3N2_2017) 
Nx <- 2 

mu <- rep(0, Nx)
tau <- diag(0.001, Nx)
data.logistic <- list("N", "Nx", "H3N2_2012", "H3N2_2017", "mu", "tau")
parameters.logistic <-c("b","OR") # Parameters to keep track of
inits.logistic <- function() {list (b= rep(0, Nx))}

set.seed(114011)
logistic.sim_2012_2017_cat0<-jags(data=data.logistic, inits=inits.logistic, parameters.logistic, n.iter=100000,
                                  n.burn=25000, model.file=logistic.model, n.thin=5, n.chains = 3)

print(logistic.sim_2012_2017_cat0, 2)
OR_H3N2_2012_2017_cat0_Bayes <- round(logistic.sim_2012_2017_cat0$BUGSoutput$summary[2,c(5,3,7)], 2)
OR_H3N2_2012_2017_cat0_Bayes

# Bayesian analysis of H3N2 2014 on H3N2 2016, over 5s
H3N2_2012 <- try7[try7$age_2012_cat==1,]$H3N2_2012
H3N2_2017 <- try7[try7$age_2012_cat==1,]$H3N2_2017

logistic.model <- function() {
  # SAMPLING DISTRIBUTION
  for (i in 1:N) {
    logit(pi[i]) <- b[1] + b[2]*H3N2_2012[i];
    H3N2_2017[i] ~ dbin(pi[i],1);  
  }
  
  # PRIORS ON BETAS
  b[1:Nx] ~ dmnorm(mu[1:Nx], tau[1:Nx,1:Nx]) # multivariate normal prior
  
  # Calculate ORs:
  for (l in 1:Nx) {
    OR[l] <- exp(b[l]);
  }
}

N <- length(H3N2_2017) 
Nx <- 2 

mu <- rep(0, Nx)
tau <- diag(0.001, Nx)
data.logistic <- list("N", "Nx", "H3N2_2012", "H3N2_2017", "mu", "tau")
parameters.logistic <-c("b","OR") # Parameters to keep track of
inits.logistic <- function() {list (b= rep(0, Nx))}

set.seed(114011)
logistic.sim_2012_2017_cat1<-jags(data=data.logistic, inits=inits.logistic, parameters.logistic, n.iter=100000,
                                  n.burn=25000, model.file=logistic.model, n.thin=5, n.chains = 3)

print(logistic.sim_2012_2017_cat1, 2)
OR_H3N2_2012_2017_cat1_Bayes <- round(logistic.sim_2012_2017_cat1$BUGSoutput$summary[2,c(5,3,7)], 2)
OR_H3N2_2012_2017_cat1_Bayes

# 2012 to 2019

try7 <- data_orig %>%
  filter(at2012==1 & at2019==1 & H3N2_2013==0 & H3N2_2014==0 & H3N2_2015==0 & H3N2_2016==0 & H3N2_2017==0 & H3N2_2018==0) %>%
  dplyr::select(Obs, H3N2_2012, H3N2_2019, age_2012_cat)

try7$H3N2_2012 <- as.numeric(try7$H3N2_2012)
try7$H3N2_2019 <- as.numeric(try7$H3N2_2019)

# 2x2 table
table(try7[try7$age_2012_cat==0,]$H3N2_2012, try7[try7$age_2012_cat==0,]$H3N2_2019) 
table(try7[try7$age_2012_cat==1,]$H3N2_2012, try7[try7$age_2012_cat==1,]$H3N2_2019) 

# Bayesian analysis of H3N2 2014 on H3N2 2016, under 5s
H3N2_2012 <- try7[try7$age_2012_cat==0,]$H3N2_2012
H3N2_2019 <- try7[try7$age_2012_cat==0,]$H3N2_2019

logistic.model <- function() {
  # SAMPLING DISTRIBUTION
  for (i in 1:N) {
    logit(pi[i]) <- b[1] + b[2]*H3N2_2012[i];
    H3N2_2019[i] ~ dbin(pi[i],1);  
  }
  
  # PRIORS ON BETAS
  b[1:Nx] ~ dmnorm(mu[1:Nx], tau[1:Nx,1:Nx]) # multivariate normal prior
  
  # Calculate ORs:
  for (l in 1:Nx) {
    OR[l] <- exp(b[l]);
  }
}

N <- length(H3N2_2019) 
Nx <- 2 

mu <- rep(0, Nx)
tau <- diag(0.001, Nx)
data.logistic <- list("N", "Nx", "H3N2_2012", "H3N2_2019", "mu", "tau")
parameters.logistic <-c("b","OR") # Parameters to keep track of
inits.logistic <- function() {list (b= rep(0, Nx))}

set.seed(114011)
logistic.sim_2012_2019_cat0<-jags(data=data.logistic, inits=inits.logistic, parameters.logistic, n.iter=100000,
                                  n.burn=25000, model.file=logistic.model, n.thin=5, n.chains = 3)

print(logistic.sim_2012_2019_cat0, 2)
OR_H3N2_2012_2019_cat0_Bayes <- round(logistic.sim_2012_2019_cat0$BUGSoutput$summary[2,c(5,3,7)], 2)
OR_H3N2_2012_2019_cat0_Bayes

# Bayesian analysis of H3N2 2014 on H3N2 2016, over 5s
H3N2_2012 <- try7[try7$age_2012_cat==1,]$H3N2_2012
H3N2_2019 <- try7[try7$age_2012_cat==1,]$H3N2_2019

logistic.model <- function() {
  # SAMPLING DISTRIBUTION
  for (i in 1:N) {
    logit(pi[i]) <- b[1] + b[2]*H3N2_2012[i];
    H3N2_2019[i] ~ dbin(pi[i],1);  
  }
  
  # PRIORS ON BETAS
  b[1:Nx] ~ dmnorm(mu[1:Nx], tau[1:Nx,1:Nx]) # multivariate normal prior
  
  # Calculate ORs:
  for (l in 1:Nx) {
    OR[l] <- exp(b[l]);
  }
}

N <- length(H3N2_2019) 
Nx <- 2 

mu <- rep(0, Nx)
tau <- diag(0.001, Nx)
data.logistic <- list("N", "Nx", "H3N2_2012", "H3N2_2019", "mu", "tau")
parameters.logistic <-c("b","OR") # Parameters to keep track of
inits.logistic <- function() {list (b= rep(0, Nx))}

set.seed(114011)
logistic.sim_2012_2019_cat1<-jags(data=data.logistic, inits=inits.logistic, parameters.logistic, n.iter=100000,
                                  n.burn=25000, model.file=logistic.model, n.thin=5, n.chains = 3)

print(logistic.sim_2012_2019_cat1, 2)
OR_H3N2_2012_2019_cat1_Bayes <- round(logistic.sim_2012_2019_cat1$BUGSoutput$summary[2,c(5,3,7)], 2)
OR_H3N2_2012_2019_cat1_Bayes

### sensitivity analysis 


try7 <- data_orig %>%
  filter(at2012==1 & at2019==1 ) %>%
  dplyr::select(Obs, H3N2_2012, H3N2_2019, age_2012_cat)

try7$H3N2_2012 <- as.numeric(try7$H3N2_2012)
try7$H3N2_2019 <- as.numeric(try7$H3N2_2019)

# 2x2 table
table(try7[try7$age_2012_cat==0,]$H3N2_2012, try7[try7$age_2012_cat==0,]$H3N2_2019) 
table(try7[try7$age_2012_cat==1,]$H3N2_2012, try7[try7$age_2012_cat==1,]$H3N2_2019) 

# Bayesian analysis of H3N2 2014 on H3N2 2016, under 5s
H3N2_2012 <- try7[try7$age_2012_cat==0,]$H3N2_2012
H3N2_2019 <- try7[try7$age_2012_cat==0,]$H3N2_2019

logistic.model <- function() {
  # SAMPLING DISTRIBUTION
  for (i in 1:N) {
    logit(pi[i]) <- b[1] + b[2]*H3N2_2012[i];
    H3N2_2019[i] ~ dbin(pi[i],1);  
  }
  
  # PRIORS ON BETAS
  b[1:Nx] ~ dmnorm(mu[1:Nx], tau[1:Nx,1:Nx]) # multivariate normal prior
  
  # Calculate ORs:
  for (l in 1:Nx) {
    OR[l] <- exp(b[l]);
  }
}

N <- length(H3N2_2019) 
Nx <- 2 

mu <- rep(0, Nx)
tau <- diag(0.001, Nx)
data.logistic <- list("N", "Nx", "H3N2_2012", "H3N2_2019", "mu", "tau")
parameters.logistic <-c("b","OR") # Parameters to keep track of
inits.logistic <- function() {list (b= rep(0, Nx))}

set.seed(114011)
logistic.sim_2012_2019_cat0<-jags(data=data.logistic, inits=inits.logistic, parameters.logistic, n.iter=100000,
                                  n.burn=25000, model.file=logistic.model, n.thin=5, n.chains = 3)

print(logistic.sim_2012_2019_cat0, 2)
OR_H3N2_2012_2019_cat0_Bayes <- round(logistic.sim_2012_2019_cat0$BUGSoutput$summary[2,c(5,3,7)], 2)
OR_H3N2_2012_2019_cat0_Bayes

# Bayesian analysis of H3N2 2014 on H3N2 2016, over 5s
H3N2_2012 <- try7[try7$age_2012_cat==1,]$H3N2_2012
H3N2_2019 <- try7[try7$age_2012_cat==1,]$H3N2_2019

logistic.model <- function() {
  # SAMPLING DISTRIBUTION
  for (i in 1:N) {
    logit(pi[i]) <- b[1] + b[2]*H3N2_2012[i];
    H3N2_2019[i] ~ dbin(pi[i],1);  
  }
  
  # PRIORS ON BETAS
  b[1:Nx] ~ dmnorm(mu[1:Nx], tau[1:Nx,1:Nx]) # multivariate normal prior
  
  # Calculate ORs:
  for (l in 1:Nx) {
    OR[l] <- exp(b[l]);
  }
}

N <- length(H3N2_2019) 
Nx <- 2 

mu <- rep(0, Nx)
tau <- diag(0.001, Nx)
data.logistic <- list("N", "Nx", "H3N2_2012", "H3N2_2019", "mu", "tau")
parameters.logistic <-c("b","OR") # Parameters to keep track of
inits.logistic <- function() {list (b= rep(0, Nx))}

set.seed(114011)
logistic.sim_2012_2019_cat1<-jags(data=data.logistic, inits=inits.logistic, parameters.logistic, n.iter=100000,
                                  n.burn=25000, model.file=logistic.model, n.thin=5, n.chains = 3)

print(logistic.sim_2012_2019_cat1, 2)
OR_H3N2_2012_2019_cat1_Bayes <- round(logistic.sim_2012_2019_cat1$BUGSoutput$summary[2,c(5,3,7)], 2)
OR_H3N2_2012_2019_cat1_Bayes

#2013 to 2016

try7 <- data_orig %>%
  filter(at2013==1 & at2016==1 & H3N2_2014==0 & H3N2_2015==0 ) %>%
  dplyr::select(Obs, H3N2_2013, H3N2_2016, age_2013_cat)

try7$H3N2_2013 <- as.numeric(try7$H3N2_2013)
try7$H3N2_2016 <- as.numeric(try7$H3N2_2016)

# 2x2 table
table(try7[try7$age_2013_cat==0,]$H3N2_2013, try7[try7$age_2013_cat==0,]$H3N2_2016) 
table(try7[try7$age_2013_cat==1,]$H3N2_2013, try7[try7$age_2013_cat==1,]$H3N2_2016) 

# Bayesian analysis of H3N2 2014 on H3N2 2016, under 5s
H3N2_2013 <- try7[try7$age_2013_cat==0,]$H3N2_2013
H3N2_2016 <- try7[try7$age_2013_cat==0,]$H3N2_2016

logistic.model <- function() {
  # SAMPLING DISTRIBUTION
  for (i in 1:N) {
    logit(pi[i]) <- b[1] + b[2]*H3N2_2013[i];
    H3N2_2016[i] ~ dbin(pi[i],1);  
  }
  
  # PRIORS ON BETAS
  b[1:Nx] ~ dmnorm(mu[1:Nx], tau[1:Nx,1:Nx]) # multivariate normal prior
  
  # Calculate ORs:
  for (l in 1:Nx) {
    OR[l] <- exp(b[l]);
  }
}

N <- length(H3N2_2016) 
Nx <- 2 

mu <- rep(0, Nx)
tau <- diag(0.001, Nx)
data.logistic <- list("N", "Nx", "H3N2_2013", "H3N2_2016", "mu", "tau")
parameters.logistic <-c("b","OR") # Parameters to keep track of
inits.logistic <- function() {list (b= rep(0, Nx))}

set.seed(114011)
logistic.sim_2013_2016_cat0<-jags(data=data.logistic, inits=inits.logistic, parameters.logistic, n.iter=100000,
                                  n.burn=25000, model.file=logistic.model, n.thin=5, n.chains = 3)

print(logistic.sim_2013_2016_cat0, 2)
OR_H3N2_2013_2016_cat0_Bayes <- round(logistic.sim_2013_2016_cat0$BUGSoutput$summary[2,c(5,3,7)], 2)
OR_H3N2_2013_2016_cat0_Bayes

# Bayesian analysis of H3N2 2014 on H3N2 2016, over 5s
H3N2_2013 <- try7[try7$age_2013_cat==1,]$H3N2_2013
H3N2_2016 <- try7[try7$age_2013_cat==1,]$H3N2_2016

logistic.model <- function() {
  # SAMPLING DISTRIBUTION
  for (i in 1:N) {
    logit(pi[i]) <- b[1] + b[2]*H3N2_2013[i];
    H3N2_2016[i] ~ dbin(pi[i],1);  
  }
  
  # PRIORS ON BETAS
  b[1:Nx] ~ dmnorm(mu[1:Nx], tau[1:Nx,1:Nx]) # multivariate normal prior
  
  # Calculate ORs:
  for (l in 1:Nx) {
    OR[l] <- exp(b[l]);
  }
}

N <- length(H3N2_2016) 
Nx <- 2 

mu <- rep(0, Nx)
tau <- diag(0.001, Nx)
data.logistic <- list("N", "Nx", "H3N2_2013", "H3N2_2016", "mu", "tau")
parameters.logistic <-c("b","OR") # Parameters to keep track of
inits.logistic <- function() {list (b= rep(0, Nx))}

set.seed(114011)
logistic.sim_2013_2016_cat1<-jags(data=data.logistic, inits=inits.logistic, parameters.logistic, n.iter=100000,
                                  n.burn=25000, model.file=logistic.model, n.thin=5, n.chains = 3)

print(logistic.sim_2013_2016_cat1, 2)
OR_H3N2_2013_2016_cat1_Bayes <- round(logistic.sim_2013_2016_cat1$BUGSoutput$summary[2,c(5,3,7)], 2)
OR_H3N2_2013_2016_cat1_Bayes

### sensitivity analysis


try7 <- data_orig %>%
  filter(at2013==1 & at2016==1) %>%
  dplyr::select(Obs, H3N2_2013, H3N2_2016, age_2013_cat)

try7$H3N2_2013 <- as.numeric(try7$H3N2_2013)
try7$H3N2_2016 <- as.numeric(try7$H3N2_2016)

# 2x2 table
table(try7[try7$age_2013_cat==0,]$H3N2_2013, try7[try7$age_2013_cat==0,]$H3N2_2016) 
table(try7[try7$age_2013_cat==1,]$H3N2_2013, try7[try7$age_2013_cat==1,]$H3N2_2016) 

# Bayesian analysis of H3N2 2014 on H3N2 2016, under 5s
H3N2_2013 <- try7[try7$age_2013_cat==0,]$H3N2_2013
H3N2_2016 <- try7[try7$age_2013_cat==0,]$H3N2_2016

logistic.model <- function() {
  # SAMPLING DISTRIBUTION
  for (i in 1:N) {
    logit(pi[i]) <- b[1] + b[2]*H3N2_2013[i];
    H3N2_2016[i] ~ dbin(pi[i],1);  
  }
  
  # PRIORS ON BETAS
  b[1:Nx] ~ dmnorm(mu[1:Nx], tau[1:Nx,1:Nx]) # multivariate normal prior
  
  # Calculate ORs:
  for (l in 1:Nx) {
    OR[l] <- exp(b[l]);
  }
}

N <- length(H3N2_2016) 
Nx <- 2 

mu <- rep(0, Nx)
tau <- diag(0.001, Nx)
data.logistic <- list("N", "Nx", "H3N2_2013", "H3N2_2016", "mu", "tau")
parameters.logistic <-c("b","OR") # Parameters to keep track of
inits.logistic <- function() {list (b= rep(0, Nx))}

set.seed(114011)
logistic.sim_2013_2016_cat0<-jags(data=data.logistic, inits=inits.logistic, parameters.logistic, n.iter=100000,
                                  n.burn=25000, model.file=logistic.model, n.thin=5, n.chains = 3)

print(logistic.sim_2013_2016_cat0, 2)
OR_H3N2_2013_2016_cat0_Bayes <- round(logistic.sim_2013_2016_cat0$BUGSoutput$summary[2,c(5,3,7)], 2)
OR_H3N2_2013_2016_cat0_Bayes

# Bayesian analysis of H3N2 2014 on H3N2 2016, over 5s
H3N2_2013 <- try7[try7$age_2013_cat==1,]$H3N2_2013
H3N2_2016 <- try7[try7$age_2013_cat==1,]$H3N2_2016

logistic.model <- function() {
  # SAMPLING DISTRIBUTION
  for (i in 1:N) {
    logit(pi[i]) <- b[1] + b[2]*H3N2_2013[i];
    H3N2_2016[i] ~ dbin(pi[i],1);  
  }
  
  # PRIORS ON BETAS
  b[1:Nx] ~ dmnorm(mu[1:Nx], tau[1:Nx,1:Nx]) # multivariate normal prior
  
  # Calculate ORs:
  for (l in 1:Nx) {
    OR[l] <- exp(b[l]);
  }
}

N <- length(H3N2_2016) 
Nx <- 2 

mu <- rep(0, Nx)
tau <- diag(0.001, Nx)
data.logistic <- list("N", "Nx", "H3N2_2013", "H3N2_2016", "mu", "tau")
parameters.logistic <-c("b","OR") # Parameters to keep track of
inits.logistic <- function() {list (b= rep(0, Nx))}

set.seed(114011)
logistic.sim_2013_2016_cat1<-jags(data=data.logistic, inits=inits.logistic, parameters.logistic, n.iter=100000,
                                  n.burn=25000, model.file=logistic.model, n.thin=5, n.chains = 3)

print(logistic.sim_2013_2016_cat1, 2)
OR_H3N2_2013_2016_cat1_Bayes <- round(logistic.sim_2013_2016_cat1$BUGSoutput$summary[2,c(5,3,7)], 2)
OR_H3N2_2013_2016_cat1_Bayes

#2013 to 2017

try7 <- data_orig %>%
  filter(at2013==1 & at2017==1 & H3N2_2014==0 & H3N2_2015==0 & H3N2_2016==0) %>%
  dplyr::select(Obs, H3N2_2013, H3N2_2017, age_2013_cat)

try7$H3N2_2013 <- as.numeric(try7$H3N2_2013)
try7$H3N2_2017 <- as.numeric(try7$H3N2_2017)

# 2x2 table
table(try7[try7$age_2013_cat==0,]$H3N2_2013, try7[try7$age_2013_cat==0,]$H3N2_2017) 
table(try7[try7$age_2013_cat==1,]$H3N2_2013, try7[try7$age_2013_cat==1,]$H3N2_2017) 

# Bayesian analysis of H3N2 2014 on H3N2 2016, under 5s
H3N2_2013 <- try7[try7$age_2013_cat==0,]$H3N2_2013
H3N2_2017 <- try7[try7$age_2013_cat==0,]$H3N2_2017

logistic.model <- function() {
  # SAMPLING DISTRIBUTION
  for (i in 1:N) {
    logit(pi[i]) <- b[1] + b[2]*H3N2_2013[i];
    H3N2_2017[i] ~ dbin(pi[i],1);  
  }
  
  # PRIORS ON BETAS
  b[1:Nx] ~ dmnorm(mu[1:Nx], tau[1:Nx,1:Nx]) # multivariate normal prior
  
  # Calculate ORs:
  for (l in 1:Nx) {
    OR[l] <- exp(b[l]);
  }
}

N <- length(H3N2_2017) 
Nx <- 2 

mu <- rep(0, Nx)
tau <- diag(0.001, Nx)
data.logistic <- list("N", "Nx", "H3N2_2013", "H3N2_2017", "mu", "tau")
parameters.logistic <-c("b","OR") # Parameters to keep track of
inits.logistic <- function() {list (b= rep(0, Nx))}

set.seed(114011)
logistic.sim_2013_2017_cat0<-jags(data=data.logistic, inits=inits.logistic, parameters.logistic, n.iter=100000,
                                  n.burn=25000, model.file=logistic.model, n.thin=5, n.chains = 3)

print(logistic.sim_2013_2017_cat0, 2)
OR_H3N2_2013_2017_cat0_Bayes <- round(logistic.sim_2013_2017_cat0$BUGSoutput$summary[2,c(5,3,7)], 2)
OR_H3N2_2013_2017_cat0_Bayes

# Bayesian analysis of H3N2 2014 on H3N2 2016, over 5s
H3N2_2013 <- try7[try7$age_2013_cat==1,]$H3N2_2013
H3N2_2017 <- try7[try7$age_2013_cat==1,]$H3N2_2017

logistic.model <- function() {
  # SAMPLING DISTRIBUTION
  for (i in 1:N) {
    logit(pi[i]) <- b[1] + b[2]*H3N2_2013[i];
    H3N2_2017[i] ~ dbin(pi[i],1);  
  }
  
  # PRIORS ON BETAS
  b[1:Nx] ~ dmnorm(mu[1:Nx], tau[1:Nx,1:Nx]) # multivariate normal prior
  
  # Calculate ORs:
  for (l in 1:Nx) {
    OR[l] <- exp(b[l]);
  }
}

N <- length(H3N2_2017) 
Nx <- 2 

mu <- rep(0, Nx)
tau <- diag(0.001, Nx)
data.logistic <- list("N", "Nx", "H3N2_2013", "H3N2_2017", "mu", "tau")
parameters.logistic <-c("b","OR") # Parameters to keep track of
inits.logistic <- function() {list (b= rep(0, Nx))}

set.seed(114011)
logistic.sim_2013_2017_cat1<-jags(data=data.logistic, inits=inits.logistic, parameters.logistic, n.iter=100000,
                                  n.burn=25000, model.file=logistic.model, n.thin=5, n.chains = 3)

print(logistic.sim_2013_2017_cat1, 2)
OR_H3N2_2013_2017_cat1_Bayes <- round(logistic.sim_2013_2017_cat1$BUGSoutput$summary[2,c(5,3,7)], 2)
OR_H3N2_2013_2017_cat1_Bayes

### sensitivity analysis


try7 <- data_orig %>%
  filter(at2013==1 & at2017==1 ) %>%
  dplyr::select(Obs, H3N2_2013, H3N2_2017, age_2013_cat)

try7$H3N2_2013 <- as.numeric(try7$H3N2_2013)
try7$H3N2_2017 <- as.numeric(try7$H3N2_2017)

# 2x2 table
table(try7[try7$age_2013_cat==0,]$H3N2_2013, try7[try7$age_2013_cat==0,]$H3N2_2017) 
table(try7[try7$age_2013_cat==1,]$H3N2_2013, try7[try7$age_2013_cat==1,]$H3N2_2017) 

# Bayesian analysis of H3N2 2014 on H3N2 2016, under 5s
H3N2_2013 <- try7[try7$age_2013_cat==0,]$H3N2_2013
H3N2_2017 <- try7[try7$age_2013_cat==0,]$H3N2_2017

logistic.model <- function() {
  # SAMPLING DISTRIBUTION
  for (i in 1:N) {
    logit(pi[i]) <- b[1] + b[2]*H3N2_2013[i];
    H3N2_2017[i] ~ dbin(pi[i],1);  
  }
  
  # PRIORS ON BETAS
  b[1:Nx] ~ dmnorm(mu[1:Nx], tau[1:Nx,1:Nx]) # multivariate normal prior
  
  # Calculate ORs:
  for (l in 1:Nx) {
    OR[l] <- exp(b[l]);
  }
}

N <- length(H3N2_2017) 
Nx <- 2 

mu <- rep(0, Nx)
tau <- diag(0.001, Nx)
data.logistic <- list("N", "Nx", "H3N2_2013", "H3N2_2017", "mu", "tau")
parameters.logistic <-c("b","OR") # Parameters to keep track of
inits.logistic <- function() {list (b= rep(0, Nx))}

set.seed(114011)
logistic.sim_2013_2017_cat0<-jags(data=data.logistic, inits=inits.logistic, parameters.logistic, n.iter=100000,
                                  n.burn=25000, model.file=logistic.model, n.thin=5, n.chains = 3)

print(logistic.sim_2013_2017_cat0, 2)
OR_H3N2_2013_2017_cat0_Bayes <- round(logistic.sim_2013_2017_cat0$BUGSoutput$summary[2,c(5,3,7)], 2)
OR_H3N2_2013_2017_cat0_Bayes

# Bayesian analysis of H3N2 2014 on H3N2 2016, over 5s
H3N2_2013 <- try7[try7$age_2013_cat==1,]$H3N2_2013
H3N2_2017 <- try7[try7$age_2013_cat==1,]$H3N2_2017

logistic.model <- function() {
  # SAMPLING DISTRIBUTION
  for (i in 1:N) {
    logit(pi[i]) <- b[1] + b[2]*H3N2_2013[i];
    H3N2_2017[i] ~ dbin(pi[i],1);  
  }
  
  # PRIORS ON BETAS
  b[1:Nx] ~ dmnorm(mu[1:Nx], tau[1:Nx,1:Nx]) # multivariate normal prior
  
  # Calculate ORs:
  for (l in 1:Nx) {
    OR[l] <- exp(b[l]);
  }
}

N <- length(H3N2_2017) 
Nx <- 2 

mu <- rep(0, Nx)
tau <- diag(0.001, Nx)
data.logistic <- list("N", "Nx", "H3N2_2013", "H3N2_2017", "mu", "tau")
parameters.logistic <-c("b","OR") # Parameters to keep track of
inits.logistic <- function() {list (b= rep(0, Nx))}

set.seed(114011)
logistic.sim_2013_2017_cat1<-jags(data=data.logistic, inits=inits.logistic, parameters.logistic, n.iter=100000,
                                  n.burn=25000, model.file=logistic.model, n.thin=5, n.chains = 3)

print(logistic.sim_2013_2017_cat1, 2)
OR_H3N2_2013_2017_cat1_Bayes <- round(logistic.sim_2013_2017_cat1$BUGSoutput$summary[2,c(5,3,7)], 2)
OR_H3N2_2013_2017_cat1_Bayes

#2013 to 2019

try7 <- data_orig %>%
  filter(at2013==1 & at2019==1 & H3N2_2014==0 & H3N2_2015==0 & H3N2_2016==0 & H3N2_2017==0 & H3N2_2018==0) %>%
  dplyr::select(Obs, H3N2_2013, H3N2_2019, age_2013_cat)

try7$H3N2_2013 <- as.numeric(try7$H3N2_2013)
try7$H3N2_2019 <- as.numeric(try7$H3N2_2019)

# 2x2 table
table(try7[try7$age_2013_cat==0,]$H3N2_2013, try7[try7$age_2013_cat==0,]$H3N2_2019) 
table(try7[try7$age_2013_cat==1,]$H3N2_2013, try7[try7$age_2013_cat==1,]$H3N2_2019) 

# Bayesian analysis of H3N2 2014 on H3N2 2016, under 5s
H3N2_2013 <- try7[try7$age_2013_cat==0,]$H3N2_2013
H3N2_2019 <- try7[try7$age_2013_cat==0,]$H3N2_2019

logistic.model <- function() {
  # SAMPLING DISTRIBUTION
  for (i in 1:N) {
    logit(pi[i]) <- b[1] + b[2]*H3N2_2013[i];
    H3N2_2019[i] ~ dbin(pi[i],1);  
  }
  
  # PRIORS ON BETAS
  b[1:Nx] ~ dmnorm(mu[1:Nx], tau[1:Nx,1:Nx]) # multivariate normal prior
  
  # Calculate ORs:
  for (l in 1:Nx) {
    OR[l] <- exp(b[l]);
  }
}

N <- length(H3N2_2019) 
Nx <- 2 

mu <- rep(0, Nx)
tau <- diag(0.001, Nx)
data.logistic <- list("N", "Nx", "H3N2_2013", "H3N2_2019", "mu", "tau")
parameters.logistic <-c("b","OR") # Parameters to keep track of
inits.logistic <- function() {list (b= rep(0, Nx))}

set.seed(114011)
logistic.sim_2013_2019_cat0<-jags(data=data.logistic, inits=inits.logistic, parameters.logistic, n.iter=100000,
                                  n.burn=25000, model.file=logistic.model, n.thin=5, n.chains = 3)

print(logistic.sim_2013_2019_cat0, 2)
OR_H3N2_2013_2019_cat0_Bayes <- round(logistic.sim_2013_2019_cat0$BUGSoutput$summary[2,c(5,3,7)], 2)
OR_H3N2_2013_2019_cat0_Bayes

# Bayesian analysis of H3N2 2014 on H3N2 2016, over 5s
H3N2_2013 <- try7[try7$age_2013_cat==1,]$H3N2_2013
H3N2_2019 <- try7[try7$age_2013_cat==1,]$H3N2_2019

logistic.model <- function() {
  # SAMPLING DISTRIBUTION
  for (i in 1:N) {
    logit(pi[i]) <- b[1] + b[2]*H3N2_2013[i];
    H3N2_2019[i] ~ dbin(pi[i],1);  
  }
  
  # PRIORS ON BETAS
  b[1:Nx] ~ dmnorm(mu[1:Nx], tau[1:Nx,1:Nx]) # multivariate normal prior
  
  # Calculate ORs:
  for (l in 1:Nx) {
    OR[l] <- exp(b[l]);
  }
}

N <- length(H3N2_2019) 
Nx <- 2 

mu <- rep(0, Nx)
tau <- diag(0.001, Nx)
data.logistic <- list("N", "Nx", "H3N2_2013", "H3N2_2019", "mu", "tau")
parameters.logistic <-c("b","OR") # Parameters to keep track of
inits.logistic <- function() {list (b= rep(0, Nx))}

set.seed(114011)
logistic.sim_2013_2019_cat1<-jags(data=data.logistic, inits=inits.logistic, parameters.logistic, n.iter=100000,
                                  n.burn=25000, model.file=logistic.model, n.thin=5, n.chains = 3)

print(logistic.sim_2013_2019_cat1, 2)
OR_H3N2_2013_2019_cat1_Bayes <- round(logistic.sim_2013_2019_cat1$BUGSoutput$summary[2,c(5,3,7)], 2)
OR_H3N2_2013_2019_cat1_Bayes

### sensitivity analysis 

try7 <- data_orig %>%
  filter(at2013==1 & at2019==1 ) %>%
  dplyr::select(Obs, H3N2_2013, H3N2_2019, age_2013_cat)

try7$H3N2_2013 <- as.numeric(try7$H3N2_2013)
try7$H3N2_2019 <- as.numeric(try7$H3N2_2019)

# 2x2 table
table(try7[try7$age_2013_cat==0,]$H3N2_2013, try7[try7$age_2013_cat==0,]$H3N2_2019) 
table(try7[try7$age_2013_cat==1,]$H3N2_2013, try7[try7$age_2013_cat==1,]$H3N2_2019) 

# Bayesian analysis of H3N2 2014 on H3N2 2016, under 5s
H3N2_2013 <- try7[try7$age_2013_cat==0,]$H3N2_2013
H3N2_2019 <- try7[try7$age_2013_cat==0,]$H3N2_2019

logistic.model <- function() {
  # SAMPLING DISTRIBUTION
  for (i in 1:N) {
    logit(pi[i]) <- b[1] + b[2]*H3N2_2013[i];
    H3N2_2019[i] ~ dbin(pi[i],1);  
  }
  
  # PRIORS ON BETAS
  b[1:Nx] ~ dmnorm(mu[1:Nx], tau[1:Nx,1:Nx]) # multivariate normal prior
  
  # Calculate ORs:
  for (l in 1:Nx) {
    OR[l] <- exp(b[l]);
  }
}

N <- length(H3N2_2019) 
Nx <- 2 

mu <- rep(0, Nx)
tau <- diag(0.001, Nx)
data.logistic <- list("N", "Nx", "H3N2_2013", "H3N2_2019", "mu", "tau")
parameters.logistic <-c("b","OR") # Parameters to keep track of
inits.logistic <- function() {list (b= rep(0, Nx))}

set.seed(114011)
logistic.sim_2013_2019_cat0<-jags(data=data.logistic, inits=inits.logistic, parameters.logistic, n.iter=100000,
                                  n.burn=25000, model.file=logistic.model, n.thin=5, n.chains = 3)

print(logistic.sim_2013_2019_cat0, 2)
OR_H3N2_2013_2019_cat0_Bayes <- round(logistic.sim_2013_2019_cat0$BUGSoutput$summary[2,c(5,3,7)], 2)
OR_H3N2_2013_2019_cat0_Bayes

# Bayesian analysis of H3N2 2014 on H3N2 2016, over 5s
H3N2_2013 <- try7[try7$age_2013_cat==1,]$H3N2_2013
H3N2_2019 <- try7[try7$age_2013_cat==1,]$H3N2_2019

logistic.model <- function() {
  # SAMPLING DISTRIBUTION
  for (i in 1:N) {
    logit(pi[i]) <- b[1] + b[2]*H3N2_2013[i];
    H3N2_2019[i] ~ dbin(pi[i],1);  
  }
  
  # PRIORS ON BETAS
  b[1:Nx] ~ dmnorm(mu[1:Nx], tau[1:Nx,1:Nx]) # multivariate normal prior
  
  # Calculate ORs:
  for (l in 1:Nx) {
    OR[l] <- exp(b[l]);
  }
}

N <- length(H3N2_2019) 
Nx <- 2 

mu <- rep(0, Nx)
tau <- diag(0.001, Nx)
data.logistic <- list("N", "Nx", "H3N2_2013", "H3N2_2019", "mu", "tau")
parameters.logistic <-c("b","OR") # Parameters to keep track of
inits.logistic <- function() {list (b= rep(0, Nx))}

set.seed(114011)
logistic.sim_2013_2019_cat1<-jags(data=data.logistic, inits=inits.logistic, parameters.logistic, n.iter=100000,
                                  n.burn=25000, model.file=logistic.model, n.thin=5, n.chains = 3)

print(logistic.sim_2013_2019_cat1, 2)
OR_H3N2_2013_2019_cat1_Bayes <- round(logistic.sim_2013_2019_cat1$BUGSoutput$summary[2,c(5,3,7)], 2)
OR_H3N2_2013_2019_cat1_Bayes

#2014 to 2017

try7 <- data_orig %>%
  filter(at2014==1 & at2017==1 & H3N2_2015==0 & H3N2_2016==0 ) %>%
  dplyr::select(Obs, H3N2_2014, H3N2_2017, age_2014_cat)

try7$H3N2_2014 <- as.numeric(try7$H3N2_2014)
try7$H3N2_2017 <- as.numeric(try7$H3N2_2017)

# 2x2 table
table(try7[try7$age_2014_cat==0,]$H3N2_2014, try7[try7$age_2014_cat==0,]$H3N2_2017) 
table(try7[try7$age_2014_cat==1,]$H3N2_2014, try7[try7$age_2014_cat==1,]$H3N2_2017) 

# Bayesian analysis of H3N2 2014 on H3N2 2016, under 5s
H3N2_2014 <- try7[try7$age_2014_cat==0,]$H3N2_2014
H3N2_2017 <- try7[try7$age_2014_cat==0,]$H3N2_2017

logistic.model <- function() {
  # SAMPLING DISTRIBUTION
  for (i in 1:N) {
    logit(pi[i]) <- b[1] + b[2]*H3N2_2014[i];
    H3N2_2017[i] ~ dbin(pi[i],1);  
  }
  
  # PRIORS ON BETAS
  b[1:Nx] ~ dmnorm(mu[1:Nx], tau[1:Nx,1:Nx]) # multivariate normal prior
  
  # Calculate ORs:
  for (l in 1:Nx) {
    OR[l] <- exp(b[l]);
  }
}

N <- length(H3N2_2017) 
Nx <- 2 

mu <- rep(0, Nx)
tau <- diag(0.001, Nx)
data.logistic <- list("N", "Nx", "H3N2_2014", "H3N2_2017", "mu", "tau")
parameters.logistic <-c("b","OR") # Parameters to keep track of
inits.logistic <- function() {list (b= rep(0, Nx))}

set.seed(114011)
logistic.sim_2014_2017_cat0<-jags(data=data.logistic, inits=inits.logistic, parameters.logistic, n.iter=100000,
                                  n.burn=25000, model.file=logistic.model, n.thin=5, n.chains = 3)

print(logistic.sim_2014_2017_cat0, 2)
OR_H3N2_2014_2017_cat0_Bayes <- round(logistic.sim_2014_2017_cat0$BUGSoutput$summary[2,c(5,3,7)], 2)
OR_H3N2_2014_2017_cat0_Bayes

# Bayesian analysis of H3N2 2014 on H3N2 2016, over 5s
H3N2_2014 <- try7[try7$age_2014_cat==1,]$H3N2_2014
H3N2_2017 <- try7[try7$age_2014_cat==1,]$H3N2_2017

logistic.model <- function() {
  # SAMPLING DISTRIBUTION
  for (i in 1:N) {
    logit(pi[i]) <- b[1] + b[2]*H3N2_2014[i];
    H3N2_2017[i] ~ dbin(pi[i],1);  
  }
  
  # PRIORS ON BETAS
  b[1:Nx] ~ dmnorm(mu[1:Nx], tau[1:Nx,1:Nx]) # multivariate normal prior
  
  # Calculate ORs:
  for (l in 1:Nx) {
    OR[l] <- exp(b[l]);
  }
}

N <- length(H3N2_2017) 
Nx <- 2 

mu <- rep(0, Nx)
tau <- diag(0.001, Nx)
data.logistic <- list("N", "Nx", "H3N2_2014", "H3N2_2017", "mu", "tau")
parameters.logistic <-c("b","OR") # Parameters to keep track of
inits.logistic <- function() {list (b= rep(0, Nx))}

set.seed(114011)
logistic.sim_2014_2017_cat1<-jags(data=data.logistic, inits=inits.logistic, parameters.logistic, n.iter=100000,
                                  n.burn=25000, model.file=logistic.model, n.thin=5, n.chains = 3)

print(logistic.sim_2014_2017_cat1, 2)
OR_H3N2_2014_2017_cat1_Bayes <- round(logistic.sim_2014_2017_cat1$BUGSoutput$summary[2,c(5,3,7)], 2)
OR_H3N2_2014_2017_cat1_Bayes

### sensitivity analysis 


try7 <- data_orig %>%
  filter(at2014==1 & at2017==1 ) %>%
  dplyr::select(Obs, H3N2_2014, H3N2_2017, age_2014_cat)

try7$H3N2_2014 <- as.numeric(try7$H3N2_2014)
try7$H3N2_2017 <- as.numeric(try7$H3N2_2017)

# 2x2 table
table(try7[try7$age_2014_cat==0,]$H3N2_2014, try7[try7$age_2014_cat==0,]$H3N2_2017) 
table(try7[try7$age_2014_cat==1,]$H3N2_2014, try7[try7$age_2014_cat==1,]$H3N2_2017) 

# Bayesian analysis of H3N2 2014 on H3N2 2016, under 5s
H3N2_2014 <- try7[try7$age_2014_cat==0,]$H3N2_2014
H3N2_2017 <- try7[try7$age_2014_cat==0,]$H3N2_2017

logistic.model <- function() {
  # SAMPLING DISTRIBUTION
  for (i in 1:N) {
    logit(pi[i]) <- b[1] + b[2]*H3N2_2014[i];
    H3N2_2017[i] ~ dbin(pi[i],1);  
  }
  
  # PRIORS ON BETAS
  b[1:Nx] ~ dmnorm(mu[1:Nx], tau[1:Nx,1:Nx]) # multivariate normal prior
  
  # Calculate ORs:
  for (l in 1:Nx) {
    OR[l] <- exp(b[l]);
  }
}

N <- length(H3N2_2017) 
Nx <- 2 

mu <- rep(0, Nx)
tau <- diag(0.001, Nx)
data.logistic <- list("N", "Nx", "H3N2_2014", "H3N2_2017", "mu", "tau")
parameters.logistic <-c("b","OR") # Parameters to keep track of
inits.logistic <- function() {list (b= rep(0, Nx))}

set.seed(114011)
logistic.sim_2014_2017_cat0<-jags(data=data.logistic, inits=inits.logistic, parameters.logistic, n.iter=100000,
                                  n.burn=25000, model.file=logistic.model, n.thin=5, n.chains = 3)

print(logistic.sim_2014_2017_cat0, 2)
OR_H3N2_2014_2017_cat0_Bayes <- round(logistic.sim_2014_2017_cat0$BUGSoutput$summary[2,c(5,3,7)], 2)
OR_H3N2_2014_2017_cat0_Bayes

# Bayesian analysis of H3N2 2014 on H3N2 2016, over 5s
H3N2_2014 <- try7[try7$age_2014_cat==1,]$H3N2_2014
H3N2_2017 <- try7[try7$age_2014_cat==1,]$H3N2_2017

logistic.model <- function() {
  # SAMPLING DISTRIBUTION
  for (i in 1:N) {
    logit(pi[i]) <- b[1] + b[2]*H3N2_2014[i];
    H3N2_2017[i] ~ dbin(pi[i],1);  
  }
  
  # PRIORS ON BETAS
  b[1:Nx] ~ dmnorm(mu[1:Nx], tau[1:Nx,1:Nx]) # multivariate normal prior
  
  # Calculate ORs:
  for (l in 1:Nx) {
    OR[l] <- exp(b[l]);
  }
}

N <- length(H3N2_2017) 
Nx <- 2 

mu <- rep(0, Nx)
tau <- diag(0.001, Nx)
data.logistic <- list("N", "Nx", "H3N2_2014", "H3N2_2017", "mu", "tau")
parameters.logistic <-c("b","OR") # Parameters to keep track of
inits.logistic <- function() {list (b= rep(0, Nx))}

set.seed(114011)
logistic.sim_2014_2017_cat1<-jags(data=data.logistic, inits=inits.logistic, parameters.logistic, n.iter=100000,
                                  n.burn=25000, model.file=logistic.model, n.thin=5, n.chains = 3)

print(logistic.sim_2014_2017_cat1, 2)
OR_H3N2_2014_2017_cat1_Bayes <- round(logistic.sim_2014_2017_cat1$BUGSoutput$summary[2,c(5,3,7)], 2)
OR_H3N2_2014_2017_cat1_Bayes

#2014 to 2019

try7 <- data_orig %>%
  filter(at2014==1 & at2019==1 & H3N2_2015==0 & H3N2_2016==0 & H3N2_2017==0 & H3N2_2018==0) %>%
  dplyr::select(Obs, H3N2_2014, H3N2_2019, age_2014_cat)

try7$H3N2_2014 <- as.numeric(try7$H3N2_2014)
try7$H3N2_2019 <- as.numeric(try7$H3N2_2019)

# 2x2 table
table(try7[try7$age_2014_cat==0,]$H3N2_2014, try7[try7$age_2014_cat==0,]$H3N2_2019) 
table(try7[try7$age_2014_cat==1,]$H3N2_2014, try7[try7$age_2014_cat==1,]$H3N2_2019) 

# Bayesian analysis of H3N2 2014 on H3N2 2016, under 5s
H3N2_2014 <- try7[try7$age_2014_cat==0,]$H3N2_2014
H3N2_2019 <- try7[try7$age_2014_cat==0,]$H3N2_2019

logistic.model <- function() {
  # SAMPLING DISTRIBUTION
  for (i in 1:N) {
    logit(pi[i]) <- b[1] + b[2]*H3N2_2014[i];
    H3N2_2019[i] ~ dbin(pi[i],1);  
  }
  
  # PRIORS ON BETAS
  b[1:Nx] ~ dmnorm(mu[1:Nx], tau[1:Nx,1:Nx]) # multivariate normal prior
  
  # Calculate ORs:
  for (l in 1:Nx) {
    OR[l] <- exp(b[l]);
  }
}

N <- length(H3N2_2019) 
Nx <- 2 

mu <- rep(0, Nx)
tau <- diag(0.001, Nx)
data.logistic <- list("N", "Nx", "H3N2_2014", "H3N2_2019", "mu", "tau")
parameters.logistic <-c("b","OR") # Parameters to keep track of
inits.logistic <- function() {list (b= rep(0, Nx))}

set.seed(114011)
logistic.sim_2014_2019_cat0<-jags(data=data.logistic, inits=inits.logistic, parameters.logistic, n.iter=100000,
                                  n.burn=25000, model.file=logistic.model, n.thin=5, n.chains = 3)

print(logistic.sim_2014_2019_cat0, 2)
OR_H3N2_2014_2019_cat0_Bayes <- round(logistic.sim_2014_2019_cat0$BUGSoutput$summary[2,c(5,3,7)], 2)
OR_H3N2_2014_2019_cat0_Bayes

# Bayesian analysis of H3N2 2014 on H3N2 2016, over 5s
H3N2_2014 <- try7[try7$age_2014_cat==1,]$H3N2_2014
H3N2_2019 <- try7[try7$age_2014_cat==1,]$H3N2_2019

logistic.model <- function() {
  # SAMPLING DISTRIBUTION
  for (i in 1:N) {
    logit(pi[i]) <- b[1] + b[2]*H3N2_2014[i];
    H3N2_2019[i] ~ dbin(pi[i],1);  
  }
  
  # PRIORS ON BETAS
  b[1:Nx] ~ dmnorm(mu[1:Nx], tau[1:Nx,1:Nx]) # multivariate normal prior
  
  # Calculate ORs:
  for (l in 1:Nx) {
    OR[l] <- exp(b[l]);
  }
}

N <- length(H3N2_2019) 
Nx <- 2 

mu <- rep(0, Nx)
tau <- diag(0.001, Nx)
data.logistic <- list("N", "Nx", "H3N2_2014", "H3N2_2019", "mu", "tau")
parameters.logistic <-c("b","OR") # Parameters to keep track of
inits.logistic <- function() {list (b= rep(0, Nx))}

set.seed(114011)
logistic.sim_2014_2019_cat1<-jags(data=data.logistic, inits=inits.logistic, parameters.logistic, n.iter=100000,
                                  n.burn=25000, model.file=logistic.model, n.thin=5, n.chains = 3)

print(logistic.sim_2014_2019_cat1, 2)
OR_H3N2_2014_2019_cat1_Bayes <- round(logistic.sim_2014_2019_cat1$BUGSoutput$summary[2,c(5,3,7)], 2)
OR_H3N2_2014_2019_cat1_Bayes

### sensitivity analysis


try7 <- data_orig %>%
  filter(at2014==1 & at2019==1 ) %>%
  dplyr::select(Obs, H3N2_2014, H3N2_2019, age_2014_cat)

try7$H3N2_2014 <- as.numeric(try7$H3N2_2014)
try7$H3N2_2019 <- as.numeric(try7$H3N2_2019)

# 2x2 table
table(try7[try7$age_2014_cat==0,]$H3N2_2014, try7[try7$age_2014_cat==0,]$H3N2_2019) 
table(try7[try7$age_2014_cat==1,]$H3N2_2014, try7[try7$age_2014_cat==1,]$H3N2_2019) 

# Bayesian analysis of H3N2 2014 on H3N2 2016, under 5s
H3N2_2014 <- try7[try7$age_2014_cat==0,]$H3N2_2014
H3N2_2019 <- try7[try7$age_2014_cat==0,]$H3N2_2019

logistic.model <- function() {
  # SAMPLING DISTRIBUTION
  for (i in 1:N) {
    logit(pi[i]) <- b[1] + b[2]*H3N2_2014[i];
    H3N2_2019[i] ~ dbin(pi[i],1);  
  }
  
  # PRIORS ON BETAS
  b[1:Nx] ~ dmnorm(mu[1:Nx], tau[1:Nx,1:Nx]) # multivariate normal prior
  
  # Calculate ORs:
  for (l in 1:Nx) {
    OR[l] <- exp(b[l]);
  }
}

N <- length(H3N2_2019) 
Nx <- 2 

mu <- rep(0, Nx)
tau <- diag(0.001, Nx)
data.logistic <- list("N", "Nx", "H3N2_2014", "H3N2_2019", "mu", "tau")
parameters.logistic <-c("b","OR") # Parameters to keep track of
inits.logistic <- function() {list (b= rep(0, Nx))}

set.seed(114011)
logistic.sim_2014_2019_cat0<-jags(data=data.logistic, inits=inits.logistic, parameters.logistic, n.iter=100000,
                                  n.burn=25000, model.file=logistic.model, n.thin=5, n.chains = 3)

print(logistic.sim_2014_2019_cat0, 2)
OR_H3N2_2014_2019_cat0_Bayes <- round(logistic.sim_2014_2019_cat0$BUGSoutput$summary[2,c(5,3,7)], 2)
OR_H3N2_2014_2019_cat0_Bayes

# Bayesian analysis of H3N2 2014 on H3N2 2016, over 5s
H3N2_2014 <- try7[try7$age_2014_cat==1,]$H3N2_2014
H3N2_2019 <- try7[try7$age_2014_cat==1,]$H3N2_2019

logistic.model <- function() {
  # SAMPLING DISTRIBUTION
  for (i in 1:N) {
    logit(pi[i]) <- b[1] + b[2]*H3N2_2014[i];
    H3N2_2019[i] ~ dbin(pi[i],1);  
  }
  
  # PRIORS ON BETAS
  b[1:Nx] ~ dmnorm(mu[1:Nx], tau[1:Nx,1:Nx]) # multivariate normal prior
  
  # Calculate ORs:
  for (l in 1:Nx) {
    OR[l] <- exp(b[l]);
  }
}

N <- length(H3N2_2019) 
Nx <- 2 

mu <- rep(0, Nx)
tau <- diag(0.001, Nx)
data.logistic <- list("N", "Nx", "H3N2_2014", "H3N2_2019", "mu", "tau")
parameters.logistic <-c("b","OR") # Parameters to keep track of
inits.logistic <- function() {list (b= rep(0, Nx))}

set.seed(114011)
logistic.sim_2014_2019_cat1<-jags(data=data.logistic, inits=inits.logistic, parameters.logistic, n.iter=100000,
                                  n.burn=25000, model.file=logistic.model, n.thin=5, n.chains = 3)

print(logistic.sim_2014_2019_cat1, 2)
OR_H3N2_2014_2019_cat1_Bayes <- round(logistic.sim_2014_2019_cat1$BUGSoutput$summary[2,c(5,3,7)], 2)
OR_H3N2_2014_2019_cat1_Bayes


#2016 to 2019

try7 <- data_orig %>%
  filter(at2016==1 & at2019==1 & H3N2_2017==0 & H3N2_2018==0) %>%
  dplyr::select(Obs, H3N2_2016, H3N2_2019, age_2016_cat)

try7$H3N2_2016 <- as.numeric(try7$H3N2_2016)
try7$H3N2_2019 <- as.numeric(try7$H3N2_2019)

# 2x2 table
table(try7[try7$age_2016_cat==0,]$H3N2_2016, try7[try7$age_2016_cat==0,]$H3N2_2019) 
table(try7[try7$age_2016_cat==1,]$H3N2_2016, try7[try7$age_2016_cat==1,]$H3N2_2019) 

# Bayesian analysis of H3N2 2014 on H3N2 2016, under 5s
H3N2_2016 <- try7[try7$age_2016_cat==0,]$H3N2_2016
H3N2_2019 <- try7[try7$age_2016_cat==0,]$H3N2_2019

logistic.model <- function() {
  # SAMPLING DISTRIBUTION
  for (i in 1:N) {
    logit(pi[i]) <- b[1] + b[2]*H3N2_2016[i];
    H3N2_2019[i] ~ dbin(pi[i],1);  
  }
  
  # PRIORS ON BETAS
  b[1:Nx] ~ dmnorm(mu[1:Nx], tau[1:Nx,1:Nx]) # multivariate normal prior
  
  # Calculate ORs:
  for (l in 1:Nx) {
    OR[l] <- exp(b[l]);
  }
}

N <- length(H3N2_2019) 
Nx <- 2 

mu <- rep(0, Nx)
tau <- diag(0.001, Nx)
data.logistic <- list("N", "Nx", "H3N2_2016", "H3N2_2019", "mu", "tau")
parameters.logistic <-c("b","OR") # Parameters to keep track of
inits.logistic <- function() {list (b= rep(0, Nx))}

set.seed(114011)
logistic.sim_2016_2019_cat0<-jags(data=data.logistic, inits=inits.logistic, parameters.logistic, n.iter=100000,
                                  n.burn=25000, model.file=logistic.model, n.thin=5, n.chains = 3)

print(logistic.sim_2016_2019_cat0, 2)
OR_H3N2_2016_2019_cat0_Bayes <- round(logistic.sim_2016_2019_cat0$BUGSoutput$summary[2,c(5,3,7)], 2)
OR_H3N2_2016_2019_cat0_Bayes

# Bayesian analysis of H3N2 2014 on H3N2 2016, over 5s
H3N2_2016 <- try7[try7$age_2016_cat==1,]$H3N2_2016
H3N2_2019 <- try7[try7$age_2016_cat==1,]$H3N2_2019

logistic.model <- function() {
  # SAMPLING DISTRIBUTION
  for (i in 1:N) {
    logit(pi[i]) <- b[1] + b[2]*H3N2_2016[i];
    H3N2_2019[i] ~ dbin(pi[i],1);  
  }
  
  # PRIORS ON BETAS
  b[1:Nx] ~ dmnorm(mu[1:Nx], tau[1:Nx,1:Nx]) # multivariate normal prior
  
  # Calculate ORs:
  for (l in 1:Nx) {
    OR[l] <- exp(b[l]);
  }
}

N <- length(H3N2_2019) 
Nx <- 2 

mu <- rep(0, Nx)
tau <- diag(0.001, Nx)
data.logistic <- list("N", "Nx", "H3N2_2016", "H3N2_2019", "mu", "tau")
parameters.logistic <-c("b","OR") # Parameters to keep track of
inits.logistic <- function() {list (b= rep(0, Nx))}

set.seed(114011)
logistic.sim_2016_2019_cat1<-jags(data=data.logistic, inits=inits.logistic, parameters.logistic, n.iter=100000,
                                  n.burn=25000, model.file=logistic.model, n.thin=5, n.chains = 3)

print(logistic.sim_2016_2019_cat1, 2)
OR_H3N2_2016_2019_cat1_Bayes <- round(logistic.sim_2016_2019_cat1$BUGSoutput$summary[2,c(5,3,7)], 2)
OR_H3N2_2016_2019_cat1_Bayes

### sensitivity analysis

try7 <- data_orig %>%
  filter(at2016==1 & at2019==1) %>%
  dplyr::select(Obs, H3N2_2016, H3N2_2019, age_2016_cat)

try7$H3N2_2016 <- as.numeric(try7$H3N2_2016)
try7$H3N2_2019 <- as.numeric(try7$H3N2_2019)

# 2x2 table
table(try7[try7$age_2016_cat==0,]$H3N2_2016, try7[try7$age_2016_cat==0,]$H3N2_2019) 
table(try7[try7$age_2016_cat==1,]$H3N2_2016, try7[try7$age_2016_cat==1,]$H3N2_2019) 

# Bayesian analysis of H3N2 2014 on H3N2 2016, under 5s
H3N2_2016 <- try7[try7$age_2016_cat==0,]$H3N2_2016
H3N2_2019 <- try7[try7$age_2016_cat==0,]$H3N2_2019

logistic.model <- function() {
  # SAMPLING DISTRIBUTION
  for (i in 1:N) {
    logit(pi[i]) <- b[1] + b[2]*H3N2_2016[i];
    H3N2_2019[i] ~ dbin(pi[i],1);  
  }
  
  # PRIORS ON BETAS
  b[1:Nx] ~ dmnorm(mu[1:Nx], tau[1:Nx,1:Nx]) # multivariate normal prior
  
  # Calculate ORs:
  for (l in 1:Nx) {
    OR[l] <- exp(b[l]);
  }
}

N <- length(H3N2_2019) 
Nx <- 2 

mu <- rep(0, Nx)
tau <- diag(0.001, Nx)
data.logistic <- list("N", "Nx", "H3N2_2016", "H3N2_2019", "mu", "tau")
parameters.logistic <-c("b","OR") # Parameters to keep track of
inits.logistic <- function() {list (b= rep(0, Nx))}

set.seed(114011)
logistic.sim_2016_2019_cat0<-jags(data=data.logistic, inits=inits.logistic, parameters.logistic, n.iter=100000,
                                  n.burn=25000, model.file=logistic.model, n.thin=5, n.chains = 3)

print(logistic.sim_2016_2019_cat0, 2)
OR_H3N2_2016_2019_cat0_Bayes <- round(logistic.sim_2016_2019_cat0$BUGSoutput$summary[2,c(5,3,7)], 2)
OR_H3N2_2016_2019_cat0_Bayes

# Bayesian analysis of H3N2 2014 on H3N2 2016, over 5s
H3N2_2016 <- try7[try7$age_2016_cat==1,]$H3N2_2016
H3N2_2019 <- try7[try7$age_2016_cat==1,]$H3N2_2019

logistic.model <- function() {
  # SAMPLING DISTRIBUTION
  for (i in 1:N) {
    logit(pi[i]) <- b[1] + b[2]*H3N2_2016[i];
    H3N2_2019[i] ~ dbin(pi[i],1);  
  }
  
  # PRIORS ON BETAS
  b[1:Nx] ~ dmnorm(mu[1:Nx], tau[1:Nx,1:Nx]) # multivariate normal prior
  
  # Calculate ORs:
  for (l in 1:Nx) {
    OR[l] <- exp(b[l]);
  }
}

N <- length(H3N2_2019) 
Nx <- 2 

mu <- rep(0, Nx)
tau <- diag(0.001, Nx)
data.logistic <- list("N", "Nx", "H3N2_2016", "H3N2_2019", "mu", "tau")
parameters.logistic <-c("b","OR") # Parameters to keep track of
inits.logistic <- function() {list (b= rep(0, Nx))}

set.seed(114011)
logistic.sim_2016_2019_cat1<-jags(data=data.logistic, inits=inits.logistic, parameters.logistic, n.iter=100000,
                                  n.burn=25000, model.file=logistic.model, n.thin=5, n.chains = 3)

print(logistic.sim_2016_2019_cat1, 2)
OR_H3N2_2016_2019_cat1_Bayes <- round(logistic.sim_2016_2019_cat1$BUGSoutput$summary[2,c(5,3,7)], 2)
OR_H3N2_2016_2019_cat1_Bayes



