library(tidyverse)
library(readxl)
library(nimble)
library(coda)

sockeye <- read_excel("data/Production Data_Detailed Format.xlsx")
sockeye <- sockeye %>% 
  filter(!is.na(total_broodyr_spawners)) %>%  
  group_by(production_stock_name, broodyr, total_broodyr_spawners) %>%
  summarize(recruits=sum(num_recruits), .groups = "drop") %>%
  group_by(production_stock_name) %>%
  mutate(scale = 10^floor(log10(max(total_broodyr_spawners, na.rm = TRUE)))) %>%
  ungroup() %>%
  rename(broodyear = broodyr, spawners = total_broodyr_spawners) %>%
  mutate(logRS = log(recruits) - log(spawners)) %>%
  filter(!is.na(spawners) & !is.na(recruits)) %>%
  mutate(spawners_scaled = spawners/scale) %>%
  mutate(recruits_scaled = recruits/scale)

# stockid <- "Stellako"
# stockid <- "Birkenhead"
stockid <- "Chilko"
sox <- sockeye %>%  
  filter(production_stock_name == stockid)

## Fit linear model:
fit.lm <- lm(logRS ~ spawners_scaled, data = sox)
logalpha <- coef(fit.lm)[1]
E <- abs(logalpha/coef(fit.lm)[2])
sigma <- sigma(fit.lm)

plot(sox$spawners_scaled, sox$recruits_scaled)
lines(seq(0, 5, 0.1), seq(0, 5, 0.1)*exp(logalpha*(1-seq(0, 5, 0.1)/E)), col = 'red')

## Confidence Intervals:
## Parametric bootstrap intervals:
Esim <- logalphasim <- betasim <- NULL
for(i in 1:10000) {
  logRS.sim <- rnorm(nrow(sox), logalpha*(1-sox$spawners_scaled/E), sd = sigma)  
  new.fit <- lm(logRS.sim~spawners_scaled, data= sox)
  logalphanew <- coef(new.fit)[1]
  logalphasim <- c(logalphasim, logalphanew)
  betasim <- c(betasim, coef(new.fit)[2])
  Esim <- c(Esim, abs(logalphanew/coef(new.fit)[2]))
}
quantile(logalphasim, c(0.025, 0.975))
quantile(betasim, c(0.025, 0.975))
quantile(Esim, c(0.025, 0.975))
confint(fit.lm)


## Now do a Bayesian analysis in R using NIMBLE:

## Fill in the nimble code to get started.
model_code <- nimbleCode({
  ## Setup Priors
  logalpha ~ 
  sigma ~ 
  E ~ 

  ## Likelihood:
  for( i in 1:nobs ){
    ## E(logRS) = logalpha * (1-S/E)
    logRS[i] ~
  }
})

data <- list(logRS = sox$logRS)
constants <- list(nobs = nrow(sox), S = sox$spawners/spscale)
inits <- function(){
  list(sigma = runif(1,0.1,3), logalpha = rnorm(1, 1.5, 3), logE = rnorm(1, 4, 2))
}

## Fill in to start a model
rmodel <- nimbleModel()

## Compile model
cmodel <- 

## configure + compile MCMC
conf <- 
mcmc <- 
cmcmc <- 

samples <- runMCMC(cmcmc, niter = 5000, nburnin = 1000, 
                    inits = inits, nchains = 3, samplesAsCodaMCMC = TRUE)

