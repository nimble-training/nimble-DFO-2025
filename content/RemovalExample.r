library(tidyverse)
library(nimble)
library(coda)

remdf <- read.csv("data/removal_data.csv")
catch <- remdf %>% 
  pivot_wider(id_cols = c(Site, Station, Year, Sweep), names_from = Species, 
              values_from = Fish.ID, values_fn = ~ length(.x), values_fill = 0) %>%
              pivot_longer(unique(remdf$Species), names_to = "Species", values_to = "Catch")

catch <- catch %>% group_by(Site, Station, Year, Species) %>%
  arrange(Sweep) %>%
  mutate(prevCatch = cumsum(Catch)-Catch) %>%
  ungroup() %>%
  mutate(Station = as.numeric(factor(Station))) %>%
  mutate(year = Year - min(Year) + 1) %>%
  filter(Species %in% c("AS", "BT")) %>%
  mutate(sppid = as.numeric(factor(Species))) %>%
  mutate(yearid = as.numeric(factor(year)))


## Fill in the nimble code to get started.
model_code <- nimbleCode({
  for( j in 1:nspp ){
    intercept[j] ~ dnorm(0, sd = 3)
    pop_trend[j] ~ dnorm(0, sd = 3) 
    logitp_0[j] ~ dnorm(0, sd = 3)
  }
  sd_re ~ dunif(0, 10)

  for( i in 1:nyear ){
    for( j in 1:nspp ){
      log(lambda[i,j]) <- intercept[j] + pop_trend[j]*year[i]  ## year species average abundance.
      logitp[i,j] ~ dnorm(logitp_0[j], sd = sd_re) ## year species detectability.
      for( s in 1:nstation ){
        N[i,j,s] ~ dpois(lambda[i,j]) ## Improve this next!
      }
    }
  }
  
  for( i in 1:nobs ){
    logit(p[i]) <- logitp[obsyear[i], spp[i]]
    caught[i] ~ dbinom(size = N[obsyear[i], spp[i], station[i] ] - total[i], prob = p[i])
  }
})

data <- list(caught = catch$Catch)
constants <- list(nobs = nrow(catch), total = catch$prevCatch, 
                  nstation = max(catch$Station), 
                  station = catch$Station, spp = catch$sppid, 
                  nspp = 2, nyear = length(unique(catch$year)),
                  obsyear = catch$yearid, 
                  year = unique(catch$year)-1)

maxn <- catch %>% 
  group_by(Station, sppid, yearid) %>% 
  summarize(n = sum(Catch), .groups = "drop")
N_init <- array(NA, c(constants$nyear, constants$nspp, constants$nstation))
for(i in 1:nrow(maxn))
  N_init[maxn$yearid[i], maxn$sppid[i], maxn$Station[i]] <- maxn$n[i] + rpois(1, 5)

inits <- function(){
  list(N = N_init, 
       lambda = matrix(50, nrow=constants$nyear, ncol = constants$nspp), intercept = rnorm(2,0,1),
       pop_trend = c(0,0), logitp_0 = c(0,0), sd_re = runif(1,0.5,2))
}

## Fill in to start a model
rmodel <- nimbleModel(model_code, data = data, constants = constants, inits = inits())

## Compile model
cmodel <- compileNimble(rmodel)

## configure + compile MCMC
conf <- configureMCMC(rmodel) ## How many samplers on N?
mcmc <- buildMCMC(conf)
cmcmc <- compileNimble(mcmc)

samples <- runMCMC(cmcmc, niter = 10000, nburnin = 2500, 
                    inits = inits, nchains = 3, samplesAsCodaMCMC = TRUE)

plot(samples[, c("pop_trend[1]", "pop_trend[2]")])

avgcatch <- catch %>% 
  group_by(Station, Species, Year) %>% 
  summarize(n = sum(Catch), .groups = "drop") %>%
  group_by(Species, Year) %>%
  summarize(count = mean(n), .groups = "drop")

ggplot(data = avgcatch, aes(x = Year, y = count, colour = Species)) +  
  geom_point() + 
  theme_bw() + 
  geom_line()