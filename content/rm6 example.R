## R/M6 model---- CSAM6
#"alpha + beta*ST[i] + theta*PS[i] + gamma*TI[i]*(1-TI[i]/delta) + epsilon*CO[i]"

m.RM6 = '
model {
# 1. Likelihood
for (i in 1:N) {
#recruitment
mu[i] <- alpha + beta*ST[i] + theta*PS[i] + gamma*TI[i]*(1-TI[i]/delta) + epsilon*CO[i]
N2[i] ~ dnorm(mu[i], sigma^-2)
N2_new[i] ~ dnorm(mu[i], sigma^-2) # #### ADB: This is simulated data   
log_lik[i] <- logdensity.norm(N2[i], mu[i], sigma^-2)

# 3. Discrepancy measures
expY[i] <- mu[i]
varY[i] <- sigma^2
Res[i] <- (N2[i] - expY[i])
PRes[i] <- (N2[i] - expY[i]) / sigma
PResNew[i] <- (N2_new[i] - expY[i]) / sigma
#Squared residuals
D[i] <- pow(PRes[i], 2) #SSQ
DNew[i] <- pow(PResNew[i], 2)
#CD[i] <- 
}
#Sum of squared Pearson residuals:
Fit <- sum(D[1:N]) # look at overdispersion
FitNew <- sum(DNew[1:N])

# 2. Priors
alpha ~ dnorm(0, 100^-2) # Int
beta ~ dnorm(0, 100^-2) # larval abun
theta ~ dnorm(0, 100^-2) # pseudocalanus
gamma ~ dunif(0, 100) # tice - max rate increase
delta ~ dgamma(11.5, 5.7) #tice-width
epsilon ~ dnorm(0, 100^-2) # condition
sigma ~ dunif(0, 100) 
}'

#gamma and delta based on Bolker pg 132 - Fig4.13 - trying for an uniformative alpha
#delta: shape(a) is mean^2/var; we used 90 days based on Ales original #work; scale(s) equal Var/mean
# gamma,delta, sigma: uninformative for condition
mean(df3$tice)^2/var(df3$tice)
var(df3$tice)/mean(df3$tice)

model_data <- list(N2 = c(pred3, rep(NA, num_forecasts)),
                   ST=c(df3$surface_tows_lag2, STpred), #from capelin_larval_indices - see df_norm
                   PS=c(df3$ps_meanTot_lag2, PSpred),
                   TI=c(df3$tice, TIpred), #made up - need new data
                   CO=c(df3$meanCond_lag, COpred),
                   N = nrow(df3) + num_forecasts)
# model_data <- list(N2 = c(pred3[1:(length(pred3)-1)]), 
#                    ST=c(df3$surface_tows_lag2), #from capelin_larval_indices - see df_norm
#                    TI=c(df3$tice), #made up - need new data
#                    CO=c(df3$meanCond_lag),
#                    N = nrow(df3))

run_RM6 <- jags(data=model_data,
                parameters.to.save = c('mu', 'sigma', 'N2', 'N2_new', 'alpha', 'beta', 'theta', 'gamma', 'delta', 'epsilon', 'Fit', 'FitNew', 'Res', 'PRes', 'expY', 'D', "log_lik"),
                model.file = textConnection(m.RM6))

#UPDATE WITH MORE BURN INS
run_RM6 <-update(run_RM6, n.iter = 1000000, n.thin = 100, n.burnin = 1)
#save(run_RM3, file=paste("Bayesian/", filepath, "/runRM3.Rdata", sep=""))


##R/M6-plot 6----
# DIAGNOSTICS
print(run_RM6, intervals=c(0.025, 0.975), digits = 3)
out <- run_RM6$BUGSoutput 
out$mean
#overdispersion - values close to 0.5 indicate a good fit pg. 77 Zuur et al. 2013 Beginners guide to GLM and GLMM
mean(out$sims.list$FitNew > out$sims.list$Fit)
# mean = 0.546333

# Asess mixing of chains to see if one MCMC goes badly Zuur et al. 2013, pg 83
vars <- c('alpha', 'beta', 'theta', 'gamma', 'delta', 'epsilon')
filepath <- paste0(filepath_gen, "/rm_6")

MyBUGSChains(out, vars)
pdf(file = paste0("Bayesian/", filepath, "/chains.pdf"), width=10, height=8)
MyBUGSChains(out, vars)
dev.off()

#autocorrelation - this looks good!!!!
MyBUGSACF(out, vars)
pdf(file = paste0("Bayesian/", filepath, "/auto_corr.pdf"), width=10, height=8)
MyBUGSACF(out, vars)
dev.off()

# Model Validation (see Zuuer et al. 2013, pg 77-79 for options for calculating Pearson residuals)
# Residual diagnostics
E1 <- out$mean$PRes # Pearson resids
F1 <- out$mean$expY # Expected values
N2 <- out$mean$N2   # N2 - observed values? Why do these fill in the NAs of df$ln_biomass_med???
D <- out$mean$D     # this is SSQ - but i'm looking for Cook'sD

pdf(paste0("Bayesian/", filepath, "/fit_obs.pdf"))
par(mfrow = c(2,2), mar = c(5,5,2,2))
plot(x=F1, y = E1, xlab = "Fitted values", ylab = "Pearson residuals")
abline(h = 0, lty = 2)
# The below is right but is the not right code.  The code is in ONeNote
#plot(D, type = "h", xlab = "Observation", ylab = "Cook distance")
plot(y = N2, x = F1, xlab = "Fitted values", ylab = "Observed data")
abline(coef = c(0,1), lty = 2)
par(mfrow = c(1,1))
dev.off()

# Residuals v covariates Zuur et al. 2013, pg 59-60: look for no patterns; patterns may indicated non-linear
pdf(paste0("Bayesian/", filepath, "/resid_covar.pdf"))
par(mfrow = c(2,2), mar = c(5,5,2,2))
MyVar <- c("tice.std", "meandCond_lag.std")
test <- as.data.frame(model_data)
test <- cbind(test, E1)
#Myxyplot(test, MyVar, "E1")
plot(test$ST, test$EI, xlab = "ST", ylab = "Pearson resids")
plot(test$PS, test$EI, xlab = "PS", ylab = "Pearson resids")
plot(test$TI, test$EI, xlab = "Tice", ylab = "Pearson resids")
plot(test$CO, test$EI, xlab = "Condition", ylab = "Pearson resids")
par(mfrow = c(1,1))
dev.off()

# CREDIBLE AND PREDICITON INTERVALS
#generate credible intervals for time series using mu
y_pred = run_RM6$BUGSoutput$sims.list$mu

#look at output
y_med = apply(y_pred,2,'median')
write.csv(y_med, paste0("Bayesian/", filepath, "/y_med.csv"))
ci_df3 <- apply(y_pred,2,'quantile', c(0.025, 0.975))
write.csv(ci_df3, paste0("Bayesian/", filepath, "/ci.csv"))
mid50_df3 <- apply(y_pred, 2, quantile, c(0.25, 0.5, 0.75))

# mariano's values... 
mariano.ci<-apply(y_pred, 2, quantile, c(0.025, 0.5, 0.975))
mariano.mean<-apply(y_pred, 2, mean, na.rm=TRUE)
mariano.combined<-rbind(mariano.ci, mariano.mean)
mariano.combined2<-t(rbind(xaxis, mariano.combined))
colnames(mariano.combined2)<-c("year", "2.5% percentile CI", "50% percentile CI", "97.5% percentile CI", "mean CI")
mariano.combined2<-as.data.frame(mariano.combined2)
mariano.combined2<-mariano.combined2 %>% filter(year<2023)
write.csv(mariano.combined2, paste0("Bayesian/", filepath, "/CI for mariano.csv"), row.names = FALSE)

#generate prediciton intevals using N2_new
y_new = run_RM6$BUGSoutput$sims.list$N2_new
#pi_df3 <- apply(y_new,2,'quantile', c(0.025, 0.975))
#pi_df3_80_p_med <- apply(y_new,2,'quantile', c(0.1, 0.25, 0.5, 0.75, 0.9))
pi_df3_80 <- apply(y_new,2,'quantile', c(0.1, 0.25, 0.5, 0.75, 0.9))
write.csv(pi_df3_80[, (ncol(pi_df3_80)-1):ncol(pi_df3_80)], paste0("Bayesian/", filepath, "/pi.csv"))

# proportion of outcomes above/below LRP
pred.year <- y_new[,(ncol(y_new)-1)]
tot.scenarios <- length(pred.year)
num.above.lrp <- length(which(exp(pred.year)>lrp))
num.below.lrp <- length(which(exp(pred.year)<lrp))

lrp.comparison <- data.frame(prop.above = num.above.lrp/tot.scenarios, prop.below = num.below.lrp/tot.scenarios)

write.csv(lrp.comparison, paste0("Bayesian/", filepath, "/lrp outcomes.csv"), row.names = FALSE)

# now need to do the Prediction intervals
mariano.pi<-apply(y_new, 2, quantile, c(0.1, 0.5, 0.9))
mariano.mean.pi<-apply(y_new, 2, mean, na.rm=TRUE)
mariano.combined.pi<-rbind(mariano.pi, mariano.mean.pi)
mariano.combined2.pi<-t(rbind(xaxis, mariano.combined.pi))
colnames(mariano.combined2.pi)<-c("year", "10% percentile PI", "50% percentile PI", "90% percentile PI", "mean PI")
mariano.combined2.pi<-as.data.frame(mariano.combined2.pi)
mariano.combined2.pi<-mariano.combined2.pi %>% filter(year>2022)
write.csv(mariano.combined2.pi, paste0("Bayesian/", filepath, "/PI for mariano.csv"), row.names = FALSE)