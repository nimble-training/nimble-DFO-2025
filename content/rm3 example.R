## R/M3 model---- CSAM3
#"alpha + beta*ST[i] + gamma*TI[i]*(1-TI[i]/delta) + epsilon*CO[i]"

# ST = larval abundance
# TI = timing of sea ice retreat (day of the year)
# CO = fall condition of age 1 and 2, male and female capelin

m.RM3 = '
model {
# 1. Likelihood
for (i in 1:N) {
#recruitment
mu[i] <- alpha + beta*ST[i] + gamma*TI[i]*(1-TI[i]/delta) + epsilon*CO[i]
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
                   TI=c(df3$tice, TIpred), #made up - need new data
                   CO=c(df3$meanCond_lag, COpred),
                   N = nrow(df3) + num_forecasts)
mod# model_data <- list(N2 = c(pred3[1:(length(pred3)-1)]), 
#                    ST=c(df3$surface_tows_lag2), #from capelin_larval_indices - see df_norm
#                    TI=c(df3$tice), #made up - need new data
#                    CO=c(df3$meanCond_lag),
#                    N = nrow(df3))

# file.for.paul <- data.frame(year=c(df3$year,2025, 2026), biomass = c(pred3, NA, NA), ST = c(df3$surface_tows_lag2, STpred), CO = c(df3$meanCond_lag, COpred), PS = c(df3$ps_meanTot_lag2, PSpred), tice=c(df3$tice, TIpred))
# write.csv(file.for.paul, file="data-2025/capelin_example.csv", row.names=FALSE)


run_RM3 <- jags(data=model_data,
                parameters.to.save = c('mu', 'sigma', 'N2', 'N2_new', 'alpha', 'beta', 'gamma', 'delta', 'epsilon', 'Fit', 'FitNew', 'Res', 'PRes', 'expY', 'D', "log_lik"),
                model.file = textConnection(m.RM3))

#UPDATE WITH MORE BURN INS
run_RM3 <-update(run_RM3, n.iter = 1000000, n.thin = 100, n.burnin = 1)
#save(run_RM3, file=paste("Bayesian/", filepath, "/runRM3.Rdata", sep=""))


##R/M3-plot 5----
# DIAGNOSTICS
print(run_RM3, intervals=c(0.025, 0.975), digits = 3)
out <- run_RM3$BUGSoutput 
out$mean
#overdispersion - values close to 0.5 indicate a good fit pg. 77 Zuur et al. 2013 Beginners guide to GLM and GLMM
mean(out$sims.list$FitNew > out$sims.list$Fit)
# mean = 0.546333

# Asess mixing of chains to see if one MCMC goes badly Zuur et al. 2013, pg 83
vars <- c('alpha', 'beta', 'gamma', 'delta', 'epsilon')
filepath <- paste0(filepath_gen, "/rm_3")

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
plot(test$TI, test$EI, xlab = "Tice", ylab = "Pearson resids")
plot(test$CO, test$EI, xlab = "Condition", ylab = "Pearson resids")
par(mfrow = c(1,1))
dev.off()

# CREDIBLE AND PREDICITON INTERVALS
#generate credible intervals for time series using mu
y_pred = run_RM3$BUGSoutput$sims.list$mu

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
y_new = run_RM3$BUGSoutput$sims.list$N2_new
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

# Feb 2023 - I don't think we can do this
# section this year as no run for last year
# Also questions about data for this year.
# Commenting out for now, but revisit in 2024
# this seems like a good place to look at the change from 2000 to 2021...
# 2020 is 2 less than max of xaxis and 2021 is 1 less than max of xaxis
# min.pi.2020<-mariano.combined2.pi[1,2]
# max.pi.2020<-mariano.combined2.pi[1,4]
# 
# prob.2021.gt.2020<-length(which(y_new[,19]>max.pi.2020))/length(y_new[,19])
# prob.2021.lt.2020<-length(which(y_new[,19]<min.pi.2020))/length(y_new[,19])
# prob.2021.eq.2020<-1 - prob.2021.gt.2020 - prob.2021.lt.2020
# 
# probs.change<-c(prob.2021.lt.2020, prob.2021.eq.2020, prob.2021.gt.2020)
# labels.change<-c("Prob. 2021 less 2020", "Prob. 2021 equals 2020", "Prob. 2021 greater 2020")
# prob.table.2020.v.2021<-data.frame(scenario=labels.change, probabilities=probs.change)
# write.csv(prob.table.2020.v.2021, paste0("Bayesian/", filepath, "/probabilities of 2020 vs 2021.csv"), row.names = FALSE)

## RM_3-Results----
dic_RM3 <- dic.samples(run_RM3$model, n.iter=1000, type="pD")
dic_RM3sum <- sum(dic_RM3$deviance)


## RM_3-Results----
# pD version of DIC
dic_RM3.pd.samples <- dic.samples(run_RM3$model, type="pD", n.iter = 1000000, thin=100, burnin = 1)
# str(x)
dic_RM3_sum.pd <- sum(dic_RM3.pd.samples$deviance)
dic_RM3_sumP.pd <- sum(dic_RM3.pd.samples$penalty)
dic_RM3.pd <- dic_RM3_sum.pd + dic_RM3_sumP.pd
dic_RM3.pd

# # popt version of DIC
# dic_RM3.popt.samples <- dic.samples(run_RM3$model, type="popt", n.iter = 20000000, thin=1000, burnin = 1)
# # str(x)
# dic_RM3_sum.popt <- sum(dic_RM3.popt.samples$deviance)
# dic_RM3_sumP.popt <- sum(dic_RM3.popt.samples$penalty)
# dic_RM3.popt <- dic_RM3_sum.popt + dic_RM3_sumP.popt
# dic_RM3.popt

RM3_r2m <- rsq_bayes(ypred = y_pred, out = run_RM3)
RM3_r2 <- c(median(RM3_r2m), sd(RM3_r2m))

# Zuur pg 85: note that MCMCSupportHighstatV2.R (line 114) says this is to look at ACF - I think that this is wrong. 

MyBUGSHist1(out, vars, transform = "yes")
pdf(file = paste0("Bayesian/", filepath, "/posteriors.pdf"), width=10, height=8)
MyBUGSHist1(out, vars, transform = "yes")
dev.off()

#plot credible and prediction intervals
test<-plotCredInt(df3, yaxis = yaxis1, xaxis=xaxis,
                  ylab = ylab1, 
                  y_line = y_med, ci=ci_df3, dpi=pi_df3_80[c(1,5),], 
                  model = txt, x = 2010, y = 8, type = "CI")
ggsave(paste0("Bayesian/", filepath, "/credInt.png"), width=10, height=8, units="in")
ggsave(paste0("Bayesian/", filepath, "/credInt.pdf"), width=10, height=8, units="in")


#plot credible and prediction intervals
test<-plotCredInt_for3(df3, yaxis = yaxis1, xaxis=xaxis,
                       ylab = ylab1, 
                       y_line = y_med, ci=ci_df3, dpi=pi_df3_80[c(1,3,5),], 
                       model = txt, x = 2010, y = 8, type = "CI")
ggsave(paste0("Bayesian/", filepath, "/newcredInt.png"), width=10, height=8, units="in")
ggsave(paste0("Bayesian/", filepath, "/newcredInt.pdf"), width=10, height=8, units="in")

# try to make a linear version of this plot...
y_med_exp<-exp(y_med)
ci_df3_exp<-exp(ci_df3)
mid50_df3_exp<-exp(mid50_df3)
dpi_exp<-exp(pi_df3_80)
yaxis.exp<-"biomass_med"
type<-"CI"
ylab.exp<-"Capelin biomass index (ktonnes)"

# this line isn't working because lwr, and upr should be pointing to independent vectors or components of data frames. 
# commenting this code out because it isn't working and I think that I deal with it far below... around line 2606
# linear.plot.data<-data.frame(med.val=y_med_exp, lwr.ci=lwr, upr.ci=upr, year=c(df3$year,c(2020,2021)))
# # same graph on a linear scale...
# all.retro.linear<-all.retro %>% mutate(y_med_l = exp(y_med), lwr_l=exp(lwr), upr_l=exp(upr))
# 
# # this is working, but I'd like to automate the last part (ribbons...)
# q <- ggplot(data=linear.plot.data, aes(x=year, y=med.val)) + ylab("Capelin biomass (ktonnes)") + xlab("Year") + geom_point() + geom_line() + geom_ribbon(aes(ymin=lwr.ci, ymax=upr.ci), alpha=0.2) + theme_classic()
# pdf(file=paste("Bayesian/", filepath, "/retro_plot_linear.pdf", sep=""), width=10, height=8)
# q

num.years<-length(xaxis)
# create initial plot
q <- ggplot()  
# plot the error bars of the fitted portion of model
q <- q + geom_ribbon(aes(x = xaxis[1:(length(xaxis)-2)], 
                         ymax = ci_df3_exp[2,1:(dim(ci_df3_exp)[2]-2)], 
                         ymin = ci_df3_exp[1,1:(dim(ci_df3_exp)[2]-2)]),
                     alpha = 0.5, 
                     fill = "grey70")
q <- q + geom_ribbon(aes(x = xaxis[1:(length(xaxis)-2)], 
                         ymax = mid50_df3_exp[3,1:(dim(mid50_df3_exp)[2]-2)], 
                         ymin = mid50_df3_exp[1,1:(dim(mid50_df3_exp)[2]-2)]),
                     alpha = 0.5, 
                     fill = "grey50")
pi_n <- dpi_exp[, (ncol(dpi_exp)-1):ncol(dpi_exp)]
# plot the error bars of the forecast portion of the model
q <- q + geom_ribbon(aes(x = c(xaxis[(num.years-1):num.years]),
                         ymax = pi_n[5, ],
                         ymin = pi_n[1, ]),
                     fill = "grey55")
# plot the error bars of the forecast portion of the model
q <- q + geom_ribbon(aes(x = c(xaxis[(num.years-1):num.years]),
                         ymax = pi_n[4, ],
                         ymin = pi_n[2, ]),
                     fill = "grey35")
# Label and format the axes
q <- q + geom_point(data = df3,
                    aes_string(y = yaxis.exp, x = "year"),
                    shape = 16, 
                    size = 1.5)
# If no type, add error bars for the observed values
if(!is.na(type)){
  q <- q + geom_errorbar(data = df3, 
                         width = 0.3, 
                         colour = "black", 
                         aes(x = year, 
                             min=bm_lci, 
                             ymax=bm_uci))
} else if (is.na(type)) {
  q
}
# Label the axes
q <- q + xlab("Year / Année") + ylab("Capelin biomass index (kt)\nIndice de la biomasse du capelan (kt)")
# Insert line for the median?
q <- q + geom_line(aes(x = xaxis[1:(length(xaxis)-2)], 
                       y = y_med_exp[1:(length(y_med_exp)-2)]))
# insert remainder of the line for the median
q <- q + geom_line(aes(x = xaxis[(length(xaxis)-1):length(xaxis)], 
                       y = y_med_exp[(length(y_med_exp)-1):length(y_med_exp)]))									   
# format the plot
q <- q + theme_bw(base_size = 30) + 
  theme(plot.margin = unit(c(0.5, 1, 0.5, 0.5), "cm"))
q <- q+ scale_y_continuous(lim=c(0,1800), expand=c(0,0))
q
ggsave(paste0("Bayesian/", filepath, "/credInt_linear.png"), width=12, height=8, units="in")
ggsave(paste0("Bayesian/", filepath, "/credInt_linear.pdf"), width=12, height=8, units="in")

# add in LRP line
q <- q + geom_hline(yintercept=155, col="darkred")
q
ggsave(paste0("Bayesian/", filepath, "/credInt_linear_lrp.png"), width=12, height=8, units="in")

# data.frame for Erin Dunne
output.forecast.model<-data.frame(year=xaxis, lower.fit.025=c(round(ci_df3_exp[1,1:(dim(ci_df3_exp)[2]-num_forecasts)],1), rep(NA, num_forecasts)), upper.fit.0975=c(round(ci_df3_exp[1,1:(dim(ci_df3_exp)[2]-num_forecasts)],1), rep(NA, num_forecasts)), lower.prediction.1=c(rep(NA, length(xaxis) - num_forecasts), round(pi_n[1,],1)), upper.prediction.9=c(rep(NA,length(xaxis) - num_forecasts), round(pi_n[2,],1)), observed.biomass=c(df3$biomass_med, rep(NA,num_forecasts)), biomass.lci=c(df3$bm_lci, rep(NA, num_forecasts)), biomass.uci=c(df3$bm_uci, rep(NA, num_forecasts)))
write.csv(output.forecast.model, file="output_values_forecast_model.csv", row.names = FALSE)

# output by parameter
OUT1 <- MyBUGSOutput(out, vars)
print(OUT1, digits =5)
write.csv(as.data.frame(OUT1), paste0("Bayesian/", filepath, "/params.csv"))

#csm3 <- out
# plot posterior against expected distribution
alpha <- posterior_fig(out$sims.list$alpha)
beta <- posterior_fig(out$sims.list$beta)
gamma <- posterior_fig1(out$sims.list$gamma, transform = "yes", parm = "slope")
delta <- posterior_fig1(out$sims.list$delta, transform = "yes", parm = "width")
epsilon <- posterior_fig(out$sims.list$epsilon)


#alpha - INT
priormean <- 0
priorsd <- 100
prior <- rnorm(n = 10000, mean = priormean, sd = priorsd)
limits <- c(min(alpha$df) * 0.99, max(alpha$df) * 0.99)
x_label <- "Intercept"
bin_1 <- (limits[2] - limits[1])/100

p1 <- postPriors(df = alpha$df, df2 = prior, df3 = alpha$df_cred, limits, x_label, priormean, priorsd, by_bin = bin_1)

#beta - larval abundance
priorsd <- 100
limits <- c(min(beta$df) * 0.99, max(beta$df) * 0.99)
x_label <- "Larval abundance"
bin_1 <- (limits[2] - limits[1])/100

p2 <- postPriors(df = beta$df, df2 = prior, df3 = beta$df_cred, limits, x_label, priormean, priorsd, by_bin = bin_1)

#gamma - tice-max rate increase
#prior <- runif(n = 100000, min = 0, max = 100)/100
prior <- runif(n = 100000, min = 0, max = 100)
limits <- c(0, max(gamma$df) * 0.99)
x_label <- expression(paste(t[italic(ice)], "MRI"))
bin_1 <- (limits[2] - limits[1])/100

p3 <- postPriors(df = gamma$df, df2 = prior, df3 = gamma$df_cred, limits, x_label, priormean, priorsd, by_bin = bin_1, vline = "no")
#p3 <- p3 + geom_hline(yintercept = 1, size = 1)
p3 <- p3 + geom_segment(aes(x = 0, y = 1, xend = 1, yend = 1), size = 1)

#delta - tice-width
priormean <- 11.5
priorsd <- 5.7
prior <- rgamma(n = 10000, shape = priormean, rate = priorsd)*100
limits <- c(min(delta$df) * 0.99, max(delta$df) * 0.99)
x_label <- expression(paste(t[italic(ice)], "width"))
bin_1 <- (limits[2] - limits[1])/100

p4 <- postPriors(df = delta$df, df2 = prior, df3 = delta$df_cred, limits, x_label, priormean, priorsd, by_bin = bin_1)

#epsilon - condition
priormean <- 0
priorsd <- 100
prior <- rnorm(n = 10000, mean = priormean, sd = priorsd)
limits <- c(min(epsilon$df) * 0.99, max(epsilon$df) * 0.99)
#limits <- c(quantile(epsilon$df, 0.05), quantile(epsilon$df, 0.8))
x_label <- "Condition"
bin_1 <- (limits[2] - limits[1])/100

quantile(epsilon$df, c(0.01, 0.1, 0.2, 0.3, 0.4, 0.50, 0.6, 0.7, 0.8, 0.9, 0.99))

p5 <- postPriors(df = epsilon$df, df2 = prior, df3 = epsilon$df_cred, limits, x_label, priormean, priorsd, by_bin = bin_1)

mm <- cowplot::plot_grid(p1, p2, p3, p4, p5, labels = c("A", "B", "C", "D", "E"), ncol=2)

ggsave(paste0("Bayesian/", filepath, "/priorPost.png"), width=10, height=8, units="in")

ggsave(paste0("Bayesian/", filepath, "/priorPost.pdf"), width=10, height=8, units="in")


# need to remove the large files:
rm(run_RM3)
rm(dic_RM3.pd.samples)
# rm(dic_RM3.popt.samples)
gc()