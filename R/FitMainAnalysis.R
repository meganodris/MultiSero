#----- Script to fit to arbovirus serology data -----#
library(rstan)
library(cmdstanr)
library(ggplot2)
library(cowplot)
library(mvtnorm)
library(loo)
library(LaplacesDemon)
library(qs)
library(data.table)
library(bayesplot)
library(posterior)
library(here)


#--- source functions
source(here('R', 'RFunctions.R'))


#--- read in data
df <- read.csv(here('data', "arbovirus_serology.csv"))
pathogens <- colnames(df)[3:13]

# if using DENV ELISA data, remove NA's
df <- df[!is.na(df$DENV_ELISA), ]

# specify which pathogens we assume to be present
present <- c('DENV1','CHIKV','JEV')
nonpres <- pathogens[!pathogens %in% present]
pathogens <- c(present, nonpres) # reorder


#--- compile data for model fitting
data <- list()
data$y <- cbind(log(df[,c(present,nonpres)])) # titer data
data$pres <- c(rep(1,length(present)), rep(0, length(nonpres))) # index present pathogens
data$N <- nrow(data$y) # N individuals
data$nP <- ncol(data$y) # N antigens
data$nPp <- sum(data$pres) # N present pathogens
data$nA <- 7 # N age groups
data$nL <- 5 # N locations
data$ageG <- df$ageG # individual age groups
data$loc <- df$locID # individual locations
data$NperL <- as.vector(table(df$locID)) # N individuals per location
data$NperLA <- t(table(df$locID, df$ageG)) # N individuals per location & age group
data$ageProp <- as.matrix(table(df$locID, df$ageG) / data$NperL) # age proportions per location
data$indELISA <- which(pathogens=='DENV_ELISA') # ELISA index
data$infM <- inf_matrix(data$nP, pres=data$pres) # infection status combinations
data$nC <- nrow(data$infM) # N status combinations

# some indexes for model fitting
npos <- rowSums(data$infM)
wpos <- wneg <- matrix(0, ncol=data$nP, nrow=data$nC)
for(c in 1:nrow(data$infM)) for(p in 1:data$nP){
  if(npos[c]>0) wpos[c,1:npos[c]] <- which(data$infM[c,]==1)
  if(npos[c]<data$nP) wneg[c,1:(data$nP-npos[c])] <- which(data$infM[c,]==0)
  
}
data$npos <- npos # N pos pathogens per infection status
data$wpos <- wpos # index pos pathogens per infection status
data$wneg <- wneg # index neg pathogens per infection status


#--- chain starting values
init <- function(data, nChains){
  
  ii <- init <- list()
  for(i in 1:nChains){
    init$sero <- array(runif(data$nPp*data$nL*data$nA, 0.05, 0.4), dim=c(data$nL, data$nA, data$nPp))
    init$sd0 <- runif(1, 0.4, 0.7)
    init$sd1 <- runif(1, 0.5, 1)
    init$mu0 <- runif(data$nP, 0, 0.2)
    init$mu1 <- runif(data$nPp, 1.5, 3.5)
    init$phi <- runif((data$nP*data$nPp-(data$nPp*2)), 0.01, 0.5)
    init$rho00 <- runif(1, 0.4, 0.7)
    init$sdE <- runif(data$nC, 0.3, 0.7)
    init$muE <- runif(data$nPp, 0.1, 4)
    ii[[i]] <- init
  } 
  
  return(ii)
}
ini <- init(data, nChains = 3)


#--- fit model
## NOTE: this can take several hours or days depending on the number of Gaussian components 
## being fit and their dimensions - we recommend using a high computing cluster where possible
check_cmdstan_toolchain()
set_cmdstan_path('C:/Users/Megan/Documents/.cmdstanr/cmdstan-2.26.1')
mod <- cmdstan_model(here('StanModels', 'MultiSero_LocAgeELISA.stan'), pedantic=F) 

# output path
folder <- paste(present,collapse='+')
dir.create(here('Results', folder))

# run model
fit <- mod$sample(data=data, chains=3, parallel_chains=3, iter_sampling=3000, refresh=100, 
                  iter_warmup=3000, init=ini, output_dir=here('Results', folder))



#--- check convergence & save chains
chains <- fit$draws(format='df')
fwrite(chains, 'Chains.csv')
color_scheme_set("mix-blue-red")
mcmc_trace(chains, regex_pars = c("seroAll","lp__"))
mcmc_trace(chains, regex_pars = c("mu0","mu1","muE"))
mcmc_trace(chains, regex_pars = c('sd0','sd1','sdE'))
mcmc_trace(chains, regex_pars = c('phi','rho00'))


#--- LogLik, WAIC, LOO, etc
ll <- fit$draws('log_lik')
loo_1 <- loo(fit$draws("log_lik"))
waicloo <- rbind(as.data.frame(loo_1$estimates), as.data.frame(waic(ll)$estimates))
lli <- data.frame(par='LogLik',med=NA,ciL=NA,ciU=NA)
lli[1,2:4] <- quantile(chains$sumloglik, c(0.5,0.025,0.975))


#--- extract pars & plots
loc <- unique(df$upazila)
chains <- as.data.frame(chains)


# prev
sero <- extract_sero(chains, data, pathogens)
seroLoc <- extract_seroLoc(chains, data, pathogens, loc)
seroAge <- extract_seroAge(chains, data, pathogens, ageG=c('0-9','10-19','20-29','30-39','40-49','50-59','60+'))
seroLocAge <- extract_seroLocAge(chains, data, pathogens, loc, ageG=c('0-9','10-19','20-29','30-39','40-49','50-59','60+'))

# cross-reactivity & rho
phi <- extract_phi(chains, data, pathogens)

# means
mu <- extract_mu(chains, data, pathogens=pathogens)

# sds
sds <- extract_sds(chains, data)
sdsE <- extract_sdsELISA(chains, data, pathogens)


#--- save results
write.csv(sero, here('Results', folder, 'Prev.csv'), row.names=F)
write.csv(seroLoc, here('Results', folder, 'PrevLoc.csv'), row.names=F)
write.csv(seroAge, here('Results', folder, 'PrevAge.csv'), row.names=F)
write.csv(seroLocAge, here('Results', folder, 'PrevLocAge.csv'), row.names=F)
write.csv(phi$phi, here('Results', folder, 'CrossReactivity.csv'), row.names=F)
write.csv(phi$rho, here('Results', folder, 'Rho.csv'), row.names=F)
write.csv(mu$mus0, here('Results', folder, 'MeanNeg.csv'), row.names=F)
write.csv(mu$mus1, here('Results', folder, 'MeanPos.csv'), row.names=F)
write.csv(sds, here('Results', folder, 'SDs.csv'), row.names=F)
write.csv(sdsE, here('Results', folder, 'SDs_ELISA.csv'), row.names=F)
write.csv(lli, here('Results', folder, 'LogLik.csv'), row.names=F)
write.csv(waicloo, here('Results', folder, 'WAICLOO.csv'))



#--- model fits (this can take a lot of time & memory if all chain iterations are used)
# thinning of the chains can speed this up if needed
distfits <- plot_distsLocAge(chains, data, pathogens=pathogens, NperLA=data$NperLA)
png(filename='DistFitsV1.png', width=20, height=20, res=400, units='cm')
plot(distfits$fit)
dev.off()
png(filename='DistFitsV2.png', width=20, height=20, res=400, units='cm')
plot(distfits$fitPN)
dev.off()
