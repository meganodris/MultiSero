#----- Simulate multi-pathogen cross-sectional serology data -----#
library(lhs)
library(matrixStats)
library(mvtnorm)
library(matrixcalc)
library(here)

#----- Function to generate infection status combination matrix
inf_matrix <- function(N_pathogen, pres){
  
  # list of possible outcomes for each pathogen
  combos <- list()
  for(c in 1:N_pathogen) combos[[c]] <- c(0,1)
  
  # matrix of all possible infection status combinations
  m <- expand.grid(combos)
  colnames(m) <- letters[1:N_pathogen]
  
  # remove positives of any absent pathogens
  if(sum(pres)<N_pathogen){
    for(abs in which(pres==0)) m <- m[m[,abs]==0, ]
  }
  
  return(m)
}


#----- Function to simulate multivariate Gaussian data
sim_multisero <- function(N, nP, nPp, sero, mu0, mu1, sd0, sd1, CR, corr00){
  
  #--- infection status matrix
  pres <- c(rep(1,nPp), rep(0,nP-nPp)) # index for which pathogens are present (1) or absent (0)
  infM <- inf_matrix(nP, pres=pres) # matrix of possible infection statuses
  
  #--- Gaussian weights
  w <- infM
  w[,which(pres==0)] <- 1
  for(p in 1:nPp) for(c in 1:nrow(infM)){
    if(infM[c,p]==0) w[c,p] <- 1 - sero[p]
    else w[c,p] <- sero[p]
  }
  ws <- rowProds(as.matrix(w)) # weights
  Nws <- round(N*ws) # N individuals per infection status
  
  #--- Gaussian means & standard deviations
  mus <- infM
  sigs <- infM
  for(c in 1:nrow(infM)){
    
    if(sum(infM[c,])==0){ # neg to all
      
      mus[c,] <- mu0 
      sigs[c,] <- sd0
      
    }else if(sum(infM[c,])==1){ # pos to just 1
      
      # which are pos/neg
      wp <- which(infM[c,]==1)
      wn <- which(infM[c,]==0)
      
      # pos
      mus[c,wp] <- mu0[wp] + mu1[wp]
      sigs[c,wp] <- sd1[wp]
      
      # neg with cr
      for(n in wn) mus[c,n] <- mu0[n]+ mu1[wp]*CR[wp,n]
      for(n in wn) sigs[c,n] <- sqrt(sd0[n]^2 + (sd1[wp]*CR[wp,n])^2)
      
    }else if(sum(infM[c,])>1){ # pos to > 1
      
      # which are pos/neg
      wp <- which(infM[c,]==1)
      npos <- sum(infM[c,])
      wn <- which(infM[c,]==0)
      
      # pos
      for(p in wp) mus[c,p] <- mu0[p] + mu1[p]
      for(p in wp) sigs[c,p] <- sd1[p]
      
      # neg with cr
      for(n in wn){
        mus[c,n] <- mu0[n]
        for(p in wp) mus[c,n] <- mus[c,n] + mu1[p]*CR[p,n]
        tmp <- rep(0,ncol(infM))
        tmp2 <- rep(0,((npos*(npos-1))/2))
        for(p in wp) tmp[p] <- (sd1[p]*CR[p,n])^2
        sigs[c,n] <- sqrt(sd0[n]^2 + sum(tmp))
      }
    }
  }
  
  #--- covariance matrices
  covs <- list()
  for(c in 1:nrow(infM)) covs[[c]] <- matrix(0, ncol=nP,nrow=nP)
  for(c in 1:nrow(infM)){
    
    for(p in 1:nP) for(p2 in p:nP){
      
      if(p==p2){ # Gaussian variance
        
        covs[[c]][p,p2] <- sigs[c,p]^2 
        
      }else if(sum(infM[c,])==0){ # negative to all pathogens
        
        covs[[c]][p,p2] <- covs[[c]][p2,p] <- corr00*sigs[c,p]*sigs[c,p2]
        
      }else if(infM[c,p]+infM[c,p2]==0){ # negative to both where others are pos
        
        npos <- sum(infM[c,])
        wp <- which(infM[c,]==1)
        
        if(npos==1){
          cv00 <- corr00*sd0[p]*sd0[p2] + CR[wp,p]*CR[wp,p2]*sd1[wp]^2
        }else{
          vrs <- rep(0,npos)
          for(x in 1:npos) vrs[x] = CR[wp[x],p]*CR[wp[x],p2]*sd1[wp[x]]^2
          cv00 <- corr00*sd0[p]*sd0[p2] + sum(vrs)
        }
        covs[[c]][p,p2] <- covs[[c]][p2,p] <- cv00
        
        
      }else if(infM[c,p]+infM[c,p2]==1){ # negative & positive
        
        pos <- c(p,p2)[which(infM[c,c(p,p2)]==1)]
        neg <- c(p,p2)[which(infM[c,c(p,p2)]==0)]
        
        covs[[c]][p,p2] <- covs[[c]][p2,p] <- CR[pos,neg]*sd1[pos]^2
        
      }
    }
  }
  
  # check the covariance matrixes are semi definite
  Qcovs <- list()
  for(c in 1:nrow(infM)) Qcovs[[c]] <- is.positive.semi.definite(covs[[c]])
  
  #--- simulate multivariate Gaussians
  gs <- list()
  for(c in 1:nrow(infM)){
    if(Nws[c]>0){
      gs[[c]] <- as.data.frame(rmvnorm(Nws[c], as.numeric(mus[c,]), covs[[c]]))
      gs[[c]]$status <- c
      gs[[c]][,LETTERS[seq(1,nP)]] <- NA
      for(p in 1:nP) gs[[c]][,nP+1+p] <- infM[c,p]
      gs[[c]]$stat <- paste('{',paste(infM[c,],collapse=','),'}',sep='')
      
    }
  }
  gs <- do.call('rbind',gs)
  
  
  #--- return data
  return(list(titers=gs, mus=mus, sigs=sigs, covs=covs, covs_check=Qcovs))
}



#----- Set parameters for simulating
set.seed(1)
nP <- 5 # N antigens 
nPp <- 3 # N antigens that are present in the population


npars <- nP + nPp*2 + ((nP*nPp)-nPp) + 3 # calculate N parameters
params <- c(rep('sero',nPp),rep('mu0',nP),rep('mu1',nPp),rep('sd0',1),rep('sd1',1),
            rep('phi',((nP*nPp)-nPp)), 'corr00') # param names


nsims <- 5 # N datasets to simulate
pars <- randomLHS(nsims,npars) # randomly draw param values for each simulation
pars <- round(pars, 2) # round the param values to 2dp
colnames(pars) <- params


# transform param values to desired scales (here considered on a log RFI scale)
parsT <- pars
parsT[,colnames(parsT)=='mu0'] <- qunif(pars[,colnames(parsT)=='mu0'], min=-0.5, max=0.5) # mu0 
parsT[,colnames(parsT)=='mu1'] <- qunif(pars[,colnames(parsT)=='mu1'], min=1, max=3.5) # mu1
parsT[,colnames(parsT)=='sd0'] <- qunif(pars[,colnames(parsT)=='sd0'], min=0.25, max=0.5) # sd0
parsT[,colnames(parsT)=='sd1'] <- qunif(pars[,colnames(parsT)=='sd1'], min=0.5, max=0.75) # sd1


#----- Simulate data
sims <- list()
truepars <- truepars2 <- list()
N <- 1500 # N samples per simulation
for(i in 1:nrow(parsT)){
  
  # cross-reactivity matrix
  CR <- matrix(0, ncol=nP, nrow=nPp)
  ind <- 1
  for(p in 1:nPp) for(p2 in 1:nP){
    if(p==p2) NULL
    else{
      CR[p,p2] <- parsT[i,which(colnames(parsT)=='phi')][ind]
      ind <- ind+1
    }
  }
  
  # simulate
  sero <- c(parsT[i,which(colnames(parsT)=='sero')], rep(0, nP-nPp))
  mu0 <- parsT[i,which(colnames(parsT)=='mu0')]
  mu1 <- c(parsT[i,which(colnames(parsT)=='mu1')], rep(0, nP-nPp))
  sd0 <- parsT[i,which(colnames(parsT)=='sd0')]
  sd1 <- parsT[i,which(colnames(parsT)=='sd1')]
  sd0 <- rep(sd0, nP)
  sd1 <- rep(sd1, nP)
  corr00 <- parsT[i,which(colnames(parsT)=='corr00')]
  sims[[i]] <- sim_multisero(N=N, nP=nP, nPp=nPp, sero=sero, mu0=mu0, mu1=mu1,
                             sd0=sd0, sd1=sd1, CR=CR, corr00=corr00)
  
  
  # true param values
  truepars[[i]] <- data.frame(pars=c(rep('sero',nP),rep('mu0',nP),rep('mu1',nP),rep('sd0',nP),rep('sd1',nP)),
                              pathogen=LETTERS[seq(1,nP)], true=c(sero,mu0,mu1,sd0,sd1))
  truepars2[[i]] <- list(CR=CR, r00=corr00)
  print(i)
  
}


#- save simulated data
saveRDS(sims, here('data', 'SimulatedData.RDS'))
saveRDS(list(truepars, truepars2), here('data', 'SimTruePars.RDS'))
