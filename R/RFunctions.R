#----- Functions for multivariate Gaussian mixture serology model -----#
library(emdbook)
library(ggplot2)
library(cowplot)
library(RColorBrewer)
library(matrixStats)
library(stringr)


#----- Generate infection status matrix
inf_matrix <- function(N_pathogen, pres=rep(1,N_pathogen)){
  
  # list of possible outcomes for each pathogen
  combos <- list()
  for(c in 1:N_pathogen) combos[[c]] <- c(0,1)
  
  # matrix of all possible infection status combinations
  m <- expand.grid(combos)
  colnames(m) <- letters[1:N_pathogen]
  
  # remove positives of absent pathogens
  if(sum(pres)<N_pathogen){
    for(abs in which(pres==0)) m <- m[m[,abs]==0, ]
  }
  
  return(m)
}


#----- Extract prevalence estimates
extract_sero <- function(chains, data, pathogens){
  
  sero <- data.frame(pathogen=pathogens, med=NA, ciL=NA, ciU=NA)
  
  for(p in 1:data$nP) sero[p,2:4] <- quantile(chains[,paste('seroAll[', paste(p,']', sep=''), sep='')], c(0.5,0.025,0.975))
  sero <- sero[!sero$med==0, ]
  
  return(sero)
}


#----- Extract prevalence estimates by location
extract_seroLoc <- function(chains, data, pathogens, loc){
  
  sero <- data.frame(pathogen=NA, loc=NA, med=NA, ciL=NA, ciU=NA)
  ind <- 1
  for(p in 1:data$nPp) for(l in 1:data$nL){
    lp <- paste(l,p,sep=',')
    sero[ind,3:5] <- quantile(chains[,paste('seroLoc[', paste(lp,']', sep=''), sep='')], c(0.5,0.025,0.975))
    sero[ind,1:2] <- c(pathogens[p], loc[l])
    ind <- ind+1
  }
  
  return(sero)
}


#----- Extract prevalence estimates by age
extract_seroAge <- function(chains, data, pathogens, ageG){
  
  sero <- data.frame(pathogen=NA, age=NA, med=NA, ciL=NA, ciU=NA)
  ind <- 1
  for(p in 1:data$nPp) for(a in 1:data$nA){
    ap <- paste(a,p,sep=',')
    sero[ind,3:5] <- quantile(chains[,paste('seroAge[', paste(ap,']', sep=''), sep='')], c(0.5,0.025,0.975))
    sero[ind,1:2] <- c(pathogens[p], ageG[a])
    ind <- ind+1
  }
  
  return(sero)
}




#----- Extract prevalence estimates by location & age
extract_seroLocAge <- function(chains, data, pathogens, loc, ageG){
  
  sero <- data.frame(loc=NA, age=NA, pathogen=NA, med=NA, ciL=NA, ciU=NA)
  ind <- 1
  for(p in 1:data$nPp) for(l in 1:data$nL) for(a in 1:data$nA){
    
    k <- paste(l,a,sep=',')
    sero[ind,1:3] <- c(loc[l],ageG[a],pathogens[p])
    sero[ind,4:6] <- quantile(chains[,paste('sero[',paste(k,paste(p,']',sep=''),sep=','),sep='')], c(0.5,0.025,0.975))
    ind <- ind+1
  }
  
  return(sero)
}



#----- Extract cross-reactivity & correlation estimates
extract_phi <- function(chains, data, pathogens){
  
  phi <- data.frame(pos=NA, neg=NA, med=NA, ciL=NA, ciU=NA)
  ind <- 1
  for(p in 1:data$nPp) for(p2 in 1:data$nP){
    if(!p==p2){
      phi[ind,1:2] <- c(pathogens[p], pathogens[p2])
      y <- str_replace_all(toString(c(p,p2))," ","")
      phi[ind,3:5] <- quantile(chains[,paste('CR[', paste(y, ']', sep=''), sep='')], c(0.5,0.025,0.975))
      ind <- ind+1
    }
  }
  
  rho <- data.frame(pars=c('rho00'), med=NA, ciL=NA, ciU=NA)
  rho[1,2:4] <- quantile(chains[,paste('rho00')], c(0.5,0.025,0.975))
  
  return(list(phi=phi,rho=rho))
  
}


#----- Extract gaussian means
extract_mu <- function(chains, data, pathogens){
  
  
  # label combination positives
  pos <- rep('neg', data$nC)
  for(c in 1:data$nC){
    np <- sum(data$infM[c,])
    if(np==1) pos[c] <- pathogens[which(data$infM[c,]==1)]
    else if(np>1) pos[c] <- paste(pathogens[which(data$infM[c,]==1)], sep='&', collapse="&")
  }
  
  # all gaussian means
  mus0 <- data.frame(pg=rep(NA,length(which(data$infM==0))), pos=NA, med=NA, ciL=NA, ciU=NA)
  mus1 <- data.frame(pg=rep(NA,length(which(data$infM==1))), pos=NA, med=NA, ciL=NA, ciU=NA)
  ix0 <- 1
  ix1 <- 1
  for(c in 1:data$nC) for(p in 1:data$nP){
    if(data$infM[c,p]==0){
      mus0$pg[ix0] <- pathogens[p]
      mus0$pos[ix0] <- pos[c] 
      y <- paste(c,p, sep=',')
      mus0[ix0,3:5] <- quantile(chains[,paste(paste('mu[',y,sep=''),']',sep='')], c(0.5,0.025,0.975))
      ix0 <- ix0+1
    }else{
      mus1$pg[ix1] <- pathogens[p]
      mus1$pos[ix1] <- pos[c] 
      y <- paste(c,p, sep=',')
      mus1[ix1,3:5] <- quantile(chains[,paste(paste('mu[',y,sep=''),']',sep='')], c(0.5,0.025,0.975))
      ix1 <- ix1+1
    }
  }
  colnames(mus0)[1] <- 'antigen'
  colnames(mus1)[1] <- 'antigen'
  
  return(list(mus0=mus0, mus1=mus1))
}


#----- Extract gaussian sds
extract_sds <- function(chains, data){
  
  sig <- data.frame(par=c('sd0','sd1'), med=NA, ciL=NA, ciU=NA)
  sig[1,2:4] <- quantile(chains$sd0, c(0.5,0.025,0.975))
  sig[2,2:4] <- quantile(chains$sd1, c(0.5,0.025,0.975))
  
  return(sig)
}


#----- Extract gaussian sds for ELISA
extract_sdsELISA <- function(chains, data, pathogens){
  
  sig <- data.frame(par=paste('sdE',seq(1,data$nC),sep='_'), med=NA, ciL=NA, ciU=NA)
  sdE <- paste('sde[', paste(seq(1,data$nC),']',sep=''),sep='')
  for(i in 1:data$nC) sig[i,2:4] <- quantile(chains[,sdE[i]], c(0.5,0.025,0.975))
  
  # label 
  sig$status <- 'negative'
  for(c in 2:nrow(sig)) sig$status[c] <- paste(pathogens[data$wpos[c,data$wpos[c,]>0]], collapse='+')
  
  return(sig)
}



#----- Extract covariance matrices per iteration
extract_covM <- function(chains, data){
  
  iter <- length(chains$lp__)
  covM <- list()
  for(i in 1:iter){
    covM[[i]] <- list()
    for(c in 1:data$nC){
      x <- matrix(NA, ncol=data$nP, nrow=data$nP)
      for(p in 1:data$nP) for(p2 in 1:data$nP){
        
        y <- paste(paste(c,p,sep=','),p2,sep=',')
        x[p,p2] <- chains[i,paste(paste('covM[',y,sep=''),']',sep='')]
        
      }
      covM[[i]][[c]] <- x
    }
  }
  return(covM)
}


#----- Plot gaussian distribution fits
plot_dists <- function(chains, data, pathogens){
  
  iter <- length(chains$lp__)
  covM <- extract_covM(chains, data)
  
  # simulate multivariate gaussians per combination
  yy <- yp <- list()
  for(p in 1:data$nP) yp[[p]] <- matrix(NA, nrow=512, ncol=iter)
  ypN <- ypP <- yp
  for(i in 1:iter){
    for(c in 1:data$nC){
      
      # simulate gaussion for combination c, iteration i
      g <- paste('mu[',c, sep='')
      nn <- ceiling(sum(data$N*chains[,paste('theta[',paste(c,']',sep=''), sep='')][i]))
      muu <- vector()
      for(p in 1:data$nP) muu[p] <- chains[i,paste(paste(g,p,sep=','),']',sep='')]
      yy[[c]] <- as.data.frame(rmvnorm(nn, mean=muu,sigma=covM[[i]][[c]]))
      yy[[c]]$C <- c
    }
    
    # density distributions per pathogen
    yc <- do.call('rbind',yy)
    for(p in 1:data$nP){
      
      yp[[p]][,i] <- density(yc[,p], bw=0.01, from=-2.5, to=5.5)$y
      
      if(data$pres[p]==1){
        z <- which(data$infM[,p]==1)
        pw <- vector()
        for(s in 1:length(z)) pw[s] <- chains[i,paste(paste('theta[',paste(z[s]),sep=''),']',sep='')] 
        propP <- sum(data$N*pw)/data$N
        ypP[[p]][,i] <- density(yc[yc$C %in% which(data$infM[,p]==1),p], bw=0.01, from=-2.5, to=5.5)$y * propP
        ypN[[p]][,i] <- density(yc[yc$C %in% which(data$infM[,p]==0),p], bw=0.01, from=-2.5, to=5.5)$y *(1-propP)
      }else{
        ypN[[p]][,i] <- density(yc[yc$C %in% which(data$infM[,p]==0),p], bw=0.01, from=-2.5, to=5.5)$y
      }
    }
  }
  
  # quantiles of density distributions
  dpq <- dpqP <- dpqN <- list()
  titer <- density(yc[,1], bw=0.01, from=-2.5, to=5.5)$x
  for(p in 1:data$nP){
    dpq[[p]] <- as.data.frame(rowQuantiles(yp[[p]], probs=c(0.5,0.025,0.975)))
    dpqN[[p]] <- as.data.frame(rowQuantiles(ypN[[p]], probs=c(0.5,0.025,0.975)))
    dpqP[[p]] <- as.data.frame(rowQuantiles(ypP[[p]], probs=c(0.5,0.025,0.975)))
    dpq[[p]]$titer <- dpqN[[p]]$titer <- dpqP[[p]]$titer <- titer
    dpq[[p]]$pathogen <- dpqN[[p]]$pathogen <- dpqP[[p]]$pathogen <- pathogens[p]
  }
  dpq <- do.call('rbind',dpq)
  dpqN <- do.call('rbind',dpqN)
  dpqP <- do.call('rbind',dpqP)
  colnames(dpq)[1:3] <- colnames(dpqN)[1:3] <- colnames(dpqP)[1:3] <- c('med','ciL','ciU')
  
  # compile data for plotting 
  dta <- as.data.frame(data$y)
  colnames(dta) <- pathogens
  dta <- tidyr::gather(dta, key='pathogen', value='t')
  
  # overall fit
  fitD <- ggplot()+ geom_histogram(data=dta, aes(t,y=..density..), bins=150, fill='grey80', col='grey70')+
    theme_minimal()+ theme(text=element_text(size=18))+ #xlim(-1.5,3.5)+
    geom_line(data=dpq, aes(titer,med),col='springgreen4')+ facet_wrap(~pathogen, scales='free_y')+ xlab('titer')+
    geom_ribbon(data=dpq, aes(x=titer,y=med,ymin=ciL,ymax=ciU), fill='springgreen3', alpha=0.4)
  
  # pos-neg fit
  fitDPN <- ggplot()+ geom_histogram(data=dta, aes(t,y=..density..), bins=150, fill='grey80', col='grey70')+
    theme_minimal()+ theme(text=element_text(size=18))+ #xlim(-1.5,3.5)+
    geom_line(data=dpqN, aes(titer,med),col='mediumblue')+
    geom_line(data=dpqP, aes(titer,med),col='violetred')+ facet_wrap(~pathogen, scales='free_y')+ xlab('titer')+
    geom_ribbon(data=dpqN, aes(x=titer,y=med,ymin=ciL,ymax=ciU), fill='mediumblue', alpha=0.3)+
    geom_ribbon(data=dpqP, aes(x=titer,y=med,ymin=ciL,ymax=ciU), fill='violetred', alpha=0.5)
  
  # return plots
  return(list(fit=fitD,fitPN=fitDPN)) 
}



#----- Plot gaussian distribution fits for location model
plot_distsLoc <- function(chains, data, pathogens){
  
  iter <- length(chains$lp__)
  covM <- extract_covM(chains, data)
  
  # simulate multivariate gaussians per combination
  yy <- yp <- list()
  for(p in 1:data$nP) yp[[p]] <- matrix(NA, nrow=512, ncol=iter)
  ypN <- ypP <- yp
  for(i in 1:iter){
    for(c in 1:data$nC){
      
      # simulate gaussion for combination c, iteration i
      g <- paste('mu[',c, sep='')
      nn <- 0
      for(l in 1:data$nL){
        lc <- paste(l,c, sep=',')
        nn <- nn + round(data$NperL[l]*chains[,paste('theta[',paste(lc,']',sep=''), sep='')][i])
      }
      
      muu <- vector()
      for(p in 1:data$nP) muu[p] <- chains[i,paste(paste(g,p,sep=','),']',sep='')]
      if(nn==0) nn <- 1
      yy[[c]] <- as.data.frame(rmvnorm(nn, mean=muu,sigma=covM[[i]][[c]]))
      yy[[c]]$C <- c
    }
    
    # density distributions per pathogen
    yc <- do.call('rbind',yy)
    for(p in 1:data$nP){
      
      yp[[p]][,i] <- density(yc[,p], bw=0.01, from=-2.5, to=5.5)$y
      
      if(data$pres[p]==1){
        z <- which(data$infM[,p]==1)
        pw <- rep(0, data$nL)
        for(l in 1:data$nL) for(s in 1:length(z)){
          lz <- paste(l,z[s], sep=',')
          pw[l] <- pw[l] + chains[i,paste(paste('theta[',lz,sep=''),']',sep='')]
        }
        
        propP <- sum(data$NperL*pw)/data$N
        ypP[[p]][,i] <- density(yc[yc$C %in% which(data$infM[,p]==1),p], bw=0.01, from=-2.5, to=5.5)$y * propP
        ypN[[p]][,i] <- density(yc[yc$C %in% which(data$infM[,p]==0),p], bw=0.01, from=-2.5, to=5.5)$y *(1-propP)
      }else{
        ypN[[p]][,i] <- density(yc[yc$C %in% which(data$infM[,p]==0),p], bw=0.01, from=-2.5, to=5.5)$y
      }
    }
  }
  
  # quantiles of density distributions
  dpq <- dpqP <- dpqN <- list()
  titer <- density(yc[,1], bw=0.01, from=-2.5, to=5.5)$x
  for(p in 1:data$nP){
    dpq[[p]] <- as.data.frame(rowQuantiles(yp[[p]], probs=c(0.5,0.025,0.975)))
    dpqN[[p]] <- as.data.frame(rowQuantiles(ypN[[p]], probs=c(0.5,0.025,0.975)))
    dpqP[[p]] <- as.data.frame(rowQuantiles(ypP[[p]], probs=c(0.5,0.025,0.975)))
    dpq[[p]]$titer <- dpqN[[p]]$titer <- dpqP[[p]]$titer <- titer
    dpq[[p]]$pathogen <- dpqN[[p]]$pathogen <- dpqP[[p]]$pathogen <- pathogens[p]
  }
  dpq <- do.call('rbind',dpq)
  dpqN <- do.call('rbind',dpqN)
  dpqP <- do.call('rbind',dpqP)
  colnames(dpq)[1:3] <- colnames(dpqN)[1:3] <- colnames(dpqP)[1:3] <- c('med','ciL','ciU')
  
  # compile data for plotting 
  dta <- as.data.frame(data$y)
  colnames(dta) <- pathogens
  dta <- tidyr::gather(dta, key='pathogen', value='t')
  
  # overall fit
  fitD <- ggplot()+ geom_histogram(data=dta, aes(t,y=..density..), bins=150, fill='grey80', col='grey70')+
    theme_minimal()+ theme(text=element_text(size=18))+ #xlim(-1.5,3.5)+
    geom_line(data=dpq, aes(titer,med),col='springgreen4')+ facet_wrap(~pathogen, scales='free_y')+ xlab('titer')+
    geom_ribbon(data=dpq, aes(x=titer,y=med,ymin=ciL,ymax=ciU), fill='springgreen3', alpha=0.4)
  
  # pos-neg fit
  fitDPN <- ggplot()+ geom_histogram(data=dta, aes(t,y=..density..), bins=150, fill='grey80', col='grey70')+
    theme_minimal()+ theme(text=element_text(size=18))+ #xlim(-1.5,3.5)+
    geom_line(data=dpqN, aes(titer,med),col='mediumblue')+
    geom_line(data=dpqP, aes(titer,med),col='violetred')+ facet_wrap(~pathogen, scales='free_y')+ xlab('titer')+
    geom_ribbon(data=dpqN, aes(x=titer,y=med,ymin=ciL,ymax=ciU), fill='mediumblue', alpha=0.3)+
    geom_ribbon(data=dpqP, aes(x=titer,y=med,ymin=ciL,ymax=ciU), fill='violetred', alpha=0.5)
  
  # return plots
  return(list(fit=fitD,fitPN=fitDPN)) 
}



#----- Plot gaussian distribution fits for location & age model 
plot_distsLocAge <- function(chains, data, pathogens, NperLA){
  
  iter <- length(chains$lp__)
  covM <- extract_covM(chains, data)
  
  # simulate multivariate gaussians per combination
  yy <- yp <- list()
  for(p in 1:data$nP) yp[[p]] <- matrix(NA, nrow=512, ncol=iter)
  ypN <- ypP <- yp
  for(i in 1:iter){
    for(c in 1:data$nC){
      
      # simulate gaussion for combination c, iteration i
      g <- paste('mu[',c, sep='')
      nn <- 0
      for(l in 1:data$nL) for(a in 1:data$nA){
        la <- paste(l,a, sep=',')
        nn <- nn + round(NperLA[a,l]*chains[,paste('theta[',paste(paste(la,c,sep=','),']',sep=''), sep='')][i])
      }
      
      muu <- vector()
      for(p in 1:data$nP) muu[p] <- chains[i,paste(paste(g,p,sep=','),']',sep='')]
      if(nn==0) nn <- 1
      yy[[c]] <- as.data.frame(rmvnorm(nn, mean=muu,sigma=covM[[i]][[c]]))
      yy[[c]]$C <- c
    }
    
    # density distributions per pathogen
    yc <- do.call('rbind',yy)
    for(p in 1:data$nP){
      
      yp[[p]][,i] <- density(yc[,p], bw=0.01, from=-2.5, to=5.5)$y
      
      if(data$pres[p]==1){
        z <- which(data$infM[,p]==1)
        pw <- matrix(0, nrow=data$nA, ncol=data$nL)
        for(l in 1:data$nL) for(a in 1:data$nA) for(s in 1:length(z)){
          la <- paste(l,a, sep=',')
          pw[a,l] <- pw[a,l] + chains[i,paste(paste('theta[',paste(la,z[s],sep=','),sep=''),']',sep='')]
        }
        
        propP <- sum(NperLA*pw)/data$N
        ypP[[p]][,i] <- density(yc[yc$C %in% which(data$infM[,p]==1),p], bw=0.01, from=-2.5, to=5.5)$y * propP
        ypN[[p]][,i] <- density(yc[yc$C %in% which(data$infM[,p]==0),p], bw=0.01, from=-2.5, to=5.5)$y *(1-propP)
      }else{
        ypN[[p]][,i] <- density(yc[yc$C %in% which(data$infM[,p]==0),p], bw=0.01, from=-2.5, to=5.5)$y
      }
    }
  }
  
  # quantiles of density distributions
  dpq <- dpqP <- dpqN <- list()
  titer <- density(yc[,1], bw=0.01, from=-2.5, to=5.5)$x
  for(p in 1:data$nP){
    dpq[[p]] <- as.data.frame(rowQuantiles(yp[[p]], probs=c(0.5,0.025,0.975)))
    dpqN[[p]] <- as.data.frame(rowQuantiles(ypN[[p]], probs=c(0.5,0.025,0.975)))
    dpqP[[p]] <- as.data.frame(rowQuantiles(ypP[[p]], probs=c(0.5,0.025,0.975)))
    dpq[[p]]$titer <- dpqN[[p]]$titer <- dpqP[[p]]$titer <- titer
    dpq[[p]]$pathogen <- dpqN[[p]]$pathogen <- dpqP[[p]]$pathogen <- pathogens[p]
  }
  dpq <- do.call('rbind',dpq)
  dpqN <- do.call('rbind',dpqN)
  dpqP <- do.call('rbind',dpqP)
  colnames(dpq)[1:3] <- colnames(dpqN)[1:3] <- colnames(dpqP)[1:3] <- c('med','ciL','ciU')
  
  # compile data for plotting 
  dta <- as.data.frame(data$y)
  colnames(dta) <- pathogens
  dta <- tidyr::gather(dta, key='pathogen', value='t')
  
  # overall fit
  fitD <- ggplot()+ geom_histogram(data=dta, aes(t,y=..density..), bins=150, fill='grey80', col='grey70')+
    theme_minimal()+ theme(text=element_text(size=18))+ #xlim(-1.5,3.5)+
    geom_line(data=dpq, aes(titer,med),col='springgreen4')+ facet_wrap(~pathogen, scales='free_y')+ xlab('titer')+
    geom_ribbon(data=dpq, aes(x=titer,y=med,ymin=ciL,ymax=ciU), fill='springgreen3', alpha=0.4)
  
  # pos-neg fit
  fitDPN <- ggplot()+ geom_histogram(data=dta, aes(t,y=..density..), bins=150, fill='grey80', col='grey70')+
    theme_minimal()+ theme(text=element_text(size=18))+ #xlim(-1.5,3.5)+
    geom_line(data=dpqN, aes(titer,med),col='mediumblue')+
    geom_line(data=dpqP, aes(titer,med),col='violetred')+ facet_wrap(~pathogen, scales='free_y')+ xlab('titer')+
    geom_ribbon(data=dpqN, aes(x=titer,y=med,ymin=ciL,ymax=ciU), fill='mediumblue', alpha=0.3)+
    geom_ribbon(data=dpqP, aes(x=titer,y=med,ymin=ciL,ymax=ciU), fill='violetred', alpha=0.5)
  
  # return plots
  return(list(fit=fitD,fitPN=fitDPN)) 
}

