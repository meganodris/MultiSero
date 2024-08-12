data {
 
 int N; // N individuals
 int nP; // N antigens
 int nPp; // N present pathogens
 int pres[nP]; // indicator for present pathogens
 int nC; // N infection status combinations
 vector[nP] y[N]; // antibody titer data
 matrix[nC,nP] infM; // infection status combination indicator
 int npos[nC]; // N pos indicator
 int wpos[nC,nP]; // which pos pathogens indicator
 int wneg[nC,nP]; // which neg pathogens indicator
 int indELISA; // index of ELISA assay

}


parameters {
 
 vector<lower=0,upper=1>[nPp] sero; // infection prevalence
 real<lower=0> sd0; // sd neg
 real<lower=0> sd1; // sd pos
 vector<lower=-1>[nP] mu0; // mean neg
 vector<lower=0>[nPp] mu1; // mean pos
 vector<lower=0>[(nP*nPp)-(nPp*2)] phi; // relative cross-reactive titer increase
 real <lower=0,upper=1> rho00; // correlation in neg titers
 real <lower=-1> muE[nPp]; // mean titers on ELISA 
 real<lower=0> sdE[nPp+1]; // sd titers on ELISA

}


transformed parameters {
  
  vector[nP] mu[nC]; // gaussian means
  vector[nP] sigma[nC]; // gaussian sds
  vector[nC] pC[N]; // probabilities per individual & component
  vector[N] log_lik; // individual likelihoods
  simplex[nC] theta; // gaussian weights
  vector[nP] seroAll = rep_vector(0,nP); // prevalence
  matrix[nP,nP] CR = rep_matrix(0,nP,nP); // relative cross-reactive titer increases
  cov_matrix[nP] covM[nC]; // covariance matrices
  real mue[nC]; // mean titers on ELISA 
  real sde[nC]; // sd titers on ELISA
  {
    
    // some temp variables
  vector[nP] W[nC];
  vector[nC] log_theta;
  int ix = 1;
  real cv;
  real sig;


  //--- gaussian weights ---//
  seroAll[1:nPp] = sero;
  for(c in 1:nC) for(p in 1:nP){
    if(infM[c,p]==0) W[c,p] = 1 - seroAll[p];
    else W[c,p] = seroAll[p];
  }
  for(c in 1:nC) theta[c] = prod(W[c,]);
  log_theta = log(theta);


  //--- cross reactivity ---//
  for(p in 1:nPp) for(p2 in 1:(nP-1)){
    if(p==p2) CR[p,p2] = 0;
    else{
      CR[p,p2] = phi[ix];
      ix = ix+1;
    }
  } 


  //--- gaussian means & sds ---//
  sigma[1,1:(nP-1)] = rep_vector(sd0, (nP-1));
  sigma[1,nP] = sdE[1]; // sd of neg elisa 
  mu[1,] = mu0;
  for(c in 2:nC){
    
    for(p in 1:npos[c]){ // positives
      sigma[c,wpos[c,p]] = sd1;
      mu[c,wpos[c,p]] = mu0[wpos[c,p]] + mu1[wpos[c,p]];
    }
    
    if(npos[c]==1){ // negatives with 1 positive
      for(p in 1:(nP-1)){
        sigma[c,wneg[c,p]] = sqrt(sigma[1,wneg[c,p]]^2 + (CR[wpos[c,1],wneg[c,p]]*sd1)^2);
        mu[c,wneg[c,p]] = mu0[wneg[c,p]] + CR[wpos[c,1],wneg[c,p]]*mu1[wpos[c,1]];
      }
    }else{ // negatives with >1 positive
      for(p in 1:(nP-npos[c])){
        
        // means
        mu[c,wneg[c,p]] = mu0[wneg[c,p]];
        for(j in 1:npos[c]) mu[c,wneg[c,p]] = mu[c,wneg[c,p]] + CR[wpos[c,j],wneg[c,p]]*mu1[wpos[c,j]];
        
        // sds
        sigma[c,wneg[c,p]] = sigma[1,wneg[c,p]]^2;
        for(j in 1:npos[c]) sigma[c,wneg[c,p]] = sigma[c,wneg[c,p]] + (CR[wpos[c,j],wneg[c,p]]*sd1)^2;
        sigma[c,wneg[c,p]] = sqrt(sigma[c,wneg[c,p]]);
      }
    }
  }
  
  // ELISA means and sds
  mue[1] = mu0[indELISA];
  for(c in 2:nC){
    if(npos[c]==1) mue[c] = muE[wpos[c,1]];
    else{
      mue[c] = max(muE[wpos[c,1:npos[c]]]);
    }
  }
  mu[1:nC,indELISA] = mue;
  sde[1] = sdE[1];
  for(c in 2:nC){
    sig = 0;
    for(p in 1:npos[c]) sig = sig + sdE[wpos[c,p]+1]^2;
    sde[c] = sqrt(sig);
  }
  sigma[,indELISA] = sde;   
  
  //--- covariance matrices ---//
  for(p in 1:nP) covM[1,p,p] = sigma[1,p]^2; // variances neg to all
  for(p in 1:(nP-1)) for(p2 in (p+1):nP){
    covM[1,p,p2] = rho00*sigma[1,p]*sigma[1,p2];
    covM[1,p2,p] = rho00*sigma[1,p]*sigma[1,p2];
  }
  for(p in 1:(nP-1)){
    covM[1,p,indELISA] = 0; // elisa corr00
    covM[1,indELISA,p] = 0; 
  } 
  for(c in 2:nC){
    
    for(p in 1:nP) covM[c,p,p] = sigma[c,p]^2; // variances
    
    for(p in 1:(nP-1)) for(p2 in (p+1):nP){ // covariances
      
      if(infM[c,p]+infM[c,p2]==2){ // pos to both
      
        covM[c,p,p2] = 0;
        covM[c,p2,p] = 0;
        
      }else if(infM[c,p]+infM[c,p2]==0){ // neg to both
        
        if(npos[c]==1){
          if(p2==indELISA){ // elisa
      
            covM[c,p,p2] = CR[wpos[c,1],p]*CR[wpos[c,1],p2]*sd1^2;
            covM[c,p2,p] = CR[wpos[c,1],p]*CR[wpos[c,1],p2]*sd1^2;
            
          }else{
            covM[c,p,p2] = rho00*sd0*sd0 + CR[wpos[c,1],p]*CR[wpos[c,1],p2]*sd1^2;
            covM[c,p2,p] = rho00*sd0*sd0 + CR[wpos[c,1],p]*CR[wpos[c,1],p2]*sd1^2;
          }
          
        }else{ // > 1 pos
          
          if(p2==indELISA){
            cv = 0;
          }else{
            cv = rho00*sd0*sd0;
          }
          for(j in 1:npos[c]) cv = cv + CR[wpos[c,j],p]*CR[wpos[c,j],p2]*sd1^2;
          covM[c,p,p2] = cv;
          covM[c,p2,p] = cv;
        }
        
      }else{ // neg and pos
        
        if(infM[c,p]==1){ // pos to p
          cv = CR[p,p2]*sd1^2;
        }else{ // pos to p2
          cv = CR[p2,p]*sd1^2;
        }
        covM[c,p,p2] = cv;
        covM[c,p2,p] = cv;
        
      }
    }
  }
    

  
  //--- likelihood calculation ---//
  for(c in 1:nC) for(n in 1:N) pC[n,c] = log_theta[c] + multi_normal_lpdf(y[n] | mu[c], covM[c]);
  for(n in 1:N) log_lik[n] = log_sum_exp(pC[n,]);
  
  }
}


model {
 
 // priors
 sero ~ beta(1,5);
 mu0 ~ normal(0,0.1);
 mu1 ~ normal(3,0.1);
 sd0 ~ normal(0.4,0.05);
 sd1 ~ normal(0.7,0.05);
 phi ~ exponential(8);
 rho00 ~ beta(2,2);
 muE ~ normal(2,1);
 sdE ~ normal(0.5,0.2);

 
 // log-likelihood
 target += sum(log_lik);
 

}


generated quantities {
  
  real sumloglik = sum(log_lik);

}

