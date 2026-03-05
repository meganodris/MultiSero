data {
  int N; // N individuals
  int nP; // N antigens
  int nPp; // N present pathogens
  array[nP] int pres; // indicator for present pathogens
  int nC; // N infection status combinations
  array[N] vector[nP] y; // antibody titer data
  matrix[nC, nP] infM; // infection status combination indicator
  array[nC] int npos; // N pos indicator
  array[nC, nP] int wpos; // which pos pathogens indicator
  array[nC, nP] int wneg; // which neg pathogens indicator
  int nL; // N locations
  array[N] int loc; // location
  vector[nL] NperL; // N per location
  int nA; // N age groups
  array[N] int ageG; // age group
  array[nL] simplex[nA] ageProp; // proportion of study pop by age per location
  array[nA] vector[nL] NperLA; // N individuals per location & age group
  
  // prior hyperparameters
  real mu0_mean;
  real<lower=0> mu0_sd;
  real mu1_mean;
  real<lower=0> mu1_sd;
  real<lower=0> sd0_mean;
  real<lower=0> sd0_sd;
  real<lower=0> sd1_mean;
  real<lower=0> sd1_sd;
  real<lower=0> phi_rate;
  real<lower=0> rho_a;
  real<lower=0> rho_b;
}

transformed data {

  matrix[nP, N] Y;
  for (n in 1:N) Y[ : , n] = y[n];
  real log_2pi = log(2 * pi());

  // Pre-calculate 1-infM for log_theta calculation
  matrix[nC, nP] infM_inv = 1.0 - infM;

  // Identify which indices in infM correspond to present pathogens
  array[nPp] int pres_idx;
  {
    int k = 1;
    for (p in 1 : nP) {
      if (pres[p] == 1) {
        pres_idx[k] = p;
        k += 1;
      }
    }
  }
}

parameters {
  array[nL, nA] vector<lower=0, upper=1>[nPp] sero; // infection prevalence
  real<lower=0> sd0; // sd neg
  real<lower=0> sd1; // sd pos
  vector<lower=-4>[nP] mu0; // mean neg
  vector<lower=0>[nPp] mu1; // mean pos
  vector<lower=0>[(nP * nPp) - (nPp)] phi; // relative cross-reactive titer increase
  real<lower=0, upper=1> rho00; // correlation in neg titers
}

transformed parameters {

  array[nC] vector[nP] mu; // gaussian means
  array[nC] vector[nP] sigma; // gaussian sds
  matrix[N, nC] pC; // probabilities per individual & component
  vector[N] log_lik; // individual likelihoods
  array[nL, nA] simplex[nC] theta; // gaussian weights
  matrix[nP, nP] CR = rep_matrix(0, nP, nP); // relative cross-reactive titer increases
  array[nC] cholesky_factor_cov[nP] L; // cholesky factors
  array[nC] matrix[nP, nP] covM;

  {
    // Local variables
    array[nL, nA] vector[nC] log_theta;
    int ix = 1;
    real cv;
    real var_sum;

    //--- 1. Vectorized log_theta calculation ---//
    for (l in 1 : nL) {
      for (a in 1 : nA) {
        vector[nP] log_sero_full = rep_vector(0, nP);
        vector[nP] log1m_sero_full = rep_vector(0, nP);
        for (k in 1 : nPp) {
          log_sero_full[pres_idx[k]] = log(sero[l, a, k]);
          log1m_sero_full[pres_idx[k]] = log1m(sero[l, a, k]);
        }
        log_theta[l, a] = infM * log_sero_full + infM_inv * log1m_sero_full;
        theta[l, a] = exp(log_theta[l, a]);
      }
    }

    //--- 2. Cross reactivity ---//
    for (p in 1 : nPp)
      for (p2 in 1 : (nP)) {
        if (p == p2)
          CR[p, p2] = 0;
        else {
          CR[p, p2] = phi[ix];
          ix = ix + 1;
        }
      }

    //--- 3. Gaussian means & sds ---//
    sigma[1] = rep_vector(sd0, nP);
    mu[1,] = mu0;

    for (c in 2 : nC) {
      for (p in 1 : npos[c]) {
        sigma[c, wpos[c, p]] = sd1;
        mu[c, wpos[c, p]] = mu0[wpos[c, p]] + mu1[wpos[c, p]];
      }

      if (npos[c] == 1) {
        for (p in 1 : (nP - npos[c])) {
          sigma[c, wneg[c, p]] = sqrt(
                                      sd0 ^ 2
                                      + (CR[wpos[c, 1], wneg[c, p]] * sd1)
                                        ^ 2);
          mu[c, wneg[c, p]] = mu0[wneg[c, p]]
                              + CR[wpos[c, 1], wneg[c, p]] * mu1[wpos[c, 1]];
        }
      } else {
        for (p in 1 : (nP - npos[c])) {
          mu[c, wneg[c, p]] = mu0[wneg[c, p]];
          var_sum = sd0 ^ 2;
          for (j in 1 : npos[c]) {
            mu[c, wneg[c, p]] += CR[wpos[c, j], wneg[c, p]] * mu1[wpos[c, j]];
            var_sum += (CR[wpos[c, j], wneg[c, p]] * sd1) ^ 2;
          }
          sigma[c, wneg[c, p]] = sqrt(var_sum);
        }
      }
    }

    //--- 4. Covariance matrices & Cholesky ---//
    for (c in 1 : nC) {
      for (p in 1 : nP)
        covM[c, p, p] = sigma[c, p] ^ 2;
      for (p in 1 : (nP - 1)) {
        for (p2 in (p + 1) : nP) {
          if (infM[c, p] + infM[c, p2] == 2) {
            cv = 0;
          } else if (infM[c, p] + infM[c, p2] == 0) {
            cv = rho00 * sd0 ^ 2;
            for (j in 1 : npos[c])
              cv += CR[wpos[c, j], p] * CR[wpos[c, j], p2] * sd1 ^ 2;
          } else {
            cv = (infM[c, p] == 1) ? CR[p, p2] * sd1 ^ 2
                 : CR[p2, p] * sd1 ^ 2;
          }
          covM[c, p, p2] = cv;
          covM[c, p2, p] = cv;
        }
      }
      L[c] = cholesky_decompose(covM[c]);
    }

    //--- 5. Optimized Likelihood calculation ---//
    for(c in 1:nC) {

      real constant_part = -0.5 * (nP * log_2pi + 2 * sum(log(diagonal(L[c]))));
      matrix[nP, N] z = mdivide_left_tri_low(L[c], Y - rep_matrix(mu[c], N));

      for(n in 1:N) pC[n, c] = log_theta[loc[n], ageG[n], c] + constant_part - 0.5 * dot_self(z[ : , n]);
    }

    for(n in 1:N) log_lik[n] = log_sum_exp(pC[n, :]);
  }
}

model {

  // priors
  for (p in 1:nPp) for (l in 1 : nL) sero[l,  : , p] ~ beta(1, 5);
  mu0 ~ normal(mu0_mean, mu0_sd); //0, 0.1
  mu1 ~ normal(mu1_mean, mu1_sd); //3, 0.1
  sd0 ~ normal(sd0_mean, sd0_sd); //0.4,0.05
  sd1 ~ normal(sd1_mean, sd1_sd); //0.7,0.05
  phi ~ exponential(phi_rate); //8
  rho00 ~ beta(rho_a, rho_b); //2,2

  target += sum(log_lik);

}

generated quantities {

  real sumloglik = sum(log_lik);
  vector[nP] seroAll = rep_vector(0, nP);
  array[nL] vector[nP] seroLoc = rep_array(rep_vector(0, nP), nL);
  array[nA] vector[nP] seroAge = rep_array(rep_vector(0, nP), nA);
  array[nL, nA] vector[nP] seroLocAge = rep_array(rep_vector(0, nP), nL, nA);

  for (p in 1 : nP) {
    if (pres[p] == 1) {
      int pp_idx = sum(pres[1 : p]);
      for (l in 1 : nL) {
        for (a in 1 : nA)
          seroLocAge[l, a, p] = sero[l, a, pp_idx];
        seroLoc[l, p] = sum(ageProp[l] .* to_vector(sero[l,  : , pp_idx]));
      }
      for (a in 1 : nA) {
        vector[nL] temp_sero;
        for (l in 1 : nL)
          temp_sero[l] = sero[l, a, pp_idx];
        seroAge[a, p] = sum(NperLA[a] .* temp_sero) / sum(NperLA[a]);
      }
      seroAll[p] = sum(to_vector(seroLoc[ : , p]) .* NperL) / N;
    }
  }
}
