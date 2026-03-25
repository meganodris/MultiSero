//---- MultiSero method for analysis of multi-pathogen serological studies
//---- cite: O'Driscoll et al., Science Translational Medicine, 2025


data {

  int N;                         // N individuals
  int nP;                        // N pathogens
  int nPp;                       // N present pathogens
  array[nP] int pres;            // indicator for present pathogens
  int nC;                        // N infection status combinations (2^nPp)
  array[N] vector[nP] y;         // antibody titer data
  matrix[nC, nP] infM;           // infection status combination indicator
  array[nC] int npos;            // N positive pathogens per combination
  array[nC, nP] int wpos;        // which pathogens are positive per combination
  array[nC, nP] int wneg;        // which pathogens are negative per combination

  // prevalence stratification
  int nS;                        // N strata (total)
  array[N] int strat;            // stratum index per individual (1-indexed, range 1..nS)
  vector[nS] NperS;              // N individuals per stratum (for overall prevalence weighting)

  // Aggregation: weighted combinations of strata
  // Set nAgg = 0 if no aggregated outputs are needed
  int nAgg;                      // N aggregated prevalence outputs
  matrix[nAgg, nS] aggW;         // aggregation weight matrix (rows sum to 1)

  // Prior parameters
  array[2] vector[nP] prior_mu0; // normal prior: [1] = mean, [2] = sd
  array[2] vector[nPp] prior_mu1;
  array[2] real sigma_scale;

}


transformed data {

  matrix[nP, N] Y;
  for (n in 1:N) Y[ : , n] = y[n];
  real log_2pi = log(2 * pi());

  // Pre-compute 1 - infM for vectorised log_theta calculation
  matrix[nC, nP] infM_inv = 1.0 - infM;

  // Indices of present pathogens in the full pathogen list
  array[nPp] int pres_idx;
  {
    int k = 1;
    for (p in 1:nP) {
      if (pres[p] == 1) {
        pres_idx[k] = p;
        k += 1;
      }
    }
  }

}


parameters {

  array[nS] vector<lower=0, upper=1>[nPp] prev; // infection prevalence per stratum
  real<lower=0> sd0;                             // SD of negative titer distribution
  real<lower=0> sd1;                             // SD of positive titer distribution
  vector[nP] mu0;                                // mean titer, negative
  vector<lower=0>[nPp] mu1;                      // mean titer increase, positive
  vector<lower=0>[(nP * nPp) - nPp] phi;         // relative cross-reactive titer increases
  real<lower=0, upper=1> rho00;                  // correlation between negative titers

}


transformed parameters {

  array[nC] vector[nP] mu;            // Gaussian means per infection combination
  array[nC] vector[nP] sigma;         // Gaussian SDs per infection combination
  matrix[N, nC] pC;                   // log-prob per individual per infection combination
  vector[N] log_lik;                  // individual log-likelihoods
  array[nS] simplex[nC] theta;        // Gaussian mixture weights per stratum
  matrix[nP, nP] CR = rep_matrix(0, nP, nP); // cross-reactive titer increase matrix
  array[nC] matrix[nP, nP] covM;              // covariance matrices per infection combination
  array[nC] cholesky_factor_cov[nP] L;        // Cholesky factors of covariance matrices

  {
    // Local variables
    array[nS] vector[nC] log_theta;
    int ix = 1;
    real cv;
    real var_sum;

    //--- 1. Vectorised log_theta per stratum ---//
    for (s in 1:nS) {
      vector[nP] log_prev_full  = rep_vector(0, nP);
      vector[nP] log1m_prev_full = rep_vector(0, nP);
      for (k in 1:nPp) {
        log_prev_full[pres_idx[k]]  = log(prev[s, k]);
        log1m_prev_full[pres_idx[k]] = log1m(prev[s, k]);
      }
      log_theta[s] = infM * log_prev_full + infM_inv * log1m_prev_full;
      theta[s]     = exp(log_theta[s]);
    }

    //--- 2. Cross-reactivity matrix ---//
    for (p in 1:nPp) {
      for (p2 in 1:nP) {
        if (p == p2)
          CR[p, p2] = 0;
        else {
          CR[p, p2] = phi[ix];
          ix += 1;
        }
      }
    }

    //--- 3. Gaussian means & SDs per infection combination ---//
    sigma[1] = rep_vector(sd0, nP);
    mu[1, ]  = mu0;

    for (c in 2:nC) {
      for (p in 1:npos[c]) {
        sigma[c, wpos[c, p]] = sd1;
        mu[c, wpos[c, p]]    = mu0[wpos[c, p]] + mu1[wpos[c, p]];
      }
      if (npos[c] == 1) {
        for (p in 1:(nP - npos[c])) {
          sigma[c, wneg[c, p]] = sqrt(sd0 ^ 2 + (CR[wpos[c, 1], wneg[c, p]] * sd1) ^ 2);
          mu[c, wneg[c, p]]    = mu0[wneg[c, p]] + CR[wpos[c, 1], wneg[c, p]] * mu1[wpos[c, 1]];
        }
      } else {
        for (p in 1:(nP - npos[c])) {
          mu[c, wneg[c, p]] = mu0[wneg[c, p]];
          var_sum = sd0 ^ 2;
          for (j in 1:npos[c]) {
            mu[c, wneg[c, p]] += CR[wpos[c, j], wneg[c, p]] * mu1[wpos[c, j]];
            var_sum           += (CR[wpos[c, j], wneg[c, p]] * sd1) ^ 2;
          }
          sigma[c, wneg[c, p]] = sqrt(var_sum);
        }
      }
    }

    //--- 4. Covariance matrices & Cholesky decomposition ---//
    for (c in 1:nC) {
      for (p in 1:nP)
        covM[c, p, p] = sigma[c, p] ^ 2;
      for (p in 1:(nP - 1)) {
        for (p2 in (p + 1):nP) {
          if (infM[c, p] + infM[c, p2] == 2) {
            cv = 0;
          } else if (infM[c, p] + infM[c, p2] == 0) {
            cv = rho00 * sd0 ^ 2;
            for (j in 1:npos[c])
              cv += CR[wpos[c, j], p] * CR[wpos[c, j], p2] * sd1 ^ 2;
          } else {
            cv = (infM[c, p] == 1) ? CR[p, p2] * sd1 ^ 2 : CR[p2, p] * sd1 ^ 2;
          }
          covM[c, p, p2] = cv;
          covM[c, p2, p] = cv;
        }
      }
      L[c] = cholesky_decompose(covM[c]);
    }

    //--- 5. Optimised multivariate normal likelihood ---//
    for (c in 1:nC) {
      real constant_part = -0.5 * (nP * log_2pi + 2 * sum(log(diagonal(L[c]))));
      matrix[nP, N] z = mdivide_left_tri_low(L[c], Y - rep_matrix(mu[c], N));
      for (n in 1:N)
        pC[n, c] = log_theta[strat[n], c] + constant_part - 0.5 * dot_self(z[ : , n]);
    }

    for (n in 1:N) log_lik[n] = log_sum_exp(pC[n, : ]);
  }
}


model {

  // Priors
  for (s in 1:nS) prev[s] ~ beta(1, 2);
  mu0   ~ normal(prior_mu0[1], prior_mu0[2]);
  mu1   ~ normal(prior_mu1[1], prior_mu1[2]);
  sd0   ~ normal(0, sigma_scale[1]);
  sd1   ~ normal(0, sigma_scale[2]);
  phi   ~ exponential(5);
  rho00 ~ beta(2, 2);

  // Sum individual likelihoods
  target += sum(log_lik);

}


generated quantities {

  real sumloglik = sum(log_lik);
  vector[nP] prevAll = rep_vector(0, nP); // overall weighted prevalence
  array[nAgg] vector[nP] prevAgg = rep_array(rep_vector(0, nP), nAgg); // weighted prevalence by strata

  for (p in 1:nP) {
    if (pres[p] == 1) {
      int pp_idx = sum(pres[1:p]);

      // collect stratum-level prevalence for pathogen p into a temporary vector
      vector[nS] prevStrat_p;
      for (s in 1:nS) prevStrat_p[s] = prev[s, pp_idx];

      // overall prevalence (weighted by stratum size)
      prevAll[p] = dot_product(NperS, prevStrat_p) / N;

      // aggregated prevalence (weighted combinations)
      for (a in 1:nAgg) prevAgg[a, p] = aggW[a] * prevStrat_p;
    }
  }

}
