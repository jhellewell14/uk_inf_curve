functions {
  // exponential quadratic kernal
  real spd_SE(real alpha, real rho, real w) {
    real S;
    S = (alpha^2) * sqrt(2*pi()) * rho * exp(-0.5*(rho^2)*(w^2));
    return S;
  }

  // basis function for approximate hilbert space gp
  // see here for details: https://arxiv.org/pdf/2004.11408.pdf
  vector phi_SE(real L, int m, vector x) {
    vector[rows(x)] fi;
    fi = 1/sqrt(L) * sin(m*pi()/(2*L) * (x+L));
    return fi;
  }

  // eigenvalues for approximate hilbert space gp
  // see here for details: https://arxiv.org/pdf/2004.11408.pdf
  real lambda(real L, int m) {
    real lam;
    lam = ((m*pi())/(2*L))^2;
    return lam;
  }

  // discretised truncated gamma pmf
  real discretised_gamma_pmf(int y, real mu, real sigma, int max_val) {
    // calculate alpha and beta for gamma distribution
    real c_sigma = sigma + 1e-5;
    real alpha = ((mu)/ c_sigma)^2;
    real beta = (mu) / (c_sigma^2);
    //account for numerical issues
    alpha = alpha <= 0 ? 1e-5 : alpha; //
    beta = beta <= 0 ? 1e-5 : beta;
    alpha = is_inf(alpha) ? 1e8 : alpha;
    beta = is_inf(beta) ? 1e8 : beta;
    return((gamma_cdf(y + 1, alpha, beta) - gamma_cdf(y, alpha, beta)) /
    (gamma_cdf(max_val + 1, alpha, beta) - gamma_cdf(1, alpha, beta)));
  }

  // discretised truncated lognormal pmf
  real discretised_lognormal_pmf(int y, real mu, real sigma, int max_val) {
    real small = 1e-5;
    real adj_y = y + small;
    return((normal_cdf((log(adj_y + 1) - mu) / sigma, 0.0, 1.0) -
    normal_cdf((log(adj_y) - mu) / sigma, 0.0, 1.0)) /
    (normal_cdf((log(max_val + small) - mu) / sigma, 0.0, 1.0) -
    normal_cdf((log(small) - mu) / sigma, 0.0, 1.0)));
  }

  // convolve a pdf and case vector
  vector convolve(vector cases, vector pdf) {
    int t = num_elements(cases);
    int max_pdf = num_elements(pdf);
    vector[t] convolved_cases = rep_vector(1e-5, t);
    for (s in 1:t) {
      convolved_cases[s] += dot_product(cases[max(1, (s - max_pdf + 1)):s], tail(pdf, min(max_pdf, s)));
    }
    return(convolved_cases);
  }
}

data {
  int t; // time steps
  int cases[t];
  int day_of_week[t];
  int <lower = 1> M;
  real L;
  vector[t] time;
  int delays;
  real delay_mean_sd[delays];  // prior sd of mean incubation period
  real delay_mean_mean[delays];// prior mean of mean incubation period
  real delay_sd_mean[delays];  // prior sd of sd of incubation period
  real delay_sd_sd[delays];    // prior sd of sd of incubation period
  int max_delay[delays];       // maximum incubation period
  real gt_mean_sd;                   // prior sd of mean generation time
  real gt_mean_mean;                 // prior mean of mean generation time
  real gt_sd_mean;                   // prior sd of sd of generation time
  real gt_sd_sd;                     // prior sd of sd of generation time
  int max_gt;                        // maximum generation time
  real lengthscale_alpha;            // alpha for gp lengthscale prior
  real lengthscale_beta;             // beta for gp lengthscale prior
  real alpha_sd;                     // standard deviation of the alpha gp kernal parameter
}

transformed data {
  matrix[t, M] PHI;
  // real gt_mean = gt_mean_mean;
  // real gt_sd = gt_sd_mean;
  // real delay_mean[delays] = delay_mean_mean;
  // real delay_sd[delays] = delay_sd_mean;

  for(m in 1:M) {
    PHI[, m] = phi_SE(L, m, time);
  }
}

parameters {
  real<lower = 0> rho; // length scale
  real<lower = 0> alpha; // scale
  vector[M] beta; // eta
  real <lower = 0> phi;
  simplex[7] day_of_week_eff_raw;
  real <lower = 0> gt_mean;
  real <lower = 0> gt_sd;
  real<lower = 0> delay_mean[delays];                 // mean of delays
  real<lower = 0> delay_sd[delays];                   // sd of delays
}

transformed parameters {
  vector[t] log_initial_inf; // f  
  vector[M] diagSPD;
  vector[M] SPD_beta;
  vector[t] infections;
  vector[t] reports_hold;
  vector[t] reports;
  vector[7] day_of_week_eff;

  for(m in 1:M) {
    // Spectral density calculation
    diagSPD[m] =  sqrt(spd_SE(alpha, rho, sqrt(lambda(L, m))));
  }

  // Linear model using the spectral densities
  SPD_beta = diagSPD .* beta;
  log_initial_inf = PHI[,] * SPD_beta;

  infections = exp(log_initial_inf);

  for (s in 1:delays) {
    // reverse the distributions to allow vectorised access
    vector[max_delay[s]] rev_delay = rep_vector(1e-5, max_delay[s]);
    for (j in 1:(max_delay[s])) {
      rev_delay[j] +=
      discretised_lognormal_pmf(max_delay[s] - j, delay_mean[s], delay_sd[s], max_delay[s]);
    }
    if (s == 1) {
      reports_hold = convolve(infections, rev_delay);
    }else{
      reports_hold = convolve(reports_hold, rev_delay);
    }
  }
  reports = reports_hold;

 day_of_week_eff = 7 * day_of_week_eff_raw;

for (s in 1:t) {
  // add reporting effects (adjust for simplex scale)
  reports[s] *= day_of_week_eff[day_of_week[s]];
 }

}

model {
  rho ~ inv_gamma(lengthscale_alpha, lengthscale_beta);
  alpha ~ normal(0, alpha_sd);
  beta ~ std_normal();
  phi ~ exponential(1);
  cases ~ neg_binomial_2(reports, phi);

  for (s in 1:delays) {
    target += normal_lpdf(delay_mean[s] | delay_mean_mean[s], delay_mean_sd[s]) * t;
    target += normal_lpdf(delay_sd[s] | delay_sd_mean[s], delay_sd_sd[s]) * t;
  }
  
  // penalised_prior on generation interval
  target += normal_lpdf(gt_mean | gt_mean_mean, gt_mean_sd) * t;
  target += normal_lpdf(gt_sd | gt_sd_mean, gt_sd_sd) * t;

  // target += normal_lpdf(sum(reports) | sum(cases), sigma);
}
