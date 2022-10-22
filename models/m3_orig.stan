// generated with brms 2.17.0
functions {
  /* cumulative-logit log-PDF for a single response
   * Args:
   *   y: response category
   *   mu: latent mean parameter
   *   disc: discrimination parameter
   *   thres: ordinal thresholds
   * Returns:
   *   a scalar to be added to the log posterior
   */
   real cumulative_logit_lpmf(int y, real mu, real disc, vector thres) {
     int nthres = num_elements(thres);
     if (y == 1) {
       return log_inv_logit(disc * (thres[1] - mu));
     } else if (y == nthres + 1) {
       return log1m_inv_logit(disc * (thres[nthres] - mu));
     } else {
       return log_diff_exp(
         log_inv_logit(disc * (thres[y] - mu)), 
         log_inv_logit(disc * (thres[y - 1] - mu))
       );
     }
   }
  /* cumulative-logit log-PDF for a single response and merged thresholds
   * Args:
   *   y: response category
   *   mu: latent mean parameter
   *   disc: discrimination parameter
   *   thres: vector of merged ordinal thresholds
   *   j: start and end index for the applid threshold within 'thres'
   * Returns:
   *   a scalar to be added to the log posterior
   */
   real cumulative_logit_merged_lpmf(int y, real mu, real disc, vector thres, int[] j) {
     return cumulative_logit_lpmf(y | mu, disc, thres[j[1]:j[2]]);
   }
  /* ordered-logistic log-PDF for a single response and merged thresholds
   * Args:
   *   y: response category
   *   mu: latent mean parameter
   *   thres: vector of merged ordinal thresholds
   *   j: start and end index for the applid threshold within 'thres'
   * Returns:
   *   a scalar to be added to the log posterior
   */
   real ordered_logistic_merged_lpmf(int y, real mu, vector thres, int[] j) {
     return ordered_logistic_lpmf(y | mu, thres[j[1]:j[2]]);
   }
}
data {
  int<lower=1> N;  // total number of observations
  int Y[N];  // response variable
  int<lower=2> nthres;  // number of thresholds
  int<lower=1> K;  // number of population-level effects
  matrix[N, K] X;  // population-level design matrix
  // data for group-level effects of ID 1
  int<lower=1> N_1;  // number of grouping levels
  int<lower=1> M_1;  // number of coefficients per level
  int<lower=1> J_1[N];  // grouping indicator per observation
  // group-level predictor values
  vector[N] Z_1_1;
  // data for group-level effects of ID 2
  int<lower=1> N_2;  // number of grouping levels
  int<lower=1> M_2;  // number of coefficients per level
  int<lower=1> J_2[N];  // grouping indicator per observation
  // group-level predictor values
  vector[N] Z_2_1;
  int prior_only;  // should the likelihood be ignored?
}
transformed data {
  int Kc = K;
  matrix[N, Kc] Xc;  // centered version of X
  vector[Kc] means_X;  // column means of X before centering
  for (i in 1:K) {
    means_X[i] = mean(X[, i]);
    Xc[, i] = X[, i] - means_X[i];
  }
}
parameters {
  vector[Kc] b;  // population-level effects
  ordered[nthres] Intercept;  // temporary thresholds for centered predictors
  vector<lower=0>[M_1] sd_1;  // group-level standard deviations
  vector[N_1] z_1[M_1];  // standardized group-level effects
  vector<lower=0>[M_2] sd_2;  // group-level standard deviations
  vector[N_2] z_2[M_2];  // standardized group-level effects
}
transformed parameters {
  real disc = 1;  // discrimination parameters
  vector[N_1] r_1_1;  // actual group-level effects
  vector[N_2] r_2_1;  // actual group-level effects
  real lprior = 0;  // prior contributions to the log posterior
  r_1_1 = (sd_1[1] * (z_1[1]));
  r_2_1 = (sd_2[1] * (z_2[1]));
  lprior += normal_lpdf(b | 0, 0.5);
  lprior += student_t_lpdf(Intercept | 3, 0, 2.5);
  lprior += student_t_lpdf(sd_1 | 3, 0, 2.5)
    - 1 * student_t_lccdf(0 | 3, 0, 2.5);
  lprior += student_t_lpdf(sd_2 | 3, 0, 2.5)
    - 1 * student_t_lccdf(0 | 3, 0, 2.5);
}
model {
  // likelihood including constants
  if (!prior_only) {
    // initialize linear predictor term
    vector[N] mu = Xc * b;
    for (n in 1:N) {
      // add more terms to the linear predictor
      mu[n] += r_1_1[J_1[n]] * Z_1_1[n] + r_2_1[J_2[n]] * Z_2_1[n];
    }
    for (n in 1:N) {
      target += ordered_logistic_lpmf(Y[n] | mu[n], Intercept);
    }
  }
  // priors including constants
  target += lprior;
  target += std_normal_lpdf(z_1[1]);
  target += std_normal_lpdf(z_2[1]);
}
generated quantities {
  // compute actual thresholds
  vector[nthres] b_Intercept = Intercept + dot_product(means_X, b);
  // additionally sample draws from priors
  real prior_b = normal_rng(0,0.5);
  real prior_Intercept = student_t_rng(3,0,2.5);
  real prior_sd_1 = student_t_rng(3,0,2.5);
  real prior_sd_2 = student_t_rng(3,0,2.5);
  // use rejection sampling for truncated priors
  while (prior_sd_1 < 0) {
    prior_sd_1 = student_t_rng(3,0,2.5);
  }
  while (prior_sd_2 < 0) {
    prior_sd_2 = student_t_rng(3,0,2.5);
  }
}
