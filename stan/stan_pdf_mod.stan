// begin file stan/stan_pdf_mod.stan
data {
  int<lower=0> N; // number of data points
  real<lower=0> y[N]; // data
  real gamma_a; // hyperparameter alpha
  real gamma_b; // hyperparameter beta
}
parameters {
  real<lower=1e-15> lambda; // strictly positive 
}
model {
  lambda ~ gamma(gamma_a, gamma_b); // prior
  for (n in 1:N)
    y[n] ~ exponential(lambda); // likelihood
}
// end file stan/stan_pdf_mod.stan
