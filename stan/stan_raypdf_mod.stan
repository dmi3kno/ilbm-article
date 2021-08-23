// begin stan/stan_raypdf_mod.stan

data {
  int<lower=0> N;
  real<lower=0> y[N];
  real ray_sigma;
}
parameters {
  real<lower=1e-15> lambda;
}
model {
  lambda ~ rayleigh(ray_sigma);
  for (n in 1:N)
    y[n] ~ exponential(lambda);
}
// end stan/stan_raypdf_mod.stan
