// begin stan/stan_raydqf_mod.stan
// direct (Rayleigh) prior, indirect (exponential) likelihood

functions{
  real exponential_ldqf_lpdf(real p, real lambda){
   if (lambda<=0) reject("lambda<=0, found lambda=", lambda);
   if (p>=1) reject("p>=1, found p=",p);
   return log(lambda)+log1m(p);
  }
}
data {
  int<lower=0> N;
  real<lower=0> y[N];
  real ray_sigma;
}
parameters {
  real<lower=1e-15> lambda;
}
model {
  real p[N];
  lambda ~ rayleigh(ray_sigma);
  // transform data into probability (given parameter) and compute likelihood
  for (i in 1:N){
    p[i] = exponential_cdf(y[i], lambda);
    p[i] ~ exponential_ldqf(lambda);
  }
}

// end stan/stan_raydqf_mod.stan
