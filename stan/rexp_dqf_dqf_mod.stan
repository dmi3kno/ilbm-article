// begin stan/rexp_dqf_dqf_mod.stan
// indirect (Rayleigh) prior, indirect (exponential) likelihood

functions{
  real exponential_ldqf_lpdf(real p, real lambda){
   if (lambda<=0) reject("lambda<=0, found lambda=", lambda);
   if (p>=1) reject("p>=1, found p=",p);
   return log(lambda)+log1m(p);
  }
  // not needed. Here for completeness 
  real rayleigh_qf(real p, real sigma){
  if (sigma<=0) reject("sigma<=0; found sigma =", sigma);
  return sigma*sqrt(-2*log1m(p));
  }
  
  real rayleigh_qdf(real p, real sigma){
    if (sigma<=0) reject("sigma<=0; found sigma =", sigma);
    return sigma/(sqrt(2)*sqrt(-log1m(p))*(1-p));
  }
  real rayleigh_s_ldqf_lpdf(real p, real sigma){
    if (sigma<=0) reject("sigma<=0; found sigma =", sigma);
    return 0.5*log(2)+log(sqrt(-log1m(p)))+log1m(p)-log(sigma);
  }

  
} // end of functions block
data {
  int<lower=0> N;
  real<lower=0> y[N];
  real ray_sigma;
}
parameters {
  real<lower=0, upper=1> v; // probability for prior
}
transformed parameters{
 real lambda = rayleigh_qf(v, ray_sigma);
}
model {
  real p[N];
  v ~ rayleigh_s_ldqf(ray_sigma);
  // Jacobian adjustment
  target+=log(abs(rayleigh_qdf(v, ray_sigma)));
  // transform data into probability (given parameter) and compute likelihood
  for (i in 1:N){
    p[i] = exponential_cdf(y[i], lambda);
    p[i] ~ exponential_ldqf(lambda);
  }
}// end of model block

// end stan/rexp_dqf_dqf_mod.stan
