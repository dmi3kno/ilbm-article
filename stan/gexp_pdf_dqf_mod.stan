// begin stan/gexp_pdf_dqf_mod.stan
functions{
  real exponential_ldqf_lpdf(real p, real lambda){
   if (p>=1) reject("p>=1, found p=",p);
   return log(lambda)+log1m(p);
  }
}
data {
  int<lower=0> N;
  real<lower=0> y[N];
  real gamma_a;
  real gamma_b;
}
parameters {
  real<lower=1e-15> lambda;
}
model {
  real u[N];
  lambda ~ gamma(gamma_a, gamma_b);
  // transform data into probability (given parameter) 
  // likelihood
  for (i in 1:N){
    u[i] = exponential_cdf(y[i], lambda);
    u[i] ~ exponential_ldqf(lambda);
  }
}
// end stan/gexp_pdf_dqf_mod.stan
