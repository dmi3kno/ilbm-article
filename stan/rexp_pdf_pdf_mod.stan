// begin stan/rexp_pdf_pdf_mod.stan
// direct (Rayleigh) prior, direct (exponential) likelihood
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
// end stan/rexp_pdf_pdf_mod.stan
