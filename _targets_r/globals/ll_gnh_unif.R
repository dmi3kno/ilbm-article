# Likelihood function for fmcmc
ll_gnh_unif <- function(pars, A, B, C, x){
  g <- pars[1]
  h <- pars[2]
  if(h<=0) return(-Inf)
  fmcmc::set_userdata(g=g, h=h)
  #if(!qpd::is_qdf_valid(qpd::fgnh, A=A,B=B,C=C,g=g, h=h)) return(-Inf)
  if(!qpd::is_gnh_valid(A=A,B=B,C=C,g=g, h=h, n_grid=1e3)) return(-Inf)
  u <- pgnh(x, A, B, C, g, h)
  log_likelihood <- dqgnh(u, A, B, C, g, h, log=TRUE) # DQF
  log_prior_g <- dnorm(g, mean=3, sd=1, log = TRUE)
  log_prior_h <- extraDistr::drayleigh(h, 0.3, log = TRUE)
  ll <- log_prior_g+log_prior_h+sum(log_likelihood)
  if(!is.finite(ll)) return(-Inf)
  ll
}
