tar_target(fit_gnh, {
  # sample from g-and-h model using fmcmc. See likelihood function ll_gnh_unif() for details
  MCMC(ll_gnh_unif, initial = make_initials_gnh(),
                     x=gnh_data,
                     A=5, B=5, C=0.8,
                     nsteps = 5000,
                     burnin = 2000,
                     multicore = FALSE, nchains = 4L,
                     kernel = kernel_ram(),
                     progress = TRUE)
})
