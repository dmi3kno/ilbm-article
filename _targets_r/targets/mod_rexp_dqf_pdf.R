tar_target(mod_rexp_dqf_pdf, {
  # compile and fit Rayleigh exponential indirect prior direct likelihood model
  stan_model("stan/rexp_dqf_pdf_mod.stan")
})
