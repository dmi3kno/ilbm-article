tar_target(mod_rexp_pdf_dqf, {
  # compile and fit Rayleigh exponential direct prior indirect likelihood model
  stan_model("stan/rexp_pdf_dqf_mod.stan")
})
