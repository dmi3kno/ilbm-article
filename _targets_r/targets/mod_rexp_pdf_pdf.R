tar_target(mod_rexp_pdf_pdf, {
  # compile and fit Rayleigh exponential direct prior direct likelihood model  
  stan_model("stan/rexp_pdf_pdf_mod.stan")
})
