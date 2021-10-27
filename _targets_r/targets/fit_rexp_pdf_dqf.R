tar_target(fit_rexp_pdf_dqf, {
  sampling(mod_rexp_pdf_dqf, data=claims_data, init = initfun_rexp, seed=42, refresh=0, iter=5000)
})
