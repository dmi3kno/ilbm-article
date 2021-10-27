tar_target(fit_rexp_dqf_pdf, {
  sampling(mod_rexp_dqf_pdf, data=claims_data, init = initfun_rexp_v, seed=42, refresh=0, iter=5000)
})
