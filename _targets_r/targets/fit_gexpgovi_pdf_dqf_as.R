tar_target(fit_gexpgovi_pdf_dqf_as, {
  sampling(mod_gexpgovi_pdf_dqf_as, data=bathtub_data, init = initfun_gexpgovi,
           iter=5000, seed=42, control=list(adapt_delta=0.9))
})
