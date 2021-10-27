tar_target(fit_rexp_pdf_pdf, {
  sampling(mod_rexp_pdf_pdf, data=claims_data, seed=42, refresh=0, iter=5000)
})
