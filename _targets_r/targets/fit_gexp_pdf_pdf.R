tar_target(fit_gexp_pdf_pdf, {
  sampling(mod_gexp_pdf_pdf, data=claims_data, seed=42, refresh=0, iter=5000)
})
