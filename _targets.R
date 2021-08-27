library(targets)

source("R/targets-functions.R")

tar_option_set(packages = c("qpd", "extraDistr", "tidyverse", "posterior", "rstan", "fmcmc"))
options(mc.cores = parallel::detectCores()-2)
options(clustermq.scheduler="multicode")

list(
  tar_target(name=claims_data, make_claims_data()),
  tar_target(name=mod_gexp_pdf_pdf, stan_model("stan/gexp_pdf_pdf_mod.stan")),
  tar_target(name=fit_gexp_pdf_pdf, sampling(mod_gexp_pdf_pdf, data=claims_data,
                                             seed=42, refresh=0, iter=5000)),
  tar_target(name=mod_gexp_pdf_dqf, stan_model("stan/gexp_pdf_dqf_mod.stan")),
  tar_target(name=fit_gexp_pdf_dqf, sampling(mod_gexp_pdf_dqf, data=claims_data,
                                             init = initfun_gexp,
                                             seed=42, refresh=0, iter=5000)),
  tar_target(name=mod_rexp_pdf_pdf, stan_model("stan/rexp_pdf_pdf_mod.stan")),
  tar_target(name=fit_rexp_pdf_pdf, sampling(mod_rexp_pdf_pdf, data=claims_data,
                                             seed=42, refresh=0, iter=5000)),
  tar_target(name=mod_rexp_pdf_dqf, stan_model("stan/rexp_pdf_dqf_mod.stan")),
  tar_target(name=fit_rexp_pdf_dqf, sampling(mod_rexp_pdf_dqf, data=claims_data,
                                             init = initfun_rexp,
                                             seed=42, refresh=0, iter=5000)),
  tar_target(name=mod_rexp_dqf_pdf, stan_model("stan/rexp_dqf_pdf_mod.stan")),
  tar_target(name=fit_rexp_dqf_pdf, sampling(mod_rexp_dqf_pdf, data=claims_data,
                                             init = initfun_rexp_v,
                                             seed=42, refresh=0, iter=5000)),
  tar_target(name=mod_rexp_dqf_dqf, stan_model("stan/rexp_dqf_dqf_mod.stan")),
  tar_target(name=fit_rexp_dqf_dqf, sampling(mod_rexp_dqf_dqf, data=claims_data,
                                             init = initfun_rexp_v,
                                             seed=42, refresh=0, iter=5000)),
  tar_target(name=bathtub_data, make_bathtub_data()),
  tar_target(name=mod_gexpgovi_pdf_dqf_as, stan_model("stan/gexpgovi_pdf_dqf_as_mod.stan", model_name="somemodel")),
  tar_target(name=fit_gexpgovi_pdf_dqf_as, sampling(mod_gexpgovi_pdf_dqf_as, data=bathtub_data,
                                                    init = initfun_gexpgovi,
                                                    iter=5000, seed=42,
                                                    control=list(adapt_delta=0.9))),
  tar_target(name=gnh_data, subset_gnh_data(100)),
  tar_target(name=fit_gnh, MCMC(ll_gnh_unif,
                                initial = make_initials_gnh(),
                                x=gnh_data,
                                A=5, B=5, C=0.8,
                                nsteps = 5000,
                                burnin = 2000,
                                multicore = FALSE, nchains = 4L,
                                kernel = kernel_ram(),
                                progress = TRUE))
  )

