
<!-- README.md is generated from README.Rmd. Please edit that file -->

# The tenets of indirect inference in Bayesian models

<!-- badges: start -->
<!-- badges: end -->

## Abstract

This paper extends the application of Bayesian inference to probability
distributions defined in terms of quantile function. We introduce the
method of indirect likelihood to be used in the Bayesian models with
sampling distributions defined by the quantile function. We provide
examples and demonstrate the equivalence of this “quantile-based”
(indirect) likelihood to the conventional “density-defined” (direct)
likelihood. We consider practical aspects of the numerical inversion of
quantile function by root-finding required by the indirect likelihood
method. In particular we consider a problem of ensuring the validity of
an arbitrary quantile function with the help of Chebyshev polynomials
and provide useful tips and implementation of these algorithms in Stan
and R. We also extend the same method to propose the definition of an
“indirect prior” and discuss the situations where it can be useful.

For full text of the article, please, refer to
<https://doi.org/10.31219/osf.io/enzgs>  
Submit your (PRE)review at the
<https://prereview.org/preprints/doi-10.31219-osf.io-enzgs>

## How this article is built

Target Markdown is a powerful R Markdown interface for reproducible
analysis pipelines. Please refer to the [respective chapter in the
`targets` user manual](https://books.ropensci.org/targets/markdown.html)
for details. This report, when rendered, will produce the paper in its
entirety. The report can also be run in interactive mode by setting
`tar_interactive` chunk options to `TRUE`.

## Setup

The paper requires several R packages. The `targets` package must be
version 0.5.0.9000 or above. Other packages required by the paper are
listed below. This chunk needs to be run interactively before the report
is rendered.

``` r
remotes::install_github("dmi3kno/qpd")
install.packages(c("targets", "tarchetypes","extraDistr", "tidyverse", "posterior", "rstan", "fmcmc", "details"))
```

First, we load `targets` to activate the specialized `knitr` engine for
Target Markdown. We also remove the `_targets_r` directory which might
have been previously written by non-interactive runs of the report.

``` r
library(targets)
tar_unscript()
```

## Globals

We first define some global options/functions common to all targets.
This loads necessary packages into the target environments and sets up a
multithread environment.

``` r
options(tidyverse.quiet = TRUE)
tar_option_set(packages = c("qpd", "extraDistr", "tidyverse", "posterior", "rstan", "fmcmc"))
options(mc.cores = parallel::detectCores()-2)
options(clustermq.scheduler="multicode")
#> Established _targets.R and _targets_r/globals/ilbm-globals.R.
```

These global variables will be used throughout the article as
hyper-parameters.

``` r
gamma_a <- 4 # shape
gamma_b <-1/0.001 # rate = 1/scale
ray_sigma <- 0.003
#> Established _targets.R and _targets_r/globals/ilbm-globals-vars.R.
```

## Targets

First we prepare a function which will create a dataset for the claims
example. Now we will define a target which will monitor the resulting
dataset. These target chunks will become the targets with the chunk name
becoming the target name.

``` r
tar_target(claims_data, {
  claims.obs <- c(100, 950, 450)
  # We will be checking against the conjugate model
  post_a <-  gamma_a+length(claims.obs)
  post_b <- gamma_b+sum(claims.obs)
  
  list(N = length(claims.obs), y = claims.obs,
      gamma_a=gamma_a, gamma_b=gamma_b, post_a=post_a, post_b=post_b, ray_sigma=ray_sigma
      )
})
#> Defined target claims_data automatically from chunk code.
#> Established _targets.R and _targets_r/targets/claims_data.R.
```

These initialization functions will be used in initializing the STAN
models

``` r
initfun_gexp <- function() list(lambda=gamma_a/gamma_b + runif(1, -2/gamma_b, +2/gamma_b))
initfun_rexp <- function() list(lambda=runif(1,1e-3,5e-3))
initfun_rexp_v <- function() list(v=runif(1,1e-1,9e-1))
#> Established _targets.R and _targets_r/globals/ilbm-globals-ifuns1.R.
```

### Gamma-Exponential models

Compile and fit gamma exponential direct prior direct likelihood model

``` r
tar_target(mod_gexp_pdf_pdf, {
  stan_model("stan/gexp_pdf_pdf_mod.stan")
})
#> Defined target mod_gexp_pdf_pdf automatically from chunk code.
#> Established _targets.R and _targets_r/targets/mod_gexp_pdf_pdf.R.
```

``` r
tar_target(fit_gexp_pdf_pdf, {
  sampling(mod_gexp_pdf_pdf, data=claims_data, seed=42, refresh=0, iter=5000)
})
#> Defined target fit_gexp_pdf_pdf automatically from chunk code.
#> Established _targets.R and _targets_r/targets/fit_gexp_pdf_pdf.R.
```

Compile and fit gamma exponential direct prior indirect likelihood model

``` r
tar_target(mod_gexp_pdf_dqf, {
  stan_model("stan/gexp_pdf_dqf_mod.stan")
})
#> Defined target mod_gexp_pdf_dqf automatically from chunk code.
#> Established _targets.R and _targets_r/targets/mod_gexp_pdf_dqf.R.
```

``` r
tar_target(fit_gexp_pdf_dqf, {
  sampling(mod_gexp_pdf_dqf, data=claims_data, init = initfun_gexp,  seed=42, refresh=0, iter=5000)
})
#> Defined target fit_gexp_pdf_dqf automatically from chunk code.
#> Established _targets.R and _targets_r/targets/fit_gexp_pdf_dqf.R.
```

### Reyleigh-Exponential models

Compile and fit Rayleigh exponential direct prior direct likelihood
model

``` r
tar_target(mod_rexp_pdf_pdf, {
  # compile and fit Rayleigh exponential direct prior direct likelihood model  
  stan_model("stan/rexp_pdf_pdf_mod.stan")
})
#> Defined target mod_rexp_pdf_pdf automatically from chunk code.
#> Established _targets.R and _targets_r/targets/mod_rexp_pdf_pdf.R.
```

``` r
tar_target(fit_rexp_pdf_pdf, {
  sampling(mod_rexp_pdf_pdf, data=claims_data, seed=42, refresh=0, iter=5000)
})
#> Defined target fit_rexp_pdf_pdf automatically from chunk code.
#> Established _targets.R and _targets_r/targets/fit_rexp_pdf_pdf.R.
```

Compile and fit Rayleigh exponential direct prior indirect likelihood
model

``` r
tar_target(mod_rexp_pdf_dqf, {
  # compile and fit Rayleigh exponential direct prior indirect likelihood model
  stan_model("stan/rexp_pdf_dqf_mod.stan")
})
#> Defined target mod_rexp_pdf_dqf automatically from chunk code.
#> Established _targets.R and _targets_r/targets/mod_rexp_pdf_dqf.R.
```

``` r
tar_target(fit_rexp_pdf_dqf, {
  sampling(mod_rexp_pdf_dqf, data=claims_data, init = initfun_rexp, seed=42, refresh=0, iter=5000)
})
#> Defined target fit_rexp_pdf_dqf automatically from chunk code.
#> Established _targets.R and _targets_r/targets/fit_rexp_pdf_dqf.R.
```

Compile and fit Rayleigh exponential indirect prior direct likelihood
model

``` r
tar_target(mod_rexp_dqf_pdf, {
  # compile and fit Rayleigh exponential indirect prior direct likelihood model
  stan_model("stan/rexp_dqf_pdf_mod.stan")
})
#> Defined target mod_rexp_dqf_pdf automatically from chunk code.
#> Established _targets.R and _targets_r/targets/mod_rexp_dqf_pdf.R.
```

``` r
tar_target(fit_rexp_dqf_pdf, {
  sampling(mod_rexp_dqf_pdf, data=claims_data, init = initfun_rexp_v, seed=42, refresh=0, iter=5000)
})
#> Defined target fit_rexp_dqf_pdf automatically from chunk code.
#> Established _targets.R and _targets_r/targets/fit_rexp_dqf_pdf.R.
```

Compile and fit Rayleigh exponential indirect prior indirect likelihood
model

``` r
tar_target(mod_rexp_dqf_dqf, {
  # compile and fit Rayleigh exponential indirect prior indirect likelihood model
  stan_model("stan/rexp_dqf_dqf_mod.stan")
})
#> Defined target mod_rexp_dqf_dqf automatically from chunk code.
#> Established _targets.R and _targets_r/targets/mod_rexp_dqf_dqf.R.
```

``` r
tar_target(fit_rexp_dqf_dqf, {
  sampling(mod_rexp_dqf_dqf, data=claims_data, init = initfun_rexp_v, seed=42, refresh=0, iter=5000)
})
#> Defined target fit_rexp_dqf_dqf automatically from chunk code.
#> Established _targets.R and _targets_r/targets/fit_rexp_dqf_dqf.R.
```

### Generalized Exponential-Govindarajulu model

``` r
initfun_gexpgovi <- function() list(gamma=runif(1,1,4), dsigma=runif(1,1e-5,1e-3))
#> Established _targets.R and _targets_r/globals/ilbm-globals-ifuns2.R.
```

Define a helper target which will prepare a bathtub data

``` r
tar_target(bathtub_data, {
  set.seed(42)
  p_grd <- qpd::make_pgrid(5000)
  vec_times <- c(0.1,7,36,67,84,0.2,11,40,67,84,1,12,45,67,
                84,1,18,46,67,85,1,18,47,72,85,1,18, 50,
                75,85,1,18,55,79,85,2,18,60,82,85,3,21,
                63,82,86,6,32,63,83,86)
  list(N=length(vec_times), x=vec_times, M=length(p_grd), ys_grd=p_grd,
         genexp_alpha=5, genexp_lambda=1, min_sigma=max(vec_times), exp_lambda=0.5,
         rel_tol=1e-12, f_tol=1e-6, max_steps=1e6, verbose=0)
})
#> Defined target bathtub_data automatically from chunk code.
#> Established _targets.R and _targets_r/targets/bathtub_data.R.
```

Compile and fit genexp-Govindarajulu direct prior indirect likelihood
model

``` r
tar_target(mod_gexpgovi_pdf_dqf_as, {
  # compile and fit genexp-Govindarajulu direct prior indirect likelihood model
  stan_model("stan/gexpgovi_pdf_dqf_as_mod.stan")
})
#> Defined target mod_gexpgovi_pdf_dqf_as automatically from chunk code.
#> Established _targets.R and _targets_r/targets/mod_gexpgovi_pdf_dqf_as.R.
```

``` r
tar_target(fit_gexpgovi_pdf_dqf_as, {
  sampling(mod_gexpgovi_pdf_dqf_as, data=bathtub_data, init = initfun_gexpgovi,
           iter=5000, seed=42, control=list(adapt_delta=0.9))
})
#> Defined target fit_gexpgovi_pdf_dqf_as automatically from chunk code.
#> Established _targets.R and _targets_r/targets/fit_gexpgovi_pdf_dqf_as.R.
```

### g-and-h example

Initialization function for the g-and-h model

``` r
# initial values for the g-and-h model
make_initials_gnh <- function(){
  set.seed(42)
  matrix(c(runif(4,2,7), runif(4, 0.1,0.5)),
                ncol = 2, byrow = FALSE, dimnames = list(NULL, c("g","h")))
}
#> Established _targets.R and _targets_r/globals/ilbm-globals-ifuns3.R.
```

The following target generates the synthethic data

``` r
tar_target(gnh_data, {
  # prepare synthetic data for g-and-h distribution
  #### g-and-h distribuition
    N <- 100
    A=5; B=5; C=0.8; g=5; h=0.25
    set.seed(42) # correct seed!
    rgnh(N, A, B, C, g, h)
})
#> Defined target gnh_data automatically from chunk code.
#> Established _targets.R and _targets_r/targets/gnh_data.R.
```

Define likelihood function for g-and-h model. We will place it into the
global target environement to be available for the MCMC target

``` r
# Likelihood function for fmcmc
ll_gnh_unif <- function(pars, A, B, C, x){
  g <- pars[1]
  h <- pars[2]
  if(h<=0) return(-Inf)
  fmcmc::set_userdata(g=g, h=h)
  #if(!qpd::is_qdf_valid(qpd::fgnh, A=A,B=B,C=C,g=g, h=h)) return(-Inf)
  if(!qpd::is_gnh_valid(A=A,B=B,C=C,g=g, h=h, n_grid=1e3)) return(-Inf)
  u <- pgnh(x, A, B, C, g, h)
  log_likelihood <- dqgnh(u, A, B, C, g, h, log=TRUE) # DQF
  log_prior_g <- dnorm(g, mean=3, sd=1, log = TRUE)
  log_prior_h <- extraDistr::drayleigh(h, 0.3, log = TRUE)
  ll <- log_prior_g+log_prior_h+sum(log_likelihood)
  if(!is.finite(ll)) return(-Inf)
  ll
}
#> Established _targets.R and _targets_r/globals/ll_gnh_unif.R.
```

Sample from g-and-h model using fmcmc. See likelihood function
`ll_gnh_unif()` above for details

``` r
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
#> Defined target fit_gnh automatically from chunk code.
#> Established _targets.R and _targets_r/targets/fit_gnh.R.
```

### Render the article

The article will be rendered automatically together with other targets.
Of course it can be also rendered manually from the `ilbm_article.Rmd`,
but this way we are guaranteed to re-render the article every time the
changes are made to any of the target components.

``` r
#tarchetypes::tar_render(report, "ilbm_article.Rmd")
tarchetypes::tar_render(report, "BA-submission.Rmd")
#> Established _targets.R and _targets_r/targets/produce-rmd.R.
```

## Execute all targets

The following code will remove the old logs, re-build the outdated
targets and re-render the final paper. As a result a new directory
`_targets` will be created on your computer and the hashed version of
the results (in compressed format) will be stored there.

If you inspect the scripts in the `_targets_r/` you will see a list of
statements like `tar_target(name=claims_data,  {...})`. These targets
will be ran each individually, as the sub-tasks (also called targets).

``` r
unlink("logs", recursive = TRUE)
targets::tar_make()
#> ✓ skip target mod_rexp_dqf_dqf
#> ✓ skip target claims_data
#> ✓ skip target mod_rexp_pdf_dqf
#> ✓ skip target mod_gexp_pdf_dqf
#> ✓ skip target bathtub_data
#> ✓ skip target gnh_data
#> ✓ skip target mod_rexp_dqf_pdf
#> ✓ skip target mod_gexpgovi_pdf_dqf_as
#> ✓ skip target mod_rexp_pdf_pdf
#> ✓ skip target mod_gexp_pdf_pdf
#> ✓ skip target fit_rexp_dqf_dqf
#> ✓ skip target fit_rexp_pdf_dqf
#> ✓ skip target fit_gexp_pdf_dqf
#> ✓ skip target fit_gnh
#> ✓ skip target fit_rexp_dqf_pdf
#> ✓ skip target fit_gexpgovi_pdf_dqf_as
#> ✓ skip target fit_rexp_pdf_pdf
#> ✓ skip target fit_gexp_pdf_pdf
#> ✓ skip target report
#> ✓ skip pipeline
```

The `targets` dependency graph helps your readers understand the steps
of the pipeline at a high level. You can review the graph of task
dependencies, including the status of individual nodes, with

``` r
targets::tar_visnetwork()
```

The nodes get invalidated if any modification is made of the node or any
of its parents (ancestors). The out-of-date nodes will be rerun next
time you run `tar_make()`.

At this point, you can go back and run `{targets}` chunks in interactive
mode without interfering with the code or data of the non-interactive
pipeline.

``` r
tar_progress_summary()
#> # A tibble: 1 × 6
#>   skipped started built errored canceled since        
#>     <int>   <int> <int>   <int>    <int> <chr>        
#> 1      19       0     0       0        0 0.046 seconds
```

The results can be retrieved by the name of the subtask. For example,
this will retrieve the fit object for the Genexp-Govindarajulu model
from the article.

``` r
tar_read(fit_gexpgovi_pdf_dqf_as)
#> Inference for Stan model: anon_model.
#> 4 chains, each with iter=5000; warmup=2500; thin=1; 
#> post-warmup draws per chain=2500, total post-warmup draws=10000.
#> 
#>           mean se_mean   sd    2.5%     25%     50%     75%   97.5% n_eff Rhat
#> gamma     2.04    0.00 0.32    1.51    1.81    2.01    2.23    2.72  5161    1
#> dsigma    0.02    0.00 0.08    0.00    0.00    0.00    0.01    0.20  6622    1
#> sigma    86.02    0.00 0.08   86.00   86.00   86.00   86.01   86.20  6622    1
#> lp__   -209.78    0.01 0.87 -212.09 -210.05 -209.48 -209.20 -209.07  3544    1
#> 
#> Samples were drawn using NUTS(diag_e) at Tue Oct 26 15:00:22 2021.
#> For each parameter, n_eff is a crude measure of effective sample size,
#> and Rhat is the potential scale reduction factor on split chains (at 
#> convergence, Rhat=1).
```

## Session Info

``` r
sessioninfo::session_info()%>%
  details::details(summary = 'Session Info')
```

<details closed>
<summary>
<span title="Click to Expand"> Session Info </span>
</summary>

``` r
─ Session info ───────────────────────────────────────────────────────────────
 setting  value                       
 version  R version 4.1.2 (2021-11-01)
 os       Ubuntu 20.04.3 LTS          
 system   x86_64, linux-gnu           
 ui       X11                         
 language (EN)                        
 collate  en_US.UTF-8                 
 ctype    en_US.UTF-8                 
 tz       Europe/Stockholm            
 date     2022-01-03                  

─ Packages ───────────────────────────────────────────────────────────────────
 package        * version    date       lib source                            
 abind            1.4-5      2016-07-21 [1] CRAN (R 4.1.0)                    
 assertthat       0.2.1      2019-03-21 [1] CRAN (R 4.1.0)                    
 backports        1.3.0      2021-10-27 [1] CRAN (R 4.1.2)                    
 broom            0.7.7      2021-06-13 [1] CRAN (R 4.1.0)                    
 callr            3.7.0      2021-04-20 [1] CRAN (R 4.1.0)                    
 cellranger       1.1.0      2016-07-27 [1] CRAN (R 4.1.0)                    
 checkmate        2.0.0      2020-02-06 [1] CRAN (R 4.1.0)                    
 cli              3.1.0      2021-10-27 [1] CRAN (R 4.1.2)                    
 clipr            0.7.1      2020-10-08 [1] CRAN (R 4.1.0)                    
 coda             0.19-4     2020-09-30 [1] CRAN (R 4.1.0)                    
 codetools        0.2-18     2020-11-04 [4] CRAN (R 4.0.3)                    
 colorspace       2.0-2      2021-06-24 [1] CRAN (R 4.1.1)                    
 crayon           1.4.2      2021-10-29 [1] CRAN (R 4.1.2)                    
 curl             4.3.2      2021-06-23 [1] CRAN (R 4.1.0)                    
 data.table       1.14.2     2021-09-27 [1] CRAN (R 4.1.2)                    
 DBI              1.1.1      2021-01-15 [1] CRAN (R 4.1.0)                    
 dbplyr           2.1.1      2021-04-06 [1] CRAN (R 4.1.0)                    
 desc             1.4.0      2021-09-28 [1] CRAN (R 4.1.1)                    
 details        * 0.3.0      2021-10-13 [1] Github (yonicd/details@04b2b4a)   
 digest           0.6.28     2021-09-23 [1] CRAN (R 4.1.1)                    
 distributional   0.2.2      2021-02-02 [1] CRAN (R 4.1.0)                    
 dplyr          * 1.0.6      2021-05-05 [1] CRAN (R 4.1.0)                    
 ellipsis         0.3.2      2021-04-29 [1] CRAN (R 4.1.0)                    
 evaluate         0.14       2019-05-28 [1] CRAN (R 4.1.0)                    
 extraDistr     * 1.9.1      2020-09-07 [1] CRAN (R 4.1.0)                    
 fansi            0.5.0      2021-05-25 [1] CRAN (R 4.1.0)                    
 farver           2.1.0      2021-02-28 [1] CRAN (R 4.1.0)                    
 fastmap          1.1.0      2021-01-25 [1] CRAN (R 4.1.0)                    
 fmcmc          * 0.5-0      2021-07-21 [1] Github (USCbiostats/fmcmc@f9e3a07)
 forcats        * 0.5.1      2021-01-27 [1] CRAN (R 4.1.0)                    
 fs               1.5.0      2020-07-31 [1] CRAN (R 4.1.0)                    
 generics         0.1.1      2021-10-25 [1] CRAN (R 4.1.2)                    
 ggplot2        * 3.3.5      2021-06-25 [1] CRAN (R 4.1.1)                    
 glue             1.6.0      2021-12-17 [1] CRAN (R 4.1.2)                    
 gridExtra        2.3        2017-09-09 [1] CRAN (R 4.1.0)                    
 gtable           0.3.0      2019-03-25 [1] CRAN (R 4.1.0)                    
 haven            2.4.1      2021-04-23 [1] CRAN (R 4.1.0)                    
 hms              1.1.0      2021-05-17 [1] CRAN (R 4.1.0)                    
 htmltools        0.5.2      2021-08-25 [1] CRAN (R 4.1.1)                    
 httr             1.4.2      2020-07-20 [1] CRAN (R 4.1.0)                    
 igraph           1.2.6      2020-10-06 [1] CRAN (R 4.1.0)                    
 inline           0.3.19     2021-05-31 [1] CRAN (R 4.1.0)                    
 jsonlite         1.7.2      2020-12-09 [1] CRAN (R 4.1.0)                    
 knitr            1.36       2021-09-29 [1] CRAN (R 4.1.1)                    
 lattice          0.20-45    2021-09-22 [4] CRAN (R 4.1.1)                    
 lifecycle        1.0.1      2021-09-24 [1] CRAN (R 4.1.2)                    
 loo              2.4.1      2020-12-09 [1] CRAN (R 4.1.0)                    
 lubridate        1.8.0      2021-10-07 [1] CRAN (R 4.1.2)                    
 magrittr       * 2.0.1      2020-11-17 [1] CRAN (R 4.1.0)                    
 MASS             7.3-54     2021-05-03 [4] CRAN (R 4.0.5)                    
 matrixStats      0.61.0     2021-09-17 [1] CRAN (R 4.1.2)                    
 modelr           0.1.8      2020-05-19 [1] CRAN (R 4.1.0)                    
 munsell          0.5.0      2018-06-12 [1] CRAN (R 4.1.0)                    
 pillar           1.6.4      2021-10-18 [1] CRAN (R 4.1.2)                    
 pkgbuild         1.2.0      2020-12-15 [1] CRAN (R 4.1.0)                    
 pkgconfig        2.0.3      2019-09-22 [1] CRAN (R 4.1.0)                    
 png              0.1-7      2013-12-03 [1] CRAN (R 4.1.0)                    
 posterior      * 1.1.0      2021-09-09 [1] CRAN (R 4.1.2)                    
 prettyunits      1.1.1      2020-01-24 [1] CRAN (R 4.1.0)                    
 processx         3.5.2      2021-04-30 [1] CRAN (R 4.1.0)                    
 ps               1.6.0      2021-02-28 [1] CRAN (R 4.1.0)                    
 purrr          * 0.3.4      2020-04-17 [1] CRAN (R 4.1.0)                    
 qpd            * 0.0.0.9000 2021-10-29 [1] local                             
 R6               2.5.1      2021-08-19 [1] CRAN (R 4.1.1)                    
 Rcpp             1.0.7      2021-07-07 [1] CRAN (R 4.1.0)                    
 RcppParallel     5.1.4      2021-05-04 [1] CRAN (R 4.1.0)                    
 readr          * 2.0.0      2021-07-20 [1] CRAN (R 4.1.0)                    
 readxl           1.3.1      2019-03-13 [1] CRAN (R 4.1.0)                    
 reprex           2.0.0      2021-04-02 [1] CRAN (R 4.1.0)                    
 rlang            0.4.12     2021-10-18 [1] CRAN (R 4.1.2)                    
 rmarkdown        2.9        2021-06-15 [1] CRAN (R 4.1.0)                    
 rprojroot        2.0.2      2020-11-15 [1] CRAN (R 4.1.0)                    
 rstan          * 2.26.1     2021-06-15 [1] local                             
 rstudioapi       0.13       2020-11-12 [1] CRAN (R 4.1.0)                    
 rvest            1.0.2      2021-10-16 [1] CRAN (R 4.1.2)                    
 scales           1.1.1      2020-05-11 [1] CRAN (R 4.1.0)                    
 sessioninfo      1.1.1      2018-11-05 [1] CRAN (R 4.1.0)                    
 StanHeaders    * 2.26.1     2021-06-15 [1] local                             
 stringi          1.7.6      2021-11-29 [1] CRAN (R 4.1.2)                    
 stringr        * 1.4.0      2019-02-10 [1] CRAN (R 4.1.0)                    
 tarchetypes    * 0.3.0      2021-08-04 [1] CRAN (R 4.1.1)                    
 targets        * 0.6.0      2021-07-21 [1] CRAN (R 4.1.1)                    
 tensorA          0.36.2     2020-11-19 [1] CRAN (R 4.1.0)                    
 tibble         * 3.1.6      2021-11-07 [1] CRAN (R 4.1.2)                    
 tidyr          * 1.1.3      2021-03-03 [1] CRAN (R 4.1.0)                    
 tidyselect       1.1.1      2021-04-30 [1] CRAN (R 4.1.0)                    
 tidyverse      * 1.3.1      2021-04-15 [1] CRAN (R 4.1.2)                    
 tzdb             0.1.2      2021-07-20 [1] CRAN (R 4.1.0)                    
 utf8             1.2.2      2021-07-24 [1] CRAN (R 4.1.0)                    
 V8               3.6.0      2021-11-10 [1] CRAN (R 4.1.2)                    
 vctrs            0.3.8      2021-04-29 [1] CRAN (R 4.1.0)                    
 withr            2.4.2      2021-04-18 [1] CRAN (R 4.1.0)                    
 xfun             0.26       2021-09-14 [1] CRAN (R 4.1.1)                    
 xml2             1.3.3      2021-11-30 [1] CRAN (R 4.1.2)                    
 yaml             2.2.1      2020-02-01 [1] CRAN (R 4.1.0)                    

[1] /home/dm0737pe/R/x86_64-pc-linux-gnu-library/4.1
[2] /usr/local/lib/R/site-library
[3] /usr/lib/R/site-library
[4] /usr/lib/R/library
```

</details>

<br>
