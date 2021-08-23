gamma_a <- 4 # shape
gamma_b <-1/0.001 # rate = 1/scale
ray_sigma <- 0.003

make_claims_data <- function(){
  claims.obs <- c(100, 950, 450)
  # We will be checking against the conjugate model
  post_a <-  gamma_a+length(claims.obs)
  post_b <- gamma_b+sum(claims.obs)

  list(N = length(claims.obs), y = claims.obs,
      gamma_a=gamma_a, gamma_b=gamma_b, post_a=post_a, post_b=post_b, ray_sigma=ray_sigma
      )
}

make_bathtub_data <- function(){
  p_grd <- qpd::make_pgrid(5000)
   vec_times <- c(0.1,7,36,67,84,0.2,11,40,67,84,1,12,45,67,
                  84,1,18,46,67,85,1,18,47,72,85,1,18, 50,
                  75,85,1,18,55,79,85,2,18,60,82,85,3,21,
                  63,82,86,6,32,63,83,86)
  list(N=length(vec_times), x=vec_times, M=length(p_grd), ys_grd=p_grd,
       genexp_alpha=5, genexp_lambda=1, min_sigma=max(vec_times), exp_lambda=0.5,
       rel_tol=1e-12, f_tol=1e-6, max_steps=1e6, verbose=0)
}

initfun_gexp <- function() {
  list(lambda=gamma_a/gamma_b + runif(1, -2/gamma_b, +2/gamma_b))
}

initfun_rexp <- function() list(lambda=runif(1,1e-3,5e-3))

initfun_rexp_v <- function() list(v=runif(1,1e-1,9e-1))

initfun_gexpgovi <- function() list(gamma=runif(1,1,4), dsigma=runif(1,1e-7,1e-1))

#### g-and-h distribuition

make_all_gnh_data <- function(){
  A=5; B=5; C=0.8; g=5; h=0.25
  tibble(p=runif(1e4),
         q=qgnh(p, A,B,C,g,h),
         p_approx=pgnh(q, A,B,C,g,h))
}

subset_gnh_data <- function(df,N){
  A=5; B=5; C=0.8; g=5; h=0.25
  set.seed(42) # correct seed!
  rgnh(N, A, B, C, g, h)
}

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

initials_gnh <- matrix(c(runif(4,2,7), runif(4, 0.1,0.5)),
                   ncol = 2, byrow = FALSE, dimnames = list(NULL, c("g","h")))
#### Wakeby distrib

make_flood_data <- function(){
  # Sources:
  # - The Flood Estimation Handbook flood peak dataset for the period 1937 – 1990 inclusive.
  # - The Highest Instantaneous Flow Series held by the National River Flow Archive CEH Wallingford
  # (supplied by the Environment Agency and its predecessors), for the period 1991 – 1998 inclusive.
  # Citation: Bayliss, A., & Reed, D. W. (2001). The use of historical data in flood frequency estimation: Report to MAFF.
  # Centre for Ecology and Hydrology. Appendix H

  flood_data <- tibble::tribble(
    ~g_n,  ~w_year, ~fl_date, ~peak_flow,
    54002L, 1937L, "01/13/38",  47.021,
    54002L, 1938L, "01/27/39", 240.382,
    54002L, 1939L, "02/08/40", 316.213,
    54002L, 1940L, "11/22/40", 187.123,
    54002L, 1941L, "01/25/42", 183.657,
    54002L, 1942L, "02/01/43", 201.259,
    54002L, 1943L, "01/24/44",   7.574,
    54002L, 1944L, "02/01/45", 103.298,
    54002L, 1945L, "12/29/45",  86.275,
    54002L, 1946L, "03/14/47", 356.187,
    54002L, 1947L, "09/13/48",   67.11,
    54002L, 1948L, "01/02/49",  91.377,
    54002L, 1949L, "02/04/50", 148.908,
    54002L, 1950L, "01/05/51", 181.934,
    54002L, 1951L, "11/09/51", 130.432,
    54002L, 1952L, "12/21/52", 130.432,
    54002L, 1953L, "02/19/54",  86.275,
    54002L, 1954L, "03/27/55", 190.617,
    54002L, 1955L, "01/31/56",  93.851,
    54002L, 1956L, "12/29/56", 138.782,
    54002L, 1957L, "02/25/58", 137.556,
    54002L, 1958L, "01/22/59", 243.687,
    54002L, 1959L, "01/24/60", 245.633,
    54002L, 1960L, "12/04/60", 215.279,
    54002L, 1961L, "01/07/62",   92.29,
    54002L, 1962L, "03/31/63",  67.913,
    54002L, 1963L, "11/19/63", 117.402,
    54002L, 1964L, "03/21/65",  41.032,
    54002L, 1965L, "12/10/65", 148.443,
    54002L, 1966L, "03/10/67",  131.49,
    54002L, 1967L, "07/11/68", 361.909,
    54002L, 1968L, "03/13/69", 198.944,
    54002L, 1969L, "02/20/70",  94.897,
    54002L, 1970L, "01/24/71",   157.4,
    54002L, 1971L, "02/04/72", 188.904,
    54002L, 1972L, "12/07/72", 112.565,
    54002L, 1973L, "02/11/74", 135.722,
    54002L, 1974L, "03/14/75", 172.612,
    54002L, 1975L, "09/26/76",  35.937,
    54002L, 1976L, "06/15/77", 176.653,
    54002L, 1977L, "01/28/78", 123.646,
    54002L, 1978L, "02/02/79", 214.387,
    54002L, 1979L, "12/28/79", 230.596,
    54002L, 1980L, "03/11/81", 215.716,
    54002L, 1981L, "12/30/81", 264.091,
    54002L, 1982L, "05/02/83", 155.035,
    54002L, 1983L, "02/07/84", 102.542,
    54002L, 1984L, "11/24/84", 174.533,
    54002L, 1985L, "01/10/86", 145.447,
    54002L, 1986L, "04/05/87", 128.578,
    54002L, 1987L, "01/24/88", 192.414,
    54002L, 1988L, "04/07/89", 115.592,
    54002L, 1989L, "02/08/90", 163.307,
    54002L, 1990L, "01/10/91", 134.179,
    54002L, 1991L, "01/09/92",   138.8,
    54002L, 1992L, "01/13/93",   212.6,
    54002L, 1993L, "01/05/94",   143.4,
    54002L, 1994L, "01/22/95",   124.3,
    54002L, 1995L, "12/22/95",   113.9,
    54002L, 1996L, "02/26/97",   31.88,
    54002L, 1997L, "04/10/98",     427,
    54002L, 1998L, "01/16/99",   149.7
  )
  flood_data
}
