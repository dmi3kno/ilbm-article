initfun_gexp <- function() list(lambda=gamma_a/gamma_b + runif(1, -2/gamma_b, +2/gamma_b))
initfun_rexp <- function() list(lambda=runif(1,1e-3,5e-3))
initfun_rexp_v <- function() list(v=runif(1,1e-1,9e-1))
