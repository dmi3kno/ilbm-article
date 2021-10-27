tar_target(claims_data, {
  claims.obs <- c(100, 950, 450)
  # We will be checking against the conjugate model
  post_a <-  gamma_a+length(claims.obs)
  post_b <- gamma_b+sum(claims.obs)
  
  list(N = length(claims.obs), y = claims.obs,
      gamma_a=gamma_a, gamma_b=gamma_b, post_a=post_a, post_b=post_b, ray_sigma=ray_sigma
      )
})
