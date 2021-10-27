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
