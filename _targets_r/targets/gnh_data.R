tar_target(gnh_data, {
  # prepare synthetic data for g-and-h distribution
  #### g-and-h distribuition
    N <- 100
    A=5; B=5; C=0.8; g=5; h=0.25
    set.seed(42) # correct seed!
    rgnh(N, A, B, C, g, h)
})
