functions{ // begin the functions block
vector vlookup(vector x, vector xs, vector ys){ //given x and a grid, returns a linearly interpolated y
  int N =rows(x); // number of ordered!! datapoints
  int M =rows(xs); // number of ordered!! grid cells
  real t;
  real w;
  vector[N] res;
  int i = 1;
  for (j in 2:M){
    while (i<=N && x[i]<=xs[j] && x[i]>=xs[j-1]) {
     res[i] = ys[j-1];
     i+=1;
    }
  }
  return res;
}
real govindarajulu_s_qf(real p, real gamma, real sigma) {
  real inv_gamma = inv(gamma);
  return sigma*gamma*pow(p, gamma)*(1+inv_gamma-p);
}
vector govindarajulu_v_qf(vector p, real gamma, real sigma) {
  int N = rows(p);
  vector[N] res=rep_vector(not_a_number(), N);
  real inv_gamma = inv(gamma);
  for (i in 1:N){
      if(p[i]>0 && p[i]<=1)
         res[i]=exp(log(sigma)+log(gamma)+gamma*log(p[i])+log1p(inv_gamma-p[i]));
    }
  return res;
}
real govindarajulu_s_qdf(real p, real gamma, real sigma) {
  return sigma*gamma*(gamma+1)*pow(p, gamma-1)*(1-p);
}
real govindarajulu_s_ldqf_lpdf(real p, real gamma, real sigma){
  return log(inv(govindarajulu_s_qdf(p, gamma, sigma)));
}
vector govindarajulu_iqf_algebra_system(vector u0, vector params, data real[] x_r, data int[] x_i){
  return [x_r[1] - govindarajulu_s_qf(u0[1], params[1], params[2])]';
}
real approx_govindarajulu_cdf_algebra(data real x, real u_guess, real gamma, real sigma, data real rel_tol, data real f_tol, data real max_steps){
  return algebra_solver(govindarajulu_iqf_algebra_system, [u_guess]', to_vector({gamma, sigma}), {x}, {0},  rel_tol, f_tol, max_steps)[1];
}
real genexp_s_lpdf(real x, real alpha, real lambda) {// generalized exponential log PDF
  return log(alpha)+log(lambda)+(alpha-1)*log1m_exp(-lambda*x)-lambda*x;
  }
} // end of block
data {
  int<lower=0> N; // size of the data
  vector[N] x; // data on variable level
  int<lower=0> M; // size of probability grid
  vector[M] ys_grd; // probability grid
  real genexp_alpha;
  real genexp_lambda;
  real min_sigma;
  real exp_lambda;
  real rel_tol; // for algebra solver
  real f_tol; // for algebra solver
  real max_steps;// for algebra solver
}
transformed data{
  vector[N] x_srt = sort_asc(x);
}
parameters {
  real<lower=1e-6> gamma; // gamma for govindarajulu
  real<lower=1e-6> dsigma; // exponential "tail" of sigma
}

transformed parameters{
   real sigma=dsigma+min_sigma;// sigma  for govindarajulu
}

model {
  vector[N] u;
  // create grid of xs given the parameter
  vector[M] xs_grd = govindarajulu_v_qf(ys_grd, gamma, sigma);
  vector[N] u_guess = vlookup(x_srt, xs_grd, ys_grd);
  //Grids are done. Sampling!
  target += genexp_s_lpdf(gamma | genexp_alpha, genexp_lambda);
  target += exponential_lpdf(dsigma | exp_lambda);
  for (i in 1:N){
   u[i] = approx_govindarajulu_cdf_algebra(x_srt[i], u_guess[i], gamma, sigma, rel_tol, f_tol, max_steps);
   target += govindarajulu_s_ldqf_lpdf(u[i] | gamma, sigma);
  }
}
