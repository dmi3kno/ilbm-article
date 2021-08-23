// begin stan/gexpgovi_pdf_dqf_mod.stan
// direct (Generalized Exponential) prior, indirect (Govindarajulu) likelihood

functions{ // begin the functions block
matrix lookup_v_bracket(vector x, vector xs, vector ys){
  int N =rows(x); // number of ordered!! datapoints
  int M =rows(xs); // number of ordered!! grid cells
  real t;
  real w;
  matrix[N,2] res; // bracket of u values for each observation
  int i = 1;
  for (j in 2:M){ // visiting each "cell" of the grid
    while (i<=N && x[i]>=xs[j-1] && x[i]<=xs[j] ) { // if there are observations falling between
     res[i,] = to_row_vector(ys[(j-1):j]); // for all observations bracketed by the pair of values
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

real genexp_s_lpdf(real x, real alpha, real lambda) {// generalized exponential log PDF
  return log(alpha)+log(lambda)+(alpha-1)*log1m_exp(-lambda*x)-lambda*x;
  }
// only correct from here
  real drootfun_du(real u, real gamma, real sigma){
    return govindarajulu_s_qdf(u, gamma, sigma);
  }
  real rootfun(real u, real x, real gamma, real sigma){
    return x-govindarajulu_s_qf(u,  gamma, sigma);
  }
// only correct down to here
  real iqf_bisection(row_vector u, real max_steps, real rel_tol, real x, real gamma, real sigma, int verbose){
    int steps_taken=0;
    real i_lower=u[1];
    real i_upper=u[2];
    while (steps_taken < max_steps){
      real m = (i_lower+i_upper) / 2.0;
      if (m==0 || fabs(i_upper-i_lower)<rel_tol){
        if(verbose) print("Steps taken ", steps_taken);
        return m;
      }
      if (rootfun(m, x,  gamma,  sigma)>0){ //rootfinding function is called once
        i_lower=m;
      } else {
        i_upper=m;
      }
      steps_taken+=1;
    }
    if(verbose) print("Steps taken ", steps_taken);
    return (i_lower+i_upper)/2.0;
  }

  real iqf_newton_raphson(row_vector u, real max_steps, real rel_tol, real x, real gamma, real sigma, int verbose){
    int steps_taken=0;
    real u_i = u[1];
    while(steps_taken<max_steps && fabs(rootfun(u_i, x, gamma,  sigma))>rel_tol){
      real drf=drootfun_du(u_i, gamma,  sigma);
      u_i+=rootfun(u_i, x, gamma,  sigma)/drf;
      steps_taken +=1;
    }
    if(verbose) print("Steps taken ", steps_taken);
    return u_i;
  }

  real iqf_secant(row_vector u, real max_steps, real rel_tol, real x, real gamma, real sigma, int verbose){
    int steps_taken=1;
    real u0=u[1];
    real u1=u[2];
    real u2=0;
    while (steps_taken<max_steps && fabs(u1-u0)>rel_tol){
      real f1=rootfun(u1, x, gamma,  sigma);
      real f0=rootfun(u0, x, gamma,  sigma);
      u2 = u1-(f1*(u1-u0)/(f1-f0));
      u1=u2;
      u0=u1;
      steps_taken+=1;
    }
    if(verbose) print("Steps taken ", steps_taken);
    return u2;
  }

  real iqf_lagrangepoly(row_vector u, real max_steps, real rel_tol, real x, real gamma, real sigma, int verbose){
    real u0=u[1];
    real u1=(u[1]+u[2])/2.0;
    real u2=u[2];
    int steps_taken=0;
    while(steps_taken<max_steps && fabs(u1-u0)>rel_tol){
      real f0=rootfun(u0, x, gamma,  sigma);
      real f1=rootfun(u1, x, gamma,  sigma);
      real f2=rootfun(u2, x, gamma,  sigma);
      real l0=(u0*f1*f2)/((f0-f1)*(f0-f2));
      real l1=(u1*f0*f2)/((f1-f0)*(f1-f2));
      real l2=(u2*f1*f0)/((f2-f0)*(f2-f1));
      real n = l0+l1+l2;
      u0=n;
      u1=u0;
      u2=u1;
      steps_taken += 1;
    }
    if(verbose) print("Steps taken ", steps_taken);
    return u0;
  }

  real iqf_brent(row_vector u, real max_steps, real rel_tol, real x, real gamma, real sigma, int verbose){
    real u0=u[1];
    real u1=u[2];
    int steps_taken=0;
    real f0=rootfun(u0, x, gamma,  sigma);
    real f1=rootfun(u1, x, gamma,  sigma);
    real t=0;
    if(f0*f1>0) reject("Root is not bracketed!");

    if (fabs(f0)<fabs(f1)){
      t=u0;
      u0=u1;
      u1=t;
      t=f0;
      f0=f1;
      f1=t;
    }
    real u2=u0;
    real f2=f0;
    int mflag=1;
    real n=0;
    real d=0;

    while(steps_taken<max_steps && fabs(u1-u0)>rel_tol){
      f0=rootfun(u0, x, gamma,  sigma);
      f1=rootfun(u1, x, gamma,  sigma);
      f2=rootfun(u2, x, gamma,  sigma);

      if(f0!=f2 && f1 !=f2){
        real l0=(u0*f1*f2)/((f0-f1)*(f0-f2));
        real l1=(u1*f0*f2)/((f1-f0)*(f1-f2));
        real l2=(u2*f1*f0)/((f2-f0)*(f2-f1));
        n = l0+l1+l2;
      } else {
        n = u1-(f1*(u1-u0)/(f1-f0));
      }

      if((n<(3*u0+u1)/4 || n>u1) ||
      (mflag==1 && fabs(n-u1)>=fabs(u1-u2)/2.0) ||
      (mflag==0 && fabs(n-u1)>=fabs(u2-d)/2.0) ||
      (mflag==1 && fabs(u1-u2)<rel_tol ) ||
      (mflag==0 && fabs(u2-d)<rel_tol)
      ){
        n=(u0+u1)/2.0;
        mflag=1;
      } else {
        mflag=0;
      }
      real fn=rootfun(n, x, gamma,  sigma);
      d=u2;
      u2=u1;
      if(f0*fn<0){
        u1=n;
      } else {
        u0=n;
      }
      if(fabs(f0)<fabs(f1)){
        t=u0;
        u0=u1;
        u1=t;
      }
      steps_taken+=1;
    }
   if(verbose) print("Steps taken ", steps_taken);
    return u1;

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
  int verbose;
}
transformed data{
  vector[N] x_srt = sort_asc(x);
}
parameters {
  real<lower=1e-9> gamma; // gamma for govindarajulu
  real<lower=1e-9> dsigma; // exponential "tail" of sigma
}
model {
  vector[N] u;
  real sigma=min_sigma+dsigma;// sigma  for govindarajulu
  // create grid of xs given the parameter
  vector[M] xs_grd = govindarajulu_v_qf(ys_grd, gamma, sigma);
  matrix[N,2] bounds = lookup_v_bracket(x_srt, xs_grd, ys_grd); // only add this line
  target += genexp_s_lpdf(gamma | genexp_alpha, genexp_lambda);
  target += exponential_lpdf(dsigma | exp_lambda);
  for (i in 1:N){
    // only add these lines from here
     row_vector[2] pair=bounds[i,]; //[1e-15, 1-1e-15];
    //u[i] = iqf_bisection(pair, max_steps, rel_tol, x_srt[i], gamma, sigma, verbose);
    //u[i] =  iqf_newton_raphson(pair, max_steps, rel_tol, x_srt[i], gamma, sigma, verbose);
    //u[i] = iqf_secant(pair, max_steps, rel_tol, x_srt[i], gamma, sigma, verbose);
    //u[i] = iqf_lagrangepoly(pair, max_steps, rel_tol, x_srt[i], gamma, sigma, verbose);
    u[i] = iqf_brent(pair, max_steps, rel_tol, x_srt[i], gamma, sigma, verbose);
    // only add these lines down to here
   target += govindarajulu_s_ldqf_lpdf(u[i] | gamma, sigma);
  }
}
// end stan/gexpgovi_pdf_dqf_mod.stan
