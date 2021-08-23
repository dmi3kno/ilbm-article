functions{
  // reimplemented from https://nickcdryan.com/2017/09/13/root-finding-algorithms-in-python-line-search-bisection-secant-newton-raphson-boydens-inverse-quadratic-interpolation-brents/
  // relevant https://math.stackexchange.com/questions/1464795/what-are-the-difference-between-some-basic-numerical-root-finding-methods
  // fortran implementation https://www.cantorsparadise.com/some-root-finding-algorithms-5c6fa8a4a165
  //  //given x and a grid, returns a linearly interpolated y
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

  real exponential_s_ldqf_lpdf(real p, real lambda){
   if (p>=1) reject("p>=1, found p=",p);
   return log(lambda)+log1m(p);
  }

  real exponential_s_qf(real u, real lambda){
    return -log1m(u)/lambda;
  }

  vector exponential_v_qf(vector p, real lambda){
   return -log1m(p)/lambda;
  }

  real drootfun_du(real u, real lambda){
    return inv(lambda*(1-u));
  }

  real rootfun(real u, real x, real lambda){
    return x-exponential_s_qf(u, lambda);
  }

  real iqf_bisection(row_vector u, real max_iter, real tol, real x, real lambda, int verbose){
    int steps_taken=0;
    real i_lower=u[1];
    real i_upper=u[2];
    while (steps_taken < max_iter){
      real m = (i_lower+i_upper) / 2.0;
      if (m==0 || fabs(i_upper-i_lower)<tol){
        if(verbose) print("Steps taken ", steps_taken);
        return m;
      }
      if (rootfun(m, x, lambda)>0){ //rootfinding function is called once
        i_lower=m;
      } else {
        i_upper=m;
      }
      steps_taken+=1;
    }
    if(verbose) print("Steps taken ", steps_taken);
    return (i_lower+i_upper)/2.0;
  }

  real iqf_newton_raphson(row_vector u, real max_iter, real tol, real x, real lambda, int verbose){
    int steps_taken=0;
    real u_i = u[1];
    while(steps_taken<max_iter && fabs(rootfun(u_i, x, lambda))>tol){
      real drf=drootfun_du(u_i, lambda);
      u_i+=rootfun(u_i, x, lambda)/drf;
      steps_taken +=1;
    }
    if(verbose) print("Steps taken ", steps_taken);
    return u_i;
  }

  real iqf_secant(row_vector u, real max_iter, real tol, real x, real lambda, int verbose){
    int steps_taken=1;
    real u0=u[1];
    real u1=u[2];
    real u2=0;
    while (steps_taken<max_iter && fabs(u1-u0)>tol){
      real f1=rootfun(u1, x, lambda);
      real f0=rootfun(u0, x, lambda);
      u2 = u1-(f1*(u1-u0)/(f1-f0));
      u1=u2;
      u0=u1;
      steps_taken+=1;
    }
    if(verbose) print("Steps taken ", steps_taken);
    return u2;
  }

  real iqf_lagrangepoly(row_vector u, real max_iter, real tol, real x, real lambda, int verbose){
    real u0=u[1];
    real u1=(u[1]+u[2])/2.0;
    real u2=u[2];
    int steps_taken=0;
    while(steps_taken<max_iter && fabs(u1-u0)>tol){
      real f0=rootfun(u0, x, lambda);
      real f1=rootfun(u1, x, lambda);
      real f2=rootfun(u2, x, lambda);
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

  real iqf_brent(row_vector u, real max_iter, real tol, real x, real lambda, int verbose){
    real u0=u[1];
    real u1=u[2];
    int steps_taken=0;
    real f0=rootfun(u0, x, lambda);
    real f1=rootfun(u1, x, lambda);
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

    while(steps_taken<max_iter && fabs(u1-u0)>tol){
      f0=rootfun(u0, x, lambda);
      f1=rootfun(u1, x, lambda);
      f2=rootfun(u2, x, lambda);

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
      (mflag==1 && fabs(u1-u2)<tol ) ||
      (mflag==0 && fabs(u2-d)<tol)
      ){
        n=(u0+u1)/2.0;
        mflag=1;
      } else {
        mflag=0;
      }
      real fn=rootfun(n, x, lambda);
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

} // end of functions block

data {
  int<lower=0> N;
  vector[N] x;
  int<lower=0> M;
  vector[M] ys_grd;
  real gamma_a;
  real gamma_b;
  real  tol;
  real max_iter;
  int verbose;
}
// The input data is a vector 'y' of length 'N'.

transformed data{
  vector[N] x_srt = sort_asc(x);
}

// The parameters accepted by the model. Our model
// accepts two parameters 'mu' and 'sigma'.
parameters {
  real<lower=1e-15> lambda;
}

model{
  vector[N] u;
  vector[M] xs_grd = exponential_v_qf(ys_grd, lambda);
  //print(xs_grd);
  matrix[N,2] bounds = lookup_v_bracket(x_srt, xs_grd, ys_grd);
  lambda ~ gamma(gamma_a, gamma_b);

   for (i in 1:N){
     row_vector[2] pair=bounds[i,]; //[1e-15, 1-1e-15];
    u[i] = iqf_bisection(pair, max_iter, tol, x_srt[i], lambda, verbose);
    //u[i] =  iqf_newton_raphson(pair, max_iter, tol, x_srt[i], lambda, verbose);
    // u[i] = iqf_secant(pair, max_iter, tol, x_srt[i], lambda, verbose);
   // u[i] = iqf_lagrangepoly(pair, max_iter, tol, x_srt[i], lambda, verbose);
    //u[i] = iqf_brent(pair, max_iter, tol, x_srt[i], lambda, verbose);
    u[i] ~ exponential_s_ldqf(lambda);
  }
}

