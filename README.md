
# The tenets of indirect inference in Bayesian models

<!-- badges: start -->
<!-- badges: end -->

This paper extends the application of Bayesian inference to probability distributions
defined in terms of quantile function. We introduce the method of indirect likelihood
to be used in the Bayesian models with sampling distributions defined by the quantile
function. We provide examples and demonstrate the equivalence of this “quantile-
based” (indirect) likelihood to the conventional “density-defined” (direct) likelihood.
We consider practical aspects of the numerical inversion of quantile function by root-
finding required by the indirect likelihood method. In particular we consider a problem
of ensuring the validity of an arbitrary quantile function with the help of Chebyshev
polynomials and provide useful tips and implementation of these algorithms in Stan
and R. We also extend the same method to propose the definition of an “indirect prior”
and discuss the situations where it can be useful.

