// Beta-binomial mixed-effects model (BBmm).
// y | u ~ BB(m, p, phi), logit(p) = X beta + Z u, u ~ N(0, D)
// D is block-diagonal: component c has variance sigma_c^2.
#include <TMB.hpp>

template <class Type>
Type dbetabinom_log(Type y, Type m, Type p, Type phi) {
  Type a = p / phi;
  Type b = (Type(1) - p) / phi;
  Type out = lgamma(m + Type(1)) - lgamma(y + Type(1)) - lgamma(m - y + Type(1));
  out += lgamma(a + y) - lgamma(a);
  out += lgamma(b + m - y) - lgamma(b);
  out += lgamma(a + b) - lgamma(a + b + m);
  return out;
}

template <class Type>
Type objective_function<Type>::operator()() {
  DATA_VECTOR(y);
  DATA_VECTOR(m);
  DATA_MATRIX(X);
  DATA_MATRIX(Z);
  // re_comp(j) in 0..nComp-1: which variance component owns u(j)
  DATA_IVECTOR(re_comp);
  // Optional multidimensional: dim_id(i) in 0..nDim-1
  DATA_IVECTOR(dim_id);
  DATA_INTEGER(nDim);

  PARAMETER_VECTOR(beta);
  PARAMETER_VECTOR(log_phi);   // length nDim
  PARAMETER_VECTOR(log_sigma); // length nComp
  PARAMETER_VECTOR(u);         // random effects

  int n = y.size();
  int nRand = u.size();
  int nComp = log_sigma.size();

  vector<Type> sigma(nComp);
  for (int c = 0; c < nComp; c++) {
    sigma(c) = exp(log_sigma(c));
  }

  vector<Type> phi(nDim);
  for (int d = 0; d < nDim; d++) {
    phi(d) = exp(log_phi(d));
  }

  vector<Type> eta = X * beta + Z * u;
  Type nll = Type(0);

  Type eps = Type(1e-8);
  for (int i = 0; i < n; i++) {
    Type p = Type(1) / (Type(1) + exp(-eta(i)));
    p = CppAD::CondExpLt(p, eps, eps, p);
    p = CppAD::CondExpGt(p, Type(1) - eps, Type(1) - eps, p);
    int d = dim_id(i);
    nll -= dbetabinom_log(y(i), m(i), p, phi(d));
  }

  for (int j = 0; j < nRand; j++) {
    int c = re_comp(j);
    nll -= dnorm(u(j), Type(0), sigma(c), true);
  }

  ADREPORT(phi);
  ADREPORT(sigma);
  REPORT(u);
  return nll;
}
