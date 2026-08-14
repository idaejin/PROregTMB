// BBreg with weak priors (for tmbstan smoke test only)
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

  PARAMETER_VECTOR(beta);
  PARAMETER(log_phi);

  Type phi = exp(log_phi);
  int n = y.size();
  Type nll = 0;

  // Weak priors (proper posterior for MCMC)
  for (int j = 0; j < beta.size(); j++) {
    nll -= dnorm(beta(j), Type(0), Type(5), true);
  }
  nll -= dnorm(log_phi, Type(log(0.3)), Type(1.0), true);

  vector<Type> eta = X * beta;
  Type eps = Type(1e-8);
  for (int i = 0; i < n; i++) {
    Type p = Type(1) / (Type(1) + exp(-eta(i)));
    p = CppAD::CondExpLt(p, eps, eps, p);
    p = CppAD::CondExpGt(p, Type(1) - eps, Type(1) - eps, p);
    nll -= dbetabinom_log(y(i), m(i), p, phi);
  }

  ADREPORT(phi);
  return nll;
}
