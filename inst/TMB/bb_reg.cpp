// Beta-binomial logistic regression (no random effects).
// Marginal MLE via TMB automatic differentiation.
#include <TMB.hpp>

template <class Type>
Type dbetabinom_log(Type y, Type m, Type p, Type phi) {
  // y | u ~ Bin(m, u), u ~ Beta(p/phi, (1-p)/phi)
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

  vector<Type> eta = X * beta;
  vector<Type> p(n);
  Type nll = 0;

  for (int i = 0; i < n; i++) {
    p(i) = Type(1) / (Type(1) + exp(-eta(i)));
    // Soft bounds to keep a, b > 0 and avoid Inf in logit extremes
    Type eps = Type(1e-8);
    p(i) = CppAD::CondExpLt(p(i), eps, eps, p(i));
    p(i) = CppAD::CondExpGt(p(i), Type(1) - eps, Type(1) - eps, p(i));
    nll -= dbetabinom_log(y(i), m(i), p(i), phi);
  }

  ADREPORT(phi);
  ADREPORT(p);
  return nll;
}
