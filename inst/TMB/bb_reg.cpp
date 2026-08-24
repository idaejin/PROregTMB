// Beta-binomial logistic regression + Eilers mixed P-splines.
// logit(p) = X beta + Zs s
// s_j ~ N(0, sigma_s[comp]^2)  (Eilers 1999; null space in X)
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

  DATA_INTEGER(n_smooth);
  DATA_MATRIX(Zs);
  DATA_IVECTOR(s_comp);

  PARAMETER_VECTOR(beta);
  PARAMETER(log_phi);
  PARAMETER_VECTOR(log_sds);
  PARAMETER_VECTOR(s);

  Type phi = exp(log_phi);
  int n = y.size();
  Type nll = Type(0);

  vector<Type> eta = X * beta;
  if (n_smooth > 0) {
    eta += Zs * s;
  }

  Type eps = Type(1e-8);
  vector<Type> p(n);
  for (int i = 0; i < n; i++) {
    p(i) = Type(1) / (Type(1) + exp(-eta(i)));
    p(i) = CppAD::CondExpLt(p(i), eps, eps, p(i));
    p(i) = CppAD::CondExpGt(p(i), Type(1) - eps, Type(1) - eps, p(i));
    nll -= dbetabinom_log(y(i), m(i), p(i), phi);
  }

  if (n_smooth > 0) {
    vector<Type> sds = exp(log_sds);
    for (int j = 0; j < s.size(); j++) {
      int c = s_comp(j);
      nll -= dnorm(s(j), Type(0), sds(c), true);
    }
    ADREPORT(sds);
    REPORT(s);
  }

  ADREPORT(phi);
  ADREPORT(p);
  return nll;
}
