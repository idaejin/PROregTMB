// Beta-binomial logistic regression (optional additive P-splines).
// logit(p) = X beta + B alpha
// alpha_j ~ GMRF(Q_j), Q_j = lambda_j S_j + kappa C_j
//   S_j = D'D (Eilers), C_j = (B_j'1)(1'B_j) sum-to-zero (additive notes)
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
  DATA_MATRIX(B);
  DATA_MATRIX(S);
  DATA_MATRIX(C);
  DATA_IVECTOR(smooth_K);
  DATA_IVECTOR(smooth_off);
  DATA_SCALAR(kappa);

  PARAMETER_VECTOR(beta);
  PARAMETER(log_phi);
  PARAMETER_VECTOR(log_lambda);
  PARAMETER_VECTOR(alpha);

  Type phi = exp(log_phi);
  int n = y.size();
  Type nll = Type(0);

  vector<Type> eta = X * beta;
  if (n_smooth > 0) {
    eta += B * alpha;
  }

  Type eps = Type(1e-8);
  vector<Type> p(n);
  for (int i = 0; i < n; i++) {
    p(i) = Type(1) / (Type(1) + exp(-eta(i)));
    p(i) = CppAD::CondExpLt(p(i), eps, eps, p(i));
    p(i) = CppAD::CondExpGt(p(i), Type(1) - eps, Type(1) - eps, p(i));
    nll -= dbetabinom_log(y(i), m(i), p(i), phi);
  }

  // Additive P-splines: difference GMRF + soft sum-to-zero
  // Q = lambda S + eps I;  nll += 0.5 * kappa * (1' B_j gamma_j)^2
  // (kappa term matches pspline_additive.R without stuffing C into Q)
  if (n_smooth > 0) {
    using namespace density;
    Type ridge = Type(1e-6);
    for (int s = 0; s < n_smooth; s++) {
      Type lam = exp(log_lambda(s));
      int K = smooth_K(s);
      int off = smooth_off(s);
      matrix<Type> Q(K, K);
      Q.setZero();
      for (int i = 0; i < K; i++) {
        for (int j = 0; j < K; j++) {
          Q(i, j) = lam * S(off + i, off + j);
        }
        Q(i, i) += ridge;
      }
      vector<Type> as(K);
      for (int i = 0; i < K; i++) as(i) = alpha(off + i);
      nll += GMRF(asSparseMatrix(Q))(as);

      // sum_i f_j(x_ij) = 1' B_j gamma_j
      Type sumf = Type(0);
      for (int i = 0; i < n; i++) {
        Type fi = Type(0);
        for (int k = 0; k < K; k++) fi += B(i, off + k) * as(k);
        sumf += fi;
      }
      nll += Type(0.5) * kappa * sumf * sumf;
    }
    vector<Type> lambda(n_smooth);
    for (int s = 0; s < n_smooth; s++) lambda(s) = exp(log_lambda(s));
    ADREPORT(lambda);
    REPORT(alpha);
  }

  ADREPORT(phi);
  ADREPORT(p);
  return nll;
}
