// Beta-binomial mixed-effects model (BBmm) + optional additive P-splines.
// logit(p) = X beta + Z u + B alpha
// Subject RE blocks as before; smooths: alpha_j ~ GMRF(lambda_j S_j + kappa C_j)
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
int n_theta_sigma(int q, int corr_type) {
  if (corr_type == 0) return q;
  return q * (q + 1) / 2;
}

template <class Type>
matrix<Type> sigma_from_theta(vector<Type> theta, int q, int corr_type,
                              vector<Type> &sd_out) {
  matrix<Type> Sigma(q, q);
  Sigma.setZero();
  if (corr_type == 0) {
    for (int t = 0; t < q; t++) {
      Type s = exp(theta(t));
      sd_out(t) = s;
      Sigma(t, t) = s * s;
    }
  } else {
    matrix<Type> L(q, q);
    L.setZero();
    int k = 0;
    for (int i = 0; i < q; i++) {
      for (int j = 0; j <= i; j++) {
        if (i == j) L(i, j) = exp(theta(k));
        else L(i, j) = theta(k);
        k++;
      }
    }
    Sigma = L * L.transpose();
    for (int t = 0; t < q; t++) sd_out(t) = sqrt(Sigma(t, t));
  }
  return Sigma;
}

template <class Type>
Type objective_function<Type>::operator()() {
  DATA_VECTOR(y);
  DATA_VECTOR(m);
  DATA_MATRIX(X);
  DATA_MATRIX(Z);
  DATA_IVECTOR(dim_id);
  DATA_INTEGER(nDim);

  DATA_INTEGER(n_blocks);
  DATA_IVECTOR(block_G);
  DATA_IVECTOR(block_q);
  DATA_IVECTOR(block_corr);
  DATA_IVECTOR(block_theta0);

  DATA_INTEGER(n_smooth);
  DATA_MATRIX(B);
  DATA_MATRIX(S);
  DATA_MATRIX(C);
  DATA_IVECTOR(smooth_K);
  DATA_IVECTOR(smooth_off);
  DATA_SCALAR(kappa);

  PARAMETER_VECTOR(beta);
  PARAMETER_VECTOR(log_phi);
  PARAMETER_VECTOR(theta_re);
  PARAMETER_VECTOR(u);
  PARAMETER_VECTOR(log_lambda);
  PARAMETER_VECTOR(alpha);

  int n = y.size();
  Type nll = Type(0);

  vector<Type> phi(nDim);
  for (int d = 0; d < nDim; d++) phi(d) = exp(log_phi(d));

  vector<Type> eta = X * beta + Z * u;
  if (n_smooth > 0) {
    eta += B * alpha;
  }

  Type eps = Type(1e-8);
  for (int i = 0; i < n; i++) {
    Type p = Type(1) / (Type(1) + exp(-eta(i)));
    p = CppAD::CondExpLt(p, eps, eps, p);
    p = CppAD::CondExpGt(p, Type(1) - eps, Type(1) - eps, p);
    int d = dim_id(i);
    nll -= dbetabinom_log(y(i), m(i), p, phi(d));
  }

  // Subject-level RE
  int u_offset = 0;
  for (int b = 0; b < n_blocks; b++) {
    int G = block_G(b);
    int q = block_q(b);
    int corr = block_corr(b);
    int nt = n_theta_sigma<Type>(q, corr);
    vector<Type> th(nt);
    for (int k = 0; k < nt; k++) th(k) = theta_re(block_theta0(b) + k);

    vector<Type> sd(q);
    matrix<Type> Sigma = sigma_from_theta(th, q, corr, sd);

    if (corr == 0) {
      for (int g = 0; g < G; g++) {
        for (int t = 0; t < q; t++) {
          nll -= dnorm(u(u_offset + g * q + t), Type(0), sd(t), true);
        }
      }
    } else {
      density::MVNORM_t<Type> mvn(Sigma);
      for (int g = 0; g < G; g++) {
        vector<Type> ug(q);
        for (int t = 0; t < q; t++) ug(t) = u(u_offset + g * q + t);
        nll += mvn(ug);
      }
    }

    if (b == 0) {
      ADREPORT(sd);
      if (q == 2 && corr == 1) {
        Type rho = Sigma(0, 1) / (sd(0) * sd(1) + Type(1e-12));
        ADREPORT(rho);
      }
    }

    u_offset += G * q;
  }

  // Additive P-splines
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
  REPORT(u);
  return nll;
}
