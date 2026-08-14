#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
IntegerVector rBB_cpp(int k, NumericVector m, NumericVector p, double phi) {
  if (k <= 0) stop("k must be positive");
  if (phi <= 0.0) stop("phi must be positive");
  const int nm = m.size();
  const int np = p.size();
  if (nm != 1 && nm != k) stop("m must be scalar or length k");
  if (np != 1 && np != k) stop("p must be scalar or length k");

  IntegerVector out(k);
  for (int i = 0; i < k; i++) {
    double mi = (nm == 1) ? m[0] : m[i];
    double pi = (np == 1) ? p[0] : p[i];
    if (mi < 1.0) stop("m must be positive");
    if (pi < 0.0 || pi > 1.0) stop("p must be in [0, 1]");
    double a = pi / phi;
    double b = (1.0 - pi) / phi;
    double u = R::rbeta(a, b);
    out[i] = R::rbinom(static_cast<int>(mi + 0.5), u);
  }
  return out;
}

// Pointwise log-PMF of BB for vectors y, m, p (scalar or length-n phi)
// [[Rcpp::export]]
NumericVector ldbb_cpp(NumericVector y, NumericVector m, NumericVector p,
                       NumericVector phi) {
  const int n = y.size();
  if (m.size() != n && m.size() != 1) stop("m length mismatch");
  if (p.size() != n && p.size() != 1) stop("p length mismatch");
  if (phi.size() != n && phi.size() != 1) stop("phi length mismatch");

  NumericVector out(n);
  for (int i = 0; i < n; i++) {
    double yi = y[i];
    double mi = (m.size() == 1) ? m[0] : m[i];
    double pi = (p.size() == 1) ? p[0] : p[i];
    double ph = (phi.size() == 1) ? phi[0] : phi[i];
    if (ph <= 0.0) stop("phi must be positive");
    // soft-bound p
    if (pi < 1e-12) pi = 1e-12;
    if (pi > 1.0 - 1e-12) pi = 1.0 - 1e-12;
    double a = pi / ph;
    double b = (1.0 - pi) / ph;
    out[i] = R::lgammafn(mi + 1.0) - R::lgammafn(yi + 1.0) - R::lgammafn(mi - yi + 1.0)
           + R::lgammafn(a + yi) - R::lgammafn(a)
           + R::lgammafn(b + mi - yi) - R::lgammafn(b)
           + R::lgammafn(a + b) - R::lgammafn(a + b + mi);
  }
  return out;
}

// Pointwise SE for f = B alpha with dense cov V: se_i = sqrt(B_i V B_i')
// [[Rcpp::export]]
NumericVector smooth_se_cpp(NumericMatrix B, NumericMatrix V) {
  const int n = B.nrow();
  const int K = B.ncol();
  if (V.nrow() != K || V.ncol() != K) stop("V must be K x K");
  NumericVector se(n);
  for (int i = 0; i < n; i++) {
    // t = V * B_i'
    NumericVector t(K);
    for (int j = 0; j < K; j++) {
      double s = 0.0;
      for (int k = 0; k < K; k++) s += V(j, k) * B(i, k);
      t[j] = s;
    }
    double q = 0.0;
    for (int k = 0; k < K; k++) q += B(i, k) * t[k];
    se[i] = (q > 0.0) ? std::sqrt(q) : 0.0;
  }
  return se;
}

// Full PMF on 0..m (scalar p, phi)
// [[Rcpp::export]]
NumericVector dBB_cpp(int m, double p, double phi) {
  if (m <= 0) stop("m must be positive");
  if (p < 0.0 || p > 1.0) stop("p must be in [0, 1]");
  if (phi <= 0.0) stop("phi must be positive");
  NumericVector out(m + 1);
  double a = p / phi;
  double b = (1.0 - p) / phi;
  for (int t = 0; t <= m; t++) {
    double lp = R::lgammafn(m + 1.0) - R::lgammafn(t + 1.0) - R::lgammafn(m - t + 1.0)
              + R::lgammafn(a + t) - R::lgammafn(a)
              + R::lgammafn(b + m - t) - R::lgammafn(b)
              + R::lgammafn(a + b) - R::lgammafn(a + b + m);
    out[t] = std::exp(lp);
  }
  return out;
}
