// One-stage joint model: Beta-Binomial mixed + Weibull PH (current-value).
// Longitudinal: y_ij | u_i ~ BB(m, p_ij, phi)
//   logit(p_ij) = (beta0 + u0_i) + (beta1 + u1_i) * t_ij
// Survival: h_i(t) = w1*w2*t^(w2-1) * exp(alpha * m_i(t))
//   m_i(t) = m * invlogit((beta0+u0_i)+(beta1+u1_i)*t)
// Random effects u0_i, u1_i ~ independent N(0, sigma0^2), N(0, sigma1^2)
// Laplace approximation via TMB.
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
Type mi_t(Type t, Type eta0, Type eta1, Type m) {
  return m * invlogit(eta0 + eta1 * t);
}

// Midpoint quadrature for cumulative hazard on (0, T]
template <class Type>
Type cumhaz_weibull_tvc(Type T, Type eta0, Type eta1, Type alpha, Type m,
                        Type w1, Type w2, int n_quad) {
  Type ch = Type(0);
  if (T <= Type(0)) return ch;
  Type ds = T / Type(n_quad);
  for (int k = 0; k < n_quad; k++) {
    Type s = ds * (Type(k) + Type(0.5));
    // s > 0 always with midpoint rule
    Type haz = w1 * w2 * pow(s, w2 - Type(1)) * exp(alpha * mi_t(s, eta0, eta1, m));
    ch += haz * ds;
  }
  return ch;
}

template <class Type>
Type objective_function<Type>::operator()() {
  // Longitudinal (long format)
  DATA_VECTOR(y);
  DATA_VECTOR(time_long);
  DATA_VECTOR(m_long);       // usually constant m
  DATA_IVECTOR(subj_long);   // 0..n_subj-1

  // Survival (one row per subject, order 0..n_subj-1)
  DATA_VECTOR(surv_time);
  DATA_IVECTOR(status);      // 1=event, 0=censored
  DATA_SCALAR(m);
  DATA_INTEGER(n_quad);

  PARAMETER_VECTOR(beta);        // length 2
  PARAMETER(log_phi);
  PARAMETER_VECTOR(log_sigma);   // length 2: RI, RS
  PARAMETER(alpha);
  PARAMETER(log_w1);
  PARAMETER(log_w2);
  PARAMETER_VECTOR(u0);          // random intercepts
  PARAMETER_VECTOR(u1);          // random slopes

  Type phi = exp(log_phi);
  Type sigma0 = exp(log_sigma(0));
  Type sigma1 = exp(log_sigma(1));
  Type w1 = exp(log_w1);
  Type w2 = exp(log_w2);

  int n_obs = y.size();
  int n_subj = u0.size();
  Type nll = Type(0);
  Type eps = Type(1e-8);

  // --- Longitudinal BB ---
  for (int j = 0; j < n_obs; j++) {
    int i = subj_long(j);
    Type eta = beta(0) + u0(i) + (beta(1) + u1(i)) * time_long(j);
    Type p = invlogit(eta);
    p = CppAD::CondExpLt(p, eps, eps, p);
    p = CppAD::CondExpGt(p, Type(1) - eps, Type(1) - eps, p);
    nll -= dbetabinom_log(y(j), m_long(j), p, phi);
  }

  // --- Survival + random effects ---
  for (int i = 0; i < n_subj; i++) {
    Type eta0 = beta(0) + u0(i);
    Type eta1 = beta(1) + u1(i);
    Type Ti = surv_time(i);
    Type ch = cumhaz_weibull_tvc(Ti, eta0, eta1, alpha, m, w1, w2, n_quad);

    // Survival contribution
    nll += ch; // -log S = H
    if (status(i) == 1) {
      Type Ti_safe = CppAD::CondExpLt(Ti, eps, eps, Ti);
      Type haz = w1 * w2 * pow(Ti_safe, w2 - Type(1)) *
                 exp(alpha * mi_t(Ti, eta0, eta1, m));
      nll -= log(haz);
    }

    // Random effects
    nll -= dnorm(u0(i), Type(0), sigma0, true);
    nll -= dnorm(u1(i), Type(0), sigma1, true);
  }

  ADREPORT(phi);
  ADREPORT(sigma0);
  ADREPORT(sigma1);
  ADREPORT(alpha);
  ADREPORT(w1);
  ADREPORT(w2);
  REPORT(u0);
  REPORT(u1);
  return nll;
}
