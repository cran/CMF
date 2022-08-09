#include <cpp11/doubles.hpp>
#include <cpp11/integers.hpp>
#include <cpp11/list.hpp>
#include <cpp11/matrix.hpp>
using namespace cpp11;

[[cpp11::register]] void p_gradUsparse(const doubles_matrix<> Xm, doubles_matrix<> Gm,
                                       const doubles_matrix<> CUm,
                                       const doubles_matrix<> OUm,
                                       const doubles_matrix<> Cm, const int idx,
                                       const double tau, const doubles Rowm,
                                       const doubles Colm) {
  double* const pGm = REAL(Gm.data());

  const int N = Xm.nrow();
  const int K = Gm.ncol();

  // If idx == 1 => update a row entity
  // If idx == 2 => update a column entity
  if (idx == 1) {
    for (int n = 0; n < N; n++) {
      const int r = Xm(n, 0) - 1;
      const int c = Xm(n, 1) - 1;

      double tmp = 0.0;
      for (int k = 0; k < K; k++) {
        tmp += CUm(r, k) * OUm(c, k);
      }
      tmp += -Xm(n, 2) + Rowm[r] + Colm[c];

      for (int k = 0; k < K; k++) {
        pGm[r + k * K] += tau * (tmp * OUm(c, k) + CUm(r, k) * Cm(c, k));
      }
    }
  } else {
    for (int n = 0; n < N; n++) {
      const int r = Xm(n, 0) - 1;
      const int c = Xm(n, 1) - 1;

      double tmp = 0.0;
      for (int k = 0; k < K; k++) {
        tmp += CUm(c, k) * OUm(r, k);
      }
      tmp += -Xm(n, 2) + Rowm[r] + Colm[c];

      for (int k = 0; k < K; k++) {
        pGm[c + k * K] += tau * (tmp * OUm(r, k) + CUm(c, k) * Cm(r, k));
      }
    }
  }
}

[[cpp11::register]] doubles p_updatePseudoData(const integers_matrix<> indices,
                                               const doubles_matrix<> U1m,
                                               const doubles_matrix<> U2m,
                                               const doubles Rv, const doubles Cv) {
  const int N = indices.nrow();
  const int K = U1m.ncol();

  writable::doubles out(N);

  for (int n = 0; n < N; n++) {
    const int r = indices(n, 0) - 1;
    const int c = indices(n, 1) - 1;

    double tmp = 0.0;
    for (int k = 0; k < K; k++) {
      tmp += U1m(r, k) * U2m(c, k);
    }
    out[n] = tmp + Rv[r] + Cv[c];
  }

  return out;
}

[[cpp11::register]] double p_updateTau(
    const doubles_matrix<> Xm, const doubles_matrix<> U1m, const doubles_matrix<> U2m,
    const doubles_matrix<> cov1m, const doubles_matrix<> cov2m, const doubles Rv,
    const doubles Cv, const doubles nu1v, const doubles nu2v) {
  const int N = Xm.nrow();
  const int K = U1m.ncol();
  double out = 0.0;

  for (int n = 0; n < N; n++) {
    const int r = static_cast<int>(Xm(n, 0)) - 1;
    const int c = static_cast<int>(Xm(n, 1)) - 1;

    double tmp = 0.0;
    for (int k = 0; k < K; k++) {
      tmp += U1m(r, k) * U2m(c, k);
    }
    tmp += Rv[r] + Cv[c];
    tmp = Xm(n, 2) - tmp;
    tmp = tmp * tmp;
    for (int k = 0; k < K; k++) {
      tmp += cov1m(r, k) * U2m(c, k) * U2m(c, k) + U1m(r, k) * U1m(r, k) * cov2m(c, k) +
             cov1m(r, k) * cov2m(c, k);
    }
    tmp += nu1v[r] + nu2v[c];

    out += tmp;
  }

  return out;
}

[[cpp11::register]] sexp p_updateMean(const doubles_matrix<> Xm,
                                      const doubles_matrix<> U1m,
                                      const doubles_matrix<> U2m, const int idx,
                                      const doubles Mv) {
  const int N = Xm.nrow();
  const int K = U1m.ncol();

  // If idx == 1 => update a row entity
  // If idx == 2 => update a column entity
  writable::doubles Nv(idx == 1 ? U1m.nrow() : U2m.nrow());
  writable::integers Cv(idx == 1 ? U1m.nrow() : U2m.nrow());
  for (int i = 0; i < Cv.size(); i++) {
    Nv[i] = 0.0;
    Cv[i] = 0;
  }

  for (int n = 0; n < N; n++) {
    const int r = Xm(n, 0) - 1;
    const int c = Xm(n, 1) - 1;

    double tmp = 0.0;
    for (int k = 0; k < K; k++) {
      tmp += U1m(r, k) * U2m(c, k);
    }
    if (idx == 1) {
      tmp = Xm(n, 2) - tmp - Mv[c];
      Nv[r] += tmp;
      Cv[r]++;
    } else {
      tmp = Xm(n, 2) - tmp - Mv[r];
      Nv[c] += tmp;
      Cv[c]++;
    }
  }

  return writable::list({"sum"_nm = Nv, "count"_nm = Cv});
}

[[cpp11::register]] void p_covUsparse(const doubles_matrix<> Xm, doubles_matrix<> Cm,
                                      const doubles_matrix<> OUm,
                                      const doubles_matrix<> OCm, const int idx,
                                      const double tau) {
  double* const pCm = REAL(Cm.data());

  const int N = Xm.nrow();
  const int K = Cm.ncol();

  if (idx == 1) {
    for (int n = 0; n < N; n++) {
      const int r = Xm(n, 0) - 1;
      const int c = Xm(n, 1) - 1;

      for (int k = 0; k < K; k++) {
        pCm[r + k * K] += tau * (OUm(c, k) * OUm(c, k) + OCm(c, k));
      }
    }
  } else {
    for (int n = 0; n < N; n++) {
      const int r = Xm(n, 0) - 1;
      const int c = Xm(n, 1) - 1;

      for (int k = 0; k < K; k++) {
        pCm[c + k * K] += tau * (OUm(r, k) * OUm(r, k) + OCm(r, k));
      }
    }
  }
}
