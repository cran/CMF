# Generated by cpp11: do not edit by hand

p_gradUsparse <- function(Xm, Gm, CUm, OUm, Cm, idx, tau, Rowm, Colm) {
  invisible(.Call(`_CMF_p_gradUsparse`, Xm, Gm, CUm, OUm, Cm, idx, tau, Rowm, Colm))
}

p_updatePseudoData <- function(indices, U1m, U2m, Rv, Cv) {
  .Call(`_CMF_p_updatePseudoData`, indices, U1m, U2m, Rv, Cv)
}

p_updateTau <- function(Xm, U1m, U2m, cov1m, cov2m, Rv, Cv, nu1v, nu2v) {
  .Call(`_CMF_p_updateTau`, Xm, U1m, U2m, cov1m, cov2m, Rv, Cv, nu1v, nu2v)
}

p_updateMean <- function(Xm, U1m, U2m, idx, Mv) {
  .Call(`_CMF_p_updateMean`, Xm, U1m, U2m, idx, Mv)
}

p_covUsparse <- function(Xm, Cm, OUm, OCm, idx, tau) {
  invisible(.Call(`_CMF_p_covUsparse`, Xm, Cm, OUm, OCm, idx, tau))
}
