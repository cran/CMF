# Generated by using Rcpp::compileAttributes() -> do not edit by hand
# Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

p_gradUsparse <- function(Xm, Gm, CUm, OUm, Cm, I, T, Rowm, Colm) {
    invisible(.Call(`_CMF_p_gradUsparse`, Xm, Gm, CUm, OUm, Cm, I, T, Rowm, Colm))
}

p_updatePseudoData <- function(Xm, U1m, U2m, Rv, Cv) {
    invisible(.Call(`_CMF_p_updatePseudoData`, Xm, U1m, U2m, Rv, Cv))
}

p_updateTau <- function(Xm, U1m, U2m, cov1m, cov2m, Rv, Cv, nu1v, nu2v) {
    .Call(`_CMF_p_updateTau`, Xm, U1m, U2m, cov1m, cov2m, Rv, Cv, nu1v, nu2v)
}

p_updateMean <- function(Xm, U1m, U2m, I, Mv) {
    .Call(`_CMF_p_updateMean`, Xm, U1m, U2m, I, Mv)
}

p_covUsparse <- function(Xm, Cm, OUm, OCm, I, T) {
    invisible(.Call(`_CMF_p_covUsparse`, Xm, Cm, OUm, OCm, I, T))
}

