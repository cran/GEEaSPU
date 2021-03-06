# This file was generated by Rcpp::compileAttributes
# Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

rowMinsC <- function(x) {
    .Call('GEEaSPU_rowMinsC', PACKAGE = 'GEEaSPU', x)
}

count_if <- function(x) {
    .Call('GEEaSPU_count_if', PACKAGE = 'GEEaSPU', x)
}

rankC <- function(x) {
    .Call('GEEaSPU_rankC', PACKAGE = 'GEEaSPU', x)
}

which_min <- function(y, x) {
    .Call('GEEaSPU_which_min', PACKAGE = 'GEEaSPU', y, x)
}

getEigen <- function(sigma) {
    .Call('GEEaSPU_getEigen', PACKAGE = 'GEEaSPU', sigma)
}

a8 <- function(U, p, k) {
    .Call('GEEaSPU_a8', PACKAGE = 'GEEaSPU', U, p, k)
}

signC <- function(x) {
    .Call('GEEaSPU_signC', PACKAGE = 'GEEaSPU', x)
}

rcpp_standPow <- function(U1, power) {
    .Call('GEEaSPU_rcpp_standPow', PACKAGE = 'GEEaSPU', U1, power)
}

InfU <- function(U) {
    .Call('GEEaSPU_InfU', PACKAGE = 'GEEaSPU', U)
}

gauss_score <- function(invR, G, res, n, k, p) {
    .Call('GEEaSPU_gauss_score', PACKAGE = 'GEEaSPU', invR, G, res, n, k, p)
}

gauss_score_cov <- function(invR, G, x, res, covres, n, k, p, tp) {
    .Call('GEEaSPU_gauss_score_cov', PACKAGE = 'GEEaSPU', invR, G, x, res, covres, n, k, p, tp)
}

bin_score <- function(va, invR, G, x, res, covres, n, k, p, tp) {
    .Call('GEEaSPU_bin_score', PACKAGE = 'GEEaSPU', va, invR, G, x, res, covres, n, k, p, tp)
}

bin_score_cov <- function(va, invR, G, x, res, covres, n, k, p, tp) {
    .Call('GEEaSPU_bin_score_cov', PACKAGE = 'GEEaSPU', va, invR, G, x, res, covres, n, k, p, tp)
}

rcpp_spuval <- function(U, po, po2, p, k) {
    .Call('GEEaSPU_rcpp_spuval', PACKAGE = 'GEEaSPU', U, po, po2, p, k)
}

perm <- function(invR, G, res, n, k, p, po, po2, nperm) {
    .Call('GEEaSPU_perm', PACKAGE = 'GEEaSPU', invR, G, res, n, k, p, po, po2, nperm)
}

perm_score <- function(invR, G, x, res, covres, n, k, p, tp, po, po2, nperm) {
    .Call('GEEaSPU_perm_score', PACKAGE = 'GEEaSPU', invR, G, x, res, covres, n, k, p, tp, po, po2, nperm)
}

permhigh <- function(invR, G, res, n, k, p, po, po2, nperm) {
    .Call('GEEaSPU_permhigh', PACKAGE = 'GEEaSPU', invR, G, res, n, k, p, po, po2, nperm)
}

permhigh_score <- function(invR, G, x, res, covres, n, k, p, tp, po, po2, nperm) {
    .Call('GEEaSPU_permhigh_score', PACKAGE = 'GEEaSPU', invR, G, x, res, covres, n, k, p, tp, po, po2, nperm)
}

sim <- function(f, va, invR, G, x, res, covres, n, k, p, tp, po, po2, nsim) {
    .Call('GEEaSPU_sim', PACKAGE = 'GEEaSPU', f, va, invR, G, x, res, covres, n, k, p, tp, po, po2, nsim)
}

sim_score <- function(f, va, invR, G, x, res, covres, n, k, p, tp, po, po2, nsim) {
    .Call('GEEaSPU_sim_score', PACKAGE = 'GEEaSPU', f, va, invR, G, x, res, covres, n, k, p, tp, po, po2, nsim)
}

simhigh <- function(f, va, invR, G, x, res, covres, n, k, p, tp, po, po2, nperm) {
    .Call('GEEaSPU_simhigh', PACKAGE = 'GEEaSPU', f, va, invR, G, x, res, covres, n, k, p, tp, po, po2, nperm)
}

simhigh_score <- function(f, va, invR, G, x, res, covres, n, k, p, tp, po, po2, nperm) {
    .Call('GEEaSPU_simhigh_score', PACKAGE = 'GEEaSPU', f, va, invR, G, x, res, covres, n, k, p, tp, po, po2, nperm)
}

pathval <- function(out, k, nSNPs, nGenes, po, po2, po3) {
    .Call('GEEaSPU_pathval', PACKAGE = 'GEEaSPU', out, k, nSNPs, nGenes, po, po2, po3)
}

permpath <- function(invR, G, res, nSNPs, nGenes, n, k, p, po, po2, po3, nperm) {
    .Call('GEEaSPU_permpath', PACKAGE = 'GEEaSPU', invR, G, res, nSNPs, nGenes, n, k, p, po, po2, po3, nperm)
}

