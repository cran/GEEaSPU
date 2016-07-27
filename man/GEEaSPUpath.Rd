\name{GEEaSPUpath}
\alias{GEEaSPUpath}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
An adaptive sum of powered score (SPU) test for gene-set (pathway) and mutiple phenotypes
}
\description{
Tests whether gene-set (pathway) is associated with multivariate outcomes;
provides a series of p-values of GEE-SPU(pow, pow2, pow3) and GEEaSPUpath tests.
}
\usage{
GEEaSPUpath(pheno, geno, nSNPs, Z = NULL, corstr = "independence", 
           pow = c(1,2,4,8), pow2 = c(1,2,4,8), pow3 = c(1,2,4,8), n.perm = 1000)}

\arguments{
\item{pheno}{a numeric phenotype matrix with each row as a different individual and each column as a separate phenotype.}
\item{geno}{a numeric genotype matrix with each row as a different individual and each column as a snp; the SNPs (with the number stored in nSNPs) from one gene are 
stored consecutively.}
\item{nSNPs}{A numeric vector, whose length matches to the total number of genes; 
each element of vector indicate the number of SNPs in each gene.}
\item{Z}{a numeric covariate matrix with each row as a different individual and each column as a covariated to be adjested.}
\item{corstr}{a character string specifying the correlation structure of phenotypes. The following are permitted: "independence", "fixed", "stat_M_dep", "non_stat_M_dep", "exchangeable", "AR-M" and "unstructured" (default = "independence").
}
\item{pow}{a vector of the power weight to be used on each SNP (default = c(1,2,4,8)).}
\item{pow2}{a vector of the power weight to be used on each phenotype (default = c(1,2,4,8)).}
\item{pow3}{a vector of the power weight to be used on each gene (default = c(1,2,4,8)).}
\item{n.perm}{a numeric value of number of null statistics (default = 1000).}
}
\details{
Adaptive association tests for gene-set (pathway) and multiple phenotypes using GEE.
}
\value{a vector of p-values from GEE-SPU(pow, pow2, pow3) tests and GEE-aSPUpath test.
}
\references{
Kim J, Zhang Y, and Wei P (2016). "Powerful and adaptive testing for multi-trait and multi-SNP associations with GWAS and sequencing data." \emph{Genetics}, 203(2), 715--731.
}
\author{Junghi Kim and Wei Pan
}
\note{
GEEaSPUpath function only supports a case for a quantitative trait (model = "gaussian") and a permutation based test (null.type = "perm"). 
}


\examples{
# -- simulating phenotypes
# -- n.subjects: number of subjects
# -- n.traits: number of phenotypes
# -- Sigma: covariance matrix of phenotypes (e.g. AR(1))

set.seed(136)
n.subjects <- 100
n.traits <- 3
sigma <- 2; rho <- 0.5
Sigma0 <- diag(n.traits)
Sigma <- sigma * rho^abs(row(Sigma0) - col(Sigma0))
eS <- eigen(Sigma, symmetric = TRUE)
ev <- eS$values
X <- matrix(rnorm(n.subjects * n.traits), n.subjects)
pheno <- X \%*\% diag(sqrt(pmax(ev, 0)), ncol(Sigma)) \%*\% eS$vectors


# -- simulating genotype
# -- Assume we have two genes each of which has 3 and 5 SNPs respectively.
# -- n.geno1: number of SNPs included in the gene1
# -- n.geno2: number of SNPs included in the gene1
# -- nSNPs <- c(3,5)

n.geno1 <- 3
n.geno2 <- 5
maf1 <- 0.2
maf2 <- 0.4
gene1 <- matrix(rbinom(n = n.subjects*n.geno1, size = 2, prob = maf1), ncol = n.geno1)
gene2 <- matrix(rbinom(n = n.subjects*n.geno2, size = 2, prob = maf2), ncol = n.geno2)
geno <- cbind(gene1, gene2)


# -- Computing the p-value of GEEaSPUpath test
Pvl <- GEEaSPUpath(pheno = pheno, geno = geno, nSNPs = c(3,5), Z = NULL, 
        corstr = "independence", pow = c(1,4,8), pow2 = c(1,4,8), 
        pow3 = c(1,4,8), n.perm = 1000)

# -- P values of GEE-SPU(pow,pow2,pow3) and GEE-aSPUpath
Pvl

# -- P value of GEE-aSPUpath test
Pvl[length(Pvl)]
}
\keyword{ gene-set/pathway based}
\keyword{ multiple phenotypes}
