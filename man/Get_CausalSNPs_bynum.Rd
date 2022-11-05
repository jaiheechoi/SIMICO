\name{Get_CausalSNPs_bynum}
\alias{Get_CausalSNPs_bynum}
\title{
Get_CausalSNPs_bynum()
}
\description{Matrix of subsetted genetic information.
}
\usage{
Get_CausalSNPs_bynum(gMat, num, Causal.MAF.Cutoff)
}
\arguments{
  \item{gMat}{Matrix of SNPs.}
  \item{num}{Number of causal variants.}
  \item{Causal.MAF.Cutoff}{Minor allele frequency value cutoff for causal SNPs.}
}
\value{
Output is a vector of indices to subset the full genetic matrix.
}
