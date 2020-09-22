\name{pmd}
\alias{dpmn}
\alias{pmatrix}
\alias{dpm_sim}
\alias{dpm_noraml}
\title{
The Poisson Multinomial Distribution.
}
\description{
The probability mass function of Poisson Multinomial Distribution and a Poisson Multinomial probability matrix generator.
}
\usage{
dpmn(kk=NULL,pp)
pmatrix(n,m)
dpm_sim(kk=NULL,pp,t)
dpm_noraml(kk=NULL,pp,x_vec)
}
\arguments{
  \item{kk}{
   Null
}
  \item{pp}{ The matrix for \eqn{p_j}'s which are the sucess probabilities matrix for indicators.
}
  \item{t}{ Simulation repetitions.
}
  \item{n}{ Number of rows of the probabilities matrix.
}
  \item{m}{ Number of columns of the probabilities matrix.
}
  \item{x_vec}{ A result vector that is used for calculating the probability.

  }
}
\details{
See the reference for computational details.
}
\value{

}
\references{
Hong, Y. (2013). On computing the distribution function for the Poisson binomial distribution. Computational Statistics & Data Analysis, Vol. 59, pp. 41-51.
}
\author{

\packageAuthor{poissonmulti}

}


\examples{
   aa=matrix(c(.2, .1, .7, .2, .4, .4, .2, .1, .7, .2, .4, .4, .2, .1, .7, .2, .4, .4), nrow=6, byrow=TRUE)
   dpmn(kk=NULL,pp=aa[1:2,])
   n <- 4
   m <- 3
   pp <- pmatrix(n,m)
   dpm_sim(kk=NULL, pp,1000)
}


