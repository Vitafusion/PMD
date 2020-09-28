\name{poissonmultinomial}
\alias{pmd}
\alias{pmatrix}
\title{
The Poisson Multinomial Distribution.
}
\description{
The probability mass function of Poisson Multinomial Distribution and a Poisson Multinomial probability matrix generator.
}
\usage{
pmd(pp,method="DFT-CF",vec=c(0,0,0,0,0),t=100)
pmatrix(n,m)
}
\arguments{
  \item{method}{The chosen method to compute probability mass.
}
  \item{pp}{ The matrix for \eqn{p_j}'s which are the sucess probabilities matrix for indicators.
}
  \item{t}{ Simulation repetitions.
}
  \item{n}{ Number of rows of the probabilities matrix.
}
  \item{m}{ Number of columns of the probabilities matrix.
}
  \item{vec}{ A result vector that is used for calculating the probability.

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

   n <- 4
   m <- 3
   pp <- pmatrix(n,m)
   pmd(pp,method="simulation",t=1000)
   vec <- c(1,2,1)
   pmd(pp,method="NA",vec,t=1000)
   pmd(pp)

}


