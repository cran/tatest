\name{rhov}
\alias{rhov}

\title{
Calculation of rho (\eqn{\rho}) value
}
\description{
The rho is used to define a gap and variation between two datasets, similar to fc(fold change).
}
\usage{
rhov(XA, XB)
}

\arguments{
  \item{XA}{
a numeric vector in condition A.
}
  \item{XB}{
a numeric vector in condition B.
}
}
\details{
The rho is defined as square root of products of \deqn{\psi} and \deqn{\zeta}, that is \deqn{\rho=\sqrt{\psi\zeta}} where \eqn{\psi} is used to define a overlap between two datasets and \eqn{\zeta} is used to define effect variation. \eqn{\rho > 1} indicates that there is gap between two datasets and noise is smaller than treatment effect, otherwise, there is overall lap and big noise in two datasets.

}
\value{
return a numeric value
}
\references{
Yuan-De Tan Anita M. Chandler, Arindam Chaudhury, and Joel R. Neilson(2015) A Powerful Statistical Approach for Large-scale Differential Transcription Analysis. \emph{Plos One}. 2015 DOI: 10.1371/journal.pone.0123658. 
}
\author{
Yuan-De Tan
\cr
\email{tanyuande@gmail.com}
}

\seealso{
\code{\link{omega}}
}
\examples{

X<-c(112,122, 108,127)
Y<-c(302, 314, 322, 328)
rhov(XA=X,XB=Y)	
#[1] 3.062173
}

\keyword{ gap }
\keyword{ noise }
