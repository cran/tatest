\name{omega}
\alias{omega}
\title{
Omega calcularion
}
\description{
Omega (\eqn{\omega}) function is a function that is used to estimate omega using simulate  null data from negative bionomial distribution. Omege is a null rho that is used as a threshold for real rho. Simulation is dependent on the original data. 
}
\usage{
omega(XA, XB, na, nb, m, alpha = 0.05, distr = "norm")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{XA}{
a numeric vector values of condition A
}
  \item{XB}{
a numeric vector values of condition b
}
  \item{na}{
a numeric value, sample size in A.
}
  \item{nb}{
a numeric value, sample size in B
}
  \item{m}{
simulation number specified.
}
  \item{alpha}{
a numeric value, statistical cutoff. Default is 0.05.
}
  \item{distr}{
data distribution specified for data. For current version, we consider three distributions: normal, negative binomial (NB) and uniform. Default distribution is "norm". If user believe that the data follow negative binomial distribution, then set distr="NB" or "negative binomial" or if the data follow uniform distribution, then set distr="unif"or "uniform". 
}
}
\details{
This function is to use  simulated null data to calculate omega value with rho  = 1. 
}
\value{
return a numeric value
}
\references{
Yuan-De Tan Anita M. Chandler, Arindam Chaudhury, and Joel R. Neilson(2015)  A Powerful Statistical Approach for Large-scale Differential Transcription Analysis. \emph{Plos One}. 2015 DOI: 10.1371/journal.pone.0123658. 
}
\author{
Yuan-De Tan
\email{tanyuande@gmail.com}
}

\seealso{
\code{\link{momega}}, \code{\link{rhov}}
}
\examples{
X<-c(112,122,108,127)
Y<-c(302, 314,322,328)
omega(XA=X, XB=Y, na=4, nb=4, m=2000, alpha = 0.05)
#[1] 0.9055152
 omega(XA=X, XB=Y, na=4, nb=4, m=2000, alpha = 0.05,distr="NB")
#[1] 0.8995424
omega(XA=X, XB=Y, na=4, nb=4, m=2000, alpha = 0.05,distr="uniform")
#[1] 0.97194
}

\keyword{thresholds}
\keyword{rho}
