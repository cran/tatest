\name{momega}
\alias{momega}

\title{
 Null rho for testing multiple null hypotheses
}
\description{
Momega is a function that is used to estimate omega for testing multiple null hypotheses. 
Omege is a null \eqn{\rho} that is used as a threshold for \eqn{\rho} from real data. Since this omega is majorly applied 
to micoarray, we named this function  as momega.
}
\usage{
momega(dat, na, nb, distr = "norm")
}

\arguments{
  \item{dat}{
a matrix dataset containing string columns and two-condition numeric columns
}
  \item{na}{
numeric value, size of sample 1 or number of replicates in condition A.
}
  \item{nb}{
numeric value, size of sample 1 or number of replicates in condition A.
}
  \item{distr}{
string, data distribution.
}
}
\details{
This function is to use  null data to calculate omega value with rho  = 1. 
}
\value{
return a numeric value
}
\references{
Yuan-De Tan Anita M. Chandler, Arindam Chaudhury, and Joel R. Neilson (2015)  A Powerful Statistical Approach for Large-scale Differential Transcription Analysis. \emph{Plos One}. 2015 DOI: 10.1371/journal.pone.0123658. 
}
\author{
Yuan-De Tan
\email{tanyuande@gmail.com}
}
\seealso{
\code{\link{omega}},\code{\link{rhov}}
}
\examples{
data(Ecadherin)
Ecad<-as.data.frame(Ecadherin)
Ecad1<-Ecad[,-1]
w<-momega(dat=Ecad1[1:2000,], na=3, nb=3, distr = "norm")
}

\keyword{ null rho}
\keyword{ multiple tests }
