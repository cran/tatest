\name{tatest-package}
\alias{tatest-package}
\docType{package}
\title{
Two-Group \eqn{T_\alpha}-Test
}
\description{
The \eqn{t_\alpha}-test is a modified two-sample or two-group t-test of Gosset(1908).
}
\details{
The package tatest consists of three t-tests: ttest, tatest and mtatest. ttest, rhov,and omega construct tatest, ttest and rhov and momega construct mtatest. \eqn{t_\alpha}-test is used to test for a single null hypothesis  and \eqn{mt_\alpha}-test is used to identify differential expressed genes in microarray data and RPPA data or proteomic data. 
}
\author{
Yuande Tan
\cr
\email{tanyuande@gmail.com}
Yuande Tan
\cr
Maintainer: \email{tanyuande@gmail.com}
}


\keyword{ tatest}
\seealso{
\code{\link[stats]{t.test}}
}
\examples{
X<-c(112,122,108,127)
Y<-c(302, 314,322,328)
dat<-cbind(X,Y)
ta.test(X=dat,na=4,nb=4, eqv=TRUE, paired=FALSE, 
alpha=0.05, alternative = "two.sided",LOG="NULL")
t.test(x=X,y=Y)
ttest(X,Y,alternative = "two.sided", LG = "NULL")
}
