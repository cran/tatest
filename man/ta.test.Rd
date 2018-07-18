\name{ta.test}
\alias{ta.test}
\title{
\eqn{t_\alpha}-test for testing difference between two conditions with small samples
}
\description{
\eqn{t_\alpha}-test is a modified t-test of Gosset(1908). 
}
\usage{
ta.test(X, nci=NULL, na, nb, alpha=0.05, paired=FALSE, eqv=TRUE,
LOG=c("NULL","LOG2","LOG","LOG10"),
 alternative = c("two.sided", "less", "greater"), distr="norm")
}
\arguments{
  \item{X}{
a numeric vector data or a numeric matrix (two columns) dataset with two conditions or groups A and B having sample size na and nb, respectively. X is also a datasheet of multiple variables (for example, genes ) with two groups A and B and information columns(see an example Ecadherin)
}
\item{nci}{
a numeric int value, column numbers for strings or information of data. If X is a numeric vector data or a numeric matrix (two columns) 
dataset, then nci=NULL, otherwise, nci is non-zeror int, meaning that at least one column is reserved as information of data
}
  \item{na}{
numeric constant, sample size in condition A
}
  \item{nb}{
numeric constant, sample size in condition B
}

  \item{alpha}{
a numeric constant, statistical cutoff for significance level. If dataset X is multiple variables (rows) (variable number is over 1000), then alpha is taken default values of 0.05 and 0.01.
}
  \item{paired}{
a logical indicating whether you want a paired t-test.
}
  \item{eqv}{
a logical variable indicating whether to treat the two variances as being equal. If TRUE then the pooled variance is used to estimate the variance otherwise the Welch (or Satterthwaite) approximation to the degrees of freedom is used.
}

  \item{LOG}{
logarithm: if data are not taken logarithm convertion, then LOG=NULL , if data are converted into logarithm with base 10, 2 or natural logarithm, then correspondingly, LOG="LOG10", "LOG2" or "LOG".
}
  \item{alternative}{
a character string specifying the alternative hypothesis, must be one of "two.sided" (default), "greater" or "less". User can specify just the initial letter.
}
  \item{distr}{
data distribution for estimating omega (\eqn{\omega}). Default is normal distribution.  In current version, we just consider three types of distributions: "negative distribution"(or "NB"), "norm"(or "normal"), and "unif"("uniform").  User can choose one of them according to user's data distribution.
}
}
\details{
The formula interface is only applicable for the 2-sample tests. 
\cr
alternative = "greater" is the alternative that x has a larger mean than y.
\cr
If paired is TRUE then both x and y must be specified and they must be the same length. Missing values are silently removed (in pairs if paired is TRUE). 
\cr
If var.equal is TRUE then the pooled estimate of the variance is used. By default, if var.equal is FALSE, then the variance is estimated separately for both groups and the Welch modification to the degrees of freedom is used.
\cr
If the input data are effectively constant (compared to the larger of the two means) an error is generated.
\cr
When pooled variance is larger, ttest is equvalent to t.test; When pooled variance is small, t-tavalue of ttest is smaller than than that of t.test; When the pooled variance is zeror, t-value of t.test is not defined but tvalue of ttest is still defined.  
\cr
Omega is dependent on alpha value.
\cr
ta.test can be used to test for differentially expressed genes in microarray data, proteomic data and RPPA data but cannot be used to count data of RNA-seq reads. 
}
\value{
If X is single-hypothesis test data, then  a list similar to t.test containing the following components:
statistic
\cr
the value of the t-statistic.
\cr
parameter
\cr
the degrees of freedom for the t-statistic.
\cr
p.value
\cr
the p-alue for the test.
\cr
conf.int
\cr
a confidence interval for the mean appropriate to the specified alternative hypothesis.
\cr
estimate
\cr
the estimated mean or difference in means depending on whether it was a one-sample test or a two-sample test.
\cr
null.value
\cr
the specified hypothesized value of the mean or mean difference depending on whether it was a one-sample test or a two-sample test.
\cr
alternative
\cr
a character string describing the alternative hypothesis.
\cr
method
\cr
a character string indicating what type of t-test was performed.
\cr
data.name
\cr
a character string giving the name(s) of the data.
\cr
omega
a threshold for rho in tatest. 
\cr
If X is microarray, RPPA data, then a list is datasheet including t01,pv01, t05 and p05 columns
}
\references{
Gosset, W.S, "Probable error of a correlation coefficient". Biometrika. 1908,6 (2/3): 302-310.  
}
\author{
Yuande Tan
\cr
\email{tanyuande@gmail.com}
}
\note{
LOG choice must be the same with log transformation of data.
}

\seealso{
\code{\link[stats]{t.test}},\code{\link{ttest}}
}
\examples{
X<-c(112,122,108,127)
Y<-c(302, 314,322,328)
dat<-cbind(X,Y)
ta.test(X=dat,nci=NULL, na=4,nb=4, eqv=TRUE, 
alpha=0.05, LOG="NULL", alternative = "two.sided")

data(Ecadherin)
res<-ta.test(X=Ecadherin[1:2000,],nci=1, na=3,nb=3, eqv=TRUE, 
LOG="LOG2", alternative = "two.sided")
}

\keyword{ rhov}
\keyword{ ttest}
