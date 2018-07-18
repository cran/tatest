\name{ttest}
\alias{ttest}
\title{
A modified t-test 
}
\description{
ttest is a modified t-test of Gosset(1098). 
}
\usage{
ttest(x, y = NULL, alternative = c("two.sided", "less", "greater"), 
mu = 0, paired = FALSE, var.equal = FALSE, conf.level = 0.95, 
LG = c("NULL", "LOG2", "LOG", "LOG10"), ...)
}

\arguments{
  \item{x}{
a (non-empty) numeric vector of data values.
}
  \item{y}{
an optional (non-empty) numeric vector of data values. For \eqn{t_\alpha}-test, y is required
}
  \item{alternative}{
a character string specifying the alternative hypothesis, must be one of "two.sided" (default), "greater" or "less". You can specify just the initial letter.
}
  \item{mu}{
a number indicating the true value of the mean (or difference in means if you are performing a two sample test).
}
  \item{paired}{
a logical indicating whether you want a paired t-test.
}
  \item{var.equal}{
a logical variable indicating whether to treat the two variances as being equal. If TRUE then the pooled variance is used to estimate the variance otherwise the Welch (or Satterthwaite) approximation to the degrees of freedom is used.
}
  \item{conf.level}{
confidence level of the interval.
}
  \item{LG}{
logarithm: if x and y are not taken logarithm convertion, the LG=NULL , if x and y are converted into logarithm with base 10, 2 or natural logarithm, then correspondingly, LG="LOG10", 
}
  \item{\dots}{
further arguments to be passed to or from methods.
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
}
\value{
A list with class "htest" containing the following components:
\cr
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
the p-value for the test.
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
a character string giving the name(s) of the data.
}
\references{
Gosset, W.S, "Probable error of a correlation coefficient". Biometrika. 1908,6(2/3): 302-310.  
}
\author{
Yuande Tan
\cr
\email{tanyuande@gmail.com}
}
\seealso{
 \code{\link{stats}{t.test}}, \code{\link{ta.test}}
}
\examples{
# The following example is log10 transformed data	
YA<-c(4.1455383,	4.7186418,4.6483316)
YB<-c(6.9184465,	5.609667,6.5920944)
t.test(x=YA, y=YB)
ttest(YA,YB,alternative = "two.sided", LG = "NULL")
ttest(YA,YB,alternative = "two.sided", LG = "LOG2")
ttest(YA,YB,alternative = "two.sided", LG = "LOG10")
}
\keyword{ t.test }
\keyword{ ttest }
