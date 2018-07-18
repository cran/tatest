\name{Ecadherin}
\alias{Ecadherin}
\docType{data}
\title{
A microarray data
}
\description{
Analysis of the effects of loss of E-cadherin and cell adhesion on human mammary epithelial cells
}
\usage{data("Ecadherin")}
\format{
  A data frame with 5000 observations on the following 7 variables.
  \describe{
    \item{\code{geneid}}{a factor with levels }
    \item{\code{shCntrl_rep1}}{a numeric vector}
    \item{\code{shCntrl_rep2}}{a numeric vector}
    \item{\code{shCntrl_rep3}}{a numeric vector}
    \item{\code{shEcad_rep1}}{a numeric vector}
    \item{\code{shEcad_rep2}}{a numeric vector}
    \item{\code{shEcad_rep3}}{a numeric vector}
  }
}
\details{
Loss of the epithelial adhesion molecule E-cadherin is thought to enable metastasis by disrupting intercellular contacts-an early step in metastatic dissemination. To further investigate the molecular basis of this notion, two methods were used to inhibit E-cadherin function that distinguish between E-cadherin's cell-cell adhesion and intracellular signaling functions. While the disruption of cell-cell contacts alone does not enable metastasis, the loss of E-cadherin protein does, through induction of an epithelial-to-mesenchymal transition, invasiveness and anoikis-resistance. In addition, gene expression analysis shows that E-cadherin loss results in the induction of multiple transcription factors, at least one of which, Twist, is necessary for E-cadherin loss-induced metastasis. These findings indicate that E-cadherin loss in tumors contributes to metastatic dissemination by inducing wide-ranging transcriptional and functional changes.
This dataset consists 5000 of 22277  gene libraries, two conditions shCntrl and shEcad, each having 3 replicates.  With access number GSE9691, you can obtain a whole data of 22277 gene libraries in GEO.
}


\examples{
data(Ecadherin)
}
\keyword{dataset}
