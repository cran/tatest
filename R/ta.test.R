ta.test<-function(X,nci=NULL,na,nb,alpha=0.05,paired=FALSE,eqv=TRUE,
LOG=c("NULL","LOG2","LOG","LOG10"),
 alternative = c("two.sided", "less", "greater"), distr="norm"){
 	
tatest <-
function(X,na,nb,eqv=FALSE,paired=FALSE,alpha=0.05,
alternative = c("two.sided", "less", "greater"), 
LOG=c("NULL","LOG2","LOG","LOG10"),distribution="norm"){
if(is.null(X)){
XA<-X[1:na]
XB<-X[(1+na):length(X)]
}else{
XA<-X[,1]
XB<-X[,2]	
}
wa<-omega(XA=XA,XB=XB,na=na,nb=nb,m=2000,alpha=alpha,distr=distribution)
wa<-round(wa,1)
conf.level<-1-alpha

rho<-rhov(XA=XA,XB=XB)
tt<-ttest(x=XA,y=XB,alternative=alternative, paired=paired,var.equal=eqv,conf.level,LG=LOG)

#print(tt)

tb<-tt$statistic
df<-tt$parameter

ta<-tb*rho/wa

alternative <- match.arg(alternative)

if(alternative == "less"){
        pval <- pt(ta, df)
        cint <- c(-Inf, ta + qt(conf.level, df))
}else if(alternative == "greater") {
 pval <- pt(ta, df, lower.tail = FALSE)
        cint <- c(ta - qt(conf.level, df), Inf)
}else{

pval<-2*pt(-abs(ta),df)
        alpha <- 1 - conf.level
        cint <- qt(1 - alpha/2, df)
        cint <- ta + c(-cint, cint)

}

tt$p.value<-pval
tt$statistic<-ta
attr(cint, "conf.level") <- conf.level
tt$conf.int<-cint
tt$omega<-wa
tt$rho<-rho

return(tt)

}

mtatest <-
function(X,nci,na,nb,paired=FALSE,eqv=TRUE,
LOG=c("NULL","LOG2","LOG","LOG10"),
 alternative= c("two.sided", "less", "greater"),distribution="norm"){
nr<-nrow(X)
nc<-ncol(X)	

XA<-X[,c((nci+1):(nci+na))]
XB<-X[,c((nci+na+1):nc)]
#print(XA[1:3,])
#print(XB[1:3,])
XA<-apply(XA,2,as.numeric)
XB<-apply(XB,2,as.numeric)

tt<-t01<-t05<-t<-rep(NA,nr)

pv01<-pv05<-pv<-rep(1,nr)
rho<-rep(0,nr)

w<-momega(dat=X[,c((nci+1):nc)],na=na,nb=nb, distr=distribution)
w01<-w[1]
w05<-w[2]

for(i in 1:nr){
tt[i]<-ttest(XA[i,],XB[i,],var.equal=eqv,paired=FALSE,alternative=alternative,LG=LOG)$statistic

rho[i]<-rhov(XA[i,],XB[i,])

t01[i]<-tt[i]*rho[i]/w01
t05[i]<-tt[i]*rho[i]/w05
	
pv01[i]<-2*pt(-abs(t01[i]),df=(na+nb-2))
pv05[i]<-2*pt(-abs(t05[i]),df=(na+nb-2))
}	
#print(w)
res<-cbind(X,t01,pv01,t05,pv05)
return(res)
}

if(is.null(nci)){
res<-tatest(X=X,na=na,nb=nb,eqv=eqv, paired=paired, alpha=alpha,
alternative= alternative, LOG=LOG,  distribution=distr)	
}else{
res<-mtatest(X=X,nci=nci,na=na,nb=nb,eqv=eqv,
LOG=LOG, alternative=alternative,distribution=distr)	
}
return(res)
}
