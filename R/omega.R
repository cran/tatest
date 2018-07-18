omega <-
function(XA,XB,na,nb,m,alpha=0.05,distr="norm"){
dat<-c(XA,XB)
MX<-mean(dat)

SDX<-var(dat)
rho<-pv<-rep(0,m)
SXA<-matrix(NA,m,na)
SXB<-matrix(NA,m,nb)

if(distr=="NB"||distr=="negtive binomial"){
for(i in 1:m){
 SXA[i,]<-rnbinom(na,mu=MX,size=SDX)
 SXB[i,]<-rnbinom(nb,mu=MX,size=SDX)
rho[i]<-rhov(SXA[i,],SXB[i,])
pv[i]<-t.test(SXA[i,],SXB[i,])$p.value
 }
} else if (distr=="norm"||distr=="normal"){
for(i in 1:m){
if(max(XA)>0 && min(XA)<0||max(XB)>0 &&min(XB)<0)	{
SXA[i,]<-(rnorm(na,mean=MX,sd=sqrt(SDX)))
SXB[i,]<-(rnorm(nb,mean=MX,sd=sqrt(SDX)))
}else{

SXA[i,]<-abs(rnorm(na,mean=MX,sd=sqrt(SDX)))
SXB[i,]<-abs(rnorm(nb,mean=MX,sd=sqrt(SDX)))	
}
#print(cbin(SXA[i,],SXB[i,]))
rho[i]<-rhov(SXA[i,],SXB[i,])
pv[i]<-t.test(SXA[i,],SXB[i,])$p.value
 }
 }else if(distr=="unif"||distr=="uniform"){
 for(i in 1:m){	
SXA[i,]<-runif(na)*MX
SXB[i,]<-runif(nb)	*MX
rho[i]<-rhov(SXA[i,],SXB[i,])
pv[i]<-t.test(SXA[i,],SXB[i,])$p.value
}
 }


#res<-rbind(pv,rho)
res<-cbind(pv,rho)
#res<-res[,order(pv)]
#n<-which(res[1,]<alpha)
nres<-subset(res, pv<alpha)
nrho<-nres[,2]
#nrho<-res[2,n]

w<-mean(nrho)
return(w)
}
