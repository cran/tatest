momega<-function(dat,na,nb, distr="norm"){

m<-nrow(dat)
MX<-apply(dat,1,mean)

SDX<-apply(dat,1,var)
rho<-pv<-rep(0,m)
SXA<-matrix(NA,m,na)
SXB<-matrix(NA,m,nb)
i<-seq(length(MX))
if(distr=="NB"||distr=="negtive binomial"){
for(i in 1:m){
 SXA[i,]<-rnbinom(na,mu=MX[i],size=SDX[i])
 SXB[i,]<-rnbinom(nb,mu=MX[i],size=SDX[i])
 rho[i]<-rhov(SXA[i,],SXB[i,])
 pv[i]<-t.test(SXA[i,],SXB[i,])$p.value
 }
} else if (distr=="norm"||distr=="normal"){
for(i in 1:m){
SXA[i,]<-abs(rnorm(na,mean=MX[i],sd=sqrt(SDX[i])))
SXB[i,]<-abs(rnorm(nb,mean=MX[i],sd=sqrt(SDX[i])))
#print(cbin(SXA[i,],SXB[i,]))
rho[i]<-rhov(SXA[i,],SXB[i,])
pv[i]<-t.test(SXA[i,],SXB[i,])$p.value
 }
} else if(distr=="unif"||distr=="uniform"){
for(i in 1:m){	
SXA[i,]<-runif(na)*MX
SXB[i,]<-runif(nb)	*MX
rho[i]<-rhov(SXA[i,],SXB[i,])
pv[i]<-t.test(SXA[i,],SXB[i,])$p.value
}
}

res<-rbind(pv,rho)
res<-res[,order(pv)]
m05<-which(res[1,]<0.05)
nrho05<-res[2,m05]
#nrho05<-sort(nrho05)
w05<-mean(nrho05)

m01<-which(res[1,]<0.01)
nrho01<-res[2,m01]
#nrho01<-sort(nrho01)
w01<-mean(nrho01)


return(c(w01,w05))
}