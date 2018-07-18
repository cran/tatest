rhov <-
function(XA,XB){
########################################################3
# caclulate maximum and minimum values of data XA and XB.
########################################################3
AM<-max(XA)
BM<-max(XB)
Am<-min(XA)
Bm<-min(XB)

#######################################################
# calculate psi=rho here
#######################################################

if(AM>0&&BM>0&&Am>=0&&Bm>=0){
rho<-max(Am/BM,Bm/AM)
}else if(AM<=0&&BM<=0&&Am<0&&Bm<0){
rho<-max(AM/(Bm),BM/(Am))
}else if(Am>=0&&BM<0){
rho<-abs((Am-Bm)/(BM+1-Bm))
}else if(Bm>=0&&AM<0){
rho<-abs((Bm-Am)/(AM+1-Am))
}	else{
C<-max(abs(min(XA)),abs(min(XB)))+1

rho<-max((C+Am)/(C+BM),(C+Bm)/(C+AM))

}


X<-c(XA,XB)

###########################################################
# calculate variances and means of XA, XB and X
############################################################
VA<-var(XA)
VB<-var(XB)
MA<-abs(mean((XA)))
MB<-abs(mean((XB)))
VX<-var(X)
MX<-abs(mean((X)))

################################################################
# caculate zeta and rho
################################################################
# zeta<-log10(1+((VX*MX)+1)/(VA*MA+VB*MB+1))
if((Am>=0&&BM<0)||(Bm>=0&&AM<0)){
zeta=log(1+(((MA+MB)*VX/2)+1)/(MA*VA+MB*VB+1))
}else{
zeta<- log(1+(VX*MX+1)/(VA*MA+VB*MB+1))
}
# print(c(rho,zeta))
rho<-sqrt(rho*zeta)

return(rho)

}
