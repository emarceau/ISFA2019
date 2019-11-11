# Cours : Lyon
# Semestre : A2019
#
# Loi gamma bivariée CRMM
#
be1<-0.1
be2<-0.1
a1<-2
a2<-2
g0<-0.5
#
set.seed(20191108)
#
mm<-1000
V0<-runif(mm)
V1<-runif(mm)
V2<-runif(mm)
Y0<-qgamma(V0,g0,1)
Y1<-qgamma(V1,a1-g0,1)
Y2<-qgamma(V2,a2-g0,1)
X1<-(Y0+Y1)/be1
X2<-(Y0+Y2)/be2
mean(X1)
mean(X2)
Fn<-(1:mm)/mm
matplot(sort(X1),cbind(Fn,pgamma(sort(X1),a1,be1)),type="s",xlab="x",ylab="F1(x)")

matplot(sort(X2),cbind(Fn,pgamma(sort(X2),a2,be2)),type="s",xlab="x",ylab="F2(x)")

S<-X1+X2
St<-sort(S)
Fn<-(1:mm)/mm
Sind<-qgamma(V1,a1,be1)+qgamma(V2,a2,be2)
Scom<-qgamma(V1,a1,be1)+qgamma(V1,a2,be2)

kap<-(1:mm)/(mm+1)
matplot(kap,cbind(St,sort(Sind),sort(Scom)),type="s",xlab="kappa",ylab="VaR(S;kappa")






