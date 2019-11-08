# Lyon
# A2017
# Cours Vendredi 2017-11-03
#
# Loi exponentielle

bet<-1/10
vu<-(1:999)/1000
VaRX<-qexp(vu,bet)
EX<-1/bet
vEX<-rep(EX,999)
TVaRX<-VaRX+EX
matplot(vu,cbind(VaRX,vEX,TVaRX),type="l",main="VaR et TVaR de la loi exponentielle")
#
# Loi gamma
#
alp<-5
bet<-1/2
vu<-(1:999)/1000
VaRX<-qgamma(vu,alp,bet)
EX<-alp/bet
vEX<-rep(EX,999)
TVaRX<-EX*(1-pgamma(VaRX,alp+1,bet))/(1-vu)
matplot(vu,cbind(VaRX,vEX,TVaRX),type="l",main="VaR et TVaR de la loi gamma",xlab="u",ylab="VaR et TVaR")
#
#
# Loi binomiale
#
qq<-0.0017
bb<-100000
EX<-bb*qq
EX
kap<-0.995
VaRX<-bb*qbinom(kap,1,qq)
VaRX
TVaRX<-EX/(1-kap)
TVaRX
vn<-c(1,10,100,1000,10000,100000,1000000)
vVaRN<-qbinom(kap,vn,qq)
vVaRN
vVaRS<-bb*vVaRN
vVaRS
nono<-length(vn)
vTVaRN<-rep(0,nono)
for (i in 1:nono)
{
  vk<-0:vn[i]
  partie1<-sum(vk*dbinom(vk,vn[i],qq)*1*(vk>vVaRN[i]))
  partie2<-vVaRN[i]*(pbinom(vVaRN[i],vn[i],qq)-kap)
  vTVaRN[i]<-(partie1+partie2)/(1-kap)  
}
cbind(kap,vVaRN,vTVaRN)
vTVaRS<-bb*vTVaRN
round(cbind(vn,vVaRS/vn,vTVaRS/vn,TVaRX-vTVaRS/vn),2)

round(cbind(vn,vVaRS/vn,vTVaRS/vn,TVaRX-vTVaRS/vn),2)

# Aggrégation des coûts par contrat

# Loi de Xi : exponentielle(bet)
# nb de contrats : n
# Sn =  coûts totaux pour n contrats
# Wn =  coûts par contrat pour un ptf de n contrats
# Loi de Sn : gamma(n,bet)
# Loi de Wn : gamma(n,bet*n)
bet<-1/10
vn<-10^(0:4)
vx<-(0:100)/5
matfWn<-matrix(0,101,5)
for (i in (1:5))
{
  matfWn[,i]<-dgamma(vx,vn[i],bet*vn[i])
}

matplot(vx,matfWn,type="l",xlab="x",ylab="fonction de densité de Wn")
#
#
# générateur de réalisations U(j) de la v.a. U \sim Unif(0,1)
#
aa<-41358 
mm<-2^31-1

x0<-2017
nn<-1000000
vx<-rep(0,nn)
vx[1]<-(aa*x0)%%mm
for (i in 2:nn)
{
  vx[i]<-(aa*vx[i-1])%%mm
}
 
#cbind(1:nn,vx,vx/mm)
vU<-vx/mm

v1<-vU[1:(nn-1)]
v2<-vU[2:nn]
#plot(v1,v2)

mean(vU)
mean(qexp(vU))
mean(qgamma(vU,2,1/5))

# 
#
#
# somme de 2 v.a. indépendantes
#
# loi de X1: gamma(a1,bet)
# loi de X2: gamma(a2,bet)
a1<-2.5
a2<-1.5
bet<-1/10
nsim<-10^6
set.seed(2017)

matU<-matrix(runif(nsim*2),nsim,2,byrow=T)
#matU

X1<-qgamma(matU[,1],a1,bet)
X2<-qgamma(matU[,2],a2,bet)

matX<-cbind(X1,X2)
S<-X1+X2
#cbind(1:nsim,X1,X2,S)

mean(S)
mean(1*(S>50))
quantile(S,c(0.5,0.9),type=1)

EX1<-a1/bet
EX2<-a2/bet
ES<-EX1+EX2
ES

xx<-50
mean(1*(S>xx))
1-pgamma(xx,a1+a2,bet)


xx<-100
mean(1*(S>xx))
1-pgamma(xx,a1+a2,bet)

kap<-c(0.5,0.9,0.99,0.999)
quantile(S,kap,type=1)
qgamma(kap,a1+a2,bet)

kap1<-0.99999
VaRSapp<-quantile(S,kap1,type=1)
TVaRSapp<-sum(S*1*(S>VaRSapp))/nsim/(1-kap1)
VaRS<-qgamma(kap1,a1+a2,bet)
TVaRS<-ES*(1-pgamma(VaRS,a1+a2+1,bet))/(1-kap1)
c(kap1,VaRSapp,VaRS,TVaRSapp,TVaRS)

VaRS















