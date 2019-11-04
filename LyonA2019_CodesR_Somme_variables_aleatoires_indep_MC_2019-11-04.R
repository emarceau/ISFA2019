# Cours : Lyon
# Semestre : A2019
#
# loi de X1: gamma
a1<-2
b1<-1/5
# loi de X2: gamma
a2<-0.5
b2<-1/20
#
EX1<-a1/b1
EX2<-a2/b2
ES<-EX1+EX2
c(EX1,EX2,ES)
mm<-100000
set.seed(20191104)
matU<-matrix(runif(2*mm),mm,2,byrow=T)

set.seed(20191104)
X1<-qgamma(matU[,1],a1,b1)
X2<-qgamma(matU[,2],a2,b2)
S<-X1+X2
ESa<-mean(S)
c(ES,ESa)
varX1<-EX1/b1
varX2<-EX2/b2
varS<-varX1+varX2
varSa<-var(S)
c(varS,varSa)
St<-sort(S)
X1t<-sort(X1)
X2t<-sort(X2)
Fn<-(1:mm)/mm
plot(St,Fn,type="s",xlab="x",ylab="Fn(x)")
vkap<-1:(mm-1)/mm
vVaRS<-St[1:(mm-1)]

vVaRX1<-X1t[1:(mm-1)]
vVaRX2<-X2t[1:(mm-1)]
matplot(vkap,cbind(vVaRS,vVaRX1+vVaRX2),type="s",xlab="kappa",ylab="VaRS")

vTVaRS<-(cumsum(St[mm:2])/(1:(mm-1)))[(mm-1):1]
vTVaRX1<-(cumsum(X1t[mm:2])/(1:(mm-1)))[(mm-1):1]
vTVaRX2<-(cumsum(X2t[mm:2])/(1:(mm-1)))[(mm-1):1]

matplot(vkap,cbind(vVaRS,vTVaRS,vTVaRX1+vTVaRX2),type="s",xlab="kappa",ylab="VaRS et TVaRS")
