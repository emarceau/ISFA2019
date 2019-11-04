# Cours : Lyon
# Semestre : A2019
#
# Loi de X : exponentielle
#
bet<-1/10
EX<-1/bet
EX
x1<-20
F<-pexp(x1,bet)
F
1-exp(-x1*bet)
vx<-0:100
vF<-pexp(vx,bet)
plot(vx,vF,type="l",xlab = "x",ylab = "F(x)")
#
vkap<-(1:9999)/10000
vVaR<-qexp(vkap,bet)
vTVaR<-vVaR+EX
vEX<-rep(EX,length(vkap))
matplot(vkap,cbind(vVaR,vTVaR,vEX),type="l",xlab = "kappa",ylab = "VaR")
##
#
#
#
# Loi de X : gamma
#
alp<-2
bet<-1/5
EX<-alp/bet
EX
x1<-20
F<-pgamma(x1,alp,bet)
F
vx<-0:100
vF<-pgamma(vx,alp,bet)
plot(vx,vF,type="l",xlab = "x",ylab = "F(x)")
#
vkap<-(1:9999)/10000
vVaR<-qgamma(vkap,alp,bet)
vTVaR<-EX*(1-pgamma(vVaR,alp+1,bet))/(1-vkap)
vEX<-rep(EX,length(vkap))
matplot(vkap,cbind(vVaR,vTVaR,vEX),type="l",xlab = "kappa",ylab = "VaR")
#
#
# Loi Poisson 
# 
lam<-5
EM<-lam
vk<-0:50
vfM<-dpois(vk,lam)
plot(vk,vfM,type="h",xlab="k",ylab="Pr(M=k)")
vkap<-(1:9999)/10000
vVaR<-qpois(vkap,lam)
plot(vkap,vVaR,type="s",xlab="kappa",ylab="VaR")
#
tvarpois<-function(u,lambd)
{
  VaR<-qpois(u,lambd)
  EX<-lambd
  vm<-0:VaR
  EXt<-sum(vm*dpois(vm,lambd))
  TVaR<-(EX-EXt+VaR*(ppois(VaR,lambd)-u))/(1-u)
  return(TVaR)
}
#
tvarpois(0.5,lam)
nb<-length(vkap)
vTVaR<-rep(0,nb)
for(j in 1 :nb)
{
  vTVaR[j]<-tvarpois(vkap[j],lam)
}  
vEM<-rep(EM,nb)
matplot(vkap,cbind(vVaR,vTVaR,vEM),type="s",xlab = "kappa",ylab = "VaR et TVaR")
#


