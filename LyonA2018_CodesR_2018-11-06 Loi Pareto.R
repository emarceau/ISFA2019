# Lyon A2018
# Date : mardi 2018-11-06
# Loi Pareto

aa<-2.5
lam<-10*(aa-1)


EX<-lam/(aa-1)
EX

u<-(1:999)/1000
VaR<-lam*(1/((1-u)^(1/aa))-1)
plot(u,VaR,type="l",xlab="kappa",ylab="VaR",main="VaR de la Loi Pareto")

VaRgamma<-qgamma(u,2,1/5)
matplot(u,cbind(VaR,VaRgamma),type="l",xlab="kappa",ylab="VaR",main="VaR de la Loi Pareto et de la loi Gamma")

set.seed(2018)
nsim<-100000
U<-runif(nsim)
X<-lam*(1/((1-U)^(1/aa))-1)
n<-1:nsim
plot(n,cumsum(X)/n,type="l")


