# Lyon A2018
#
# Simulation loi binomiale negative composée
#
#
# loi de M : binomiale negative
qq<-0.5
rr<-2
EM<-rr*(1-qq)/qq
VarM<-EM/qq

# loi de B : lognormale
mu<-log(10)-0.405
sig<-0.9
EB<-exp(mu+(sig^2)/2)
EB2<-exp(2*mu+2*(sig^2))
VarB<-EB2-EB^2

EX<-EM*EB
VarX<-EM*VarB+VarM*(EB^2)

c(EM,EB,EX)
c(VarM,VarB,VarX)


# simulation de réalisations la loi binomiale négative composée

# fixer la valeur source 
set.seed(2018)

# nb de simulation
nsim<-100000

# vecteurs des réalisations de M et X
vM<-rep(0,nsim)
vX<-rep(0,nsim)

# algo de simulation des réalisations de M et X
for (i in 1:nsim)
{
  U<-runif(1)
  vM[i]<-qnbinom(U,rr,qq)
  if (vM[i]>0) vX[i]<-sum(qlnorm(runif(vM[i]),mu,sig))
}

# valeurs exactes
c(EM,VarM,EB,VarB,EX,VarX)

# valeurs approximatives
mean(vM)
var(vM)
mean(vX)
var(vX)

# courbe des valeurs approximatives de FX(x)
plot.ecdf(vX,xlab="x")

# 100 premières valeurs des réalisations de (M,X)
cbind(1:100,vM[1:100],vX[1:100])

# approximation de FX(x)
x<-0
Fx<-sum((vX<=x))/nsim
c(x,Fx)


# approximation de VaR et TVaR de X
vXs<-sort(vX)

#vkappa<-c(0.5,0.9,0.99,0.999)
vkappa<-(1:999)/1000
vVaRX<-quantile(vX,probs=vkappa,type=1)
#cbind(vkappa,vVaRX,vXs[vkappa*nsim])

nkap<-length(vkappa)
TVaRX<-rep(0,nkap)


for (i in 1:nkap)
{
  dum<-vkappa[i]*nsim+1
  TVaRX[i]<-sum(vXs[dum:nsim])/(nsim*(1-vkappa[i]))
}

#cbind(vkappa,vVaRX,vXs[vkappa*nsim],TVaRX)

matplot(vkappa,cbind(vVaRX,TVaRX),type="l")


# courbe des valeurs approximatives de VaR

vu<-(1:999)/1000
vquanX<-quantile(vX,probs=vu,type=1)
plot(vu,vquanX,type="l")
