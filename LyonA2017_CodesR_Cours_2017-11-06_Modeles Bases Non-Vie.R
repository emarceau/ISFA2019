# Lyon
# A2017
# Cours Lundi 2017-11-06

# Loi Poisson composée avec sinistre individuel de loi gamma
# Loi de M : Poisson 
lambda=0.5
EM<-lambda
VarM<-lambda
# Loi de B : Gamma
alp<-5
bet<-1/200
EB<-alp/bet
VarB<-EB/bet
#
EX<-EM*EB
VarX<-EM*VarB+VarM*(EB^2)
EX
VarX
#
# Fonction de répartition de X
#
Fpoisgamma<-function(x,la,aa,bb,kmax=1000)
{
  p0<-dpois(0,la)
  vk<-1:kmax
  pk<-dpois(vk,la)
  vprob<-pgamma(x,aa*vk,bb)
  FX<-p0+sum(pk*vprob)
  return(FX)
}

Fpoisgamma(x=3200,la=lambda,aa=alp,bb=bet,kmax=1000)

vx<-(0:20)*500
long<-length(vx)
vFx<-rep(0,long)
for(i in 1:long)
{
  vFx[i]<-Fpoisgamma(x=vx[i],la=lambda,aa=alp,bb=bet,kmax=1000)
}

plot(c(0,vx),c(0,vFx),type="l",xlab="x",ylab="F_X(x)",main="Fonction de répartition de X - Loi Poisson composée sinistres de loi gamma")

#
# VaR et TVaR
#
# on utilise cette approche pour kappa > F_X(0)

Fpoisgamma(0,la=lambda,aa=alp,bb=bet,kmax=1000)

kappa<-0.9999
f<-function(x) abs(Fpoisgamma(x,la=lambda,aa=alp,bb=bet,kmax=1000)-kappa)
res<-optimize(f, c(0,10000),tol=0.000000001)
res
VaRX<-res$minimum
Fpoisgamma(VaRX,la=lambda,aa=alp,bb=bet,kmax=1000)

TVaRpoisgamma<-function(u,la,aa,bb,kmax=1000,bornes=c(0,10000))
{
  kappa<-u
  f<-function(x) abs(Fpoisgamma(x,la=lambda,aa=alp,bb=bet,kmax=1000)-kappa)
  res<-optimize(f, bornes,tol=0.000000001)
  VaR<-res$minimum
  vk<-1:kmax
  pk<-dpois(vk,la)
  vEtronc<-(1-pgamma(VaR,aa*vk+1,bb))*aa/bb*(vk)
  TVaR<-sum(pk*vEtronc)/(1-u)
  return(c(VaR,TVaR))
}

TVaRpoisgamma(0.9999,la=lambda,aa=alp,bb=bet,kmax=1000,bornes=c(0,10000))

# ________________________________________________________
#
# Loi Binomiale négative composée (sinistres de loi lognormale)
# loi de M: binomiale négative
rr<-0.5
qq<-0.5
EM<-rr*(1-qq)/qq
VarM<-EM/qq

# Loi de B : LNormale
mu<-log(1000)-0.32
sig<-0.8
EB<-exp(mu+(sig^2)/2)
EB2<-exp(2*mu+2*(sig^2))
VarB<-EB2-(EB^2)
# X:
EX<-EM*EB
VarX<-EM*VarB+VarM*(EB^2)
EX
VarX

# Simulons :
set.seed(2017)
nsim<-100000
vM<-rep(0,nsim)
vX<-rep(0,nsim)
for(i in 1:nsim)
{
  U<-runif(1)
  vM[i]<-qnbinom(U,rr,qq)
  if (vM[i]>0) 
  {
    vU<-runif(vM[i])
    vX[i]<-sum(qlnorm(vU,mu,sig))
  }
}
#cbind(vM,vX)

plot.ecdf(vX)

vkap<-(1:9999)/10000

VaRX<-quantile(vX,prob=vkap,type=1)
plot(vkap,VaRX,type="l",xlab="kappa",ylab="VaRX")



