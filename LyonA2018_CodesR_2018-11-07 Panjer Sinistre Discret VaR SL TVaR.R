# Lyon A2018
#
# Exemple : Algo Panjer + discrétisation
#
# Cas particulier : Algo Panjer Poisson composée

panjer.poisson<-function(lam,fB,xmax)
{
  # Algorithme de Panjer
  # Cas Poisson
  # Loi discrete pour B
  aa<-0
  bb<-lam
  ll<-length(fB)
  fX<-exp(lam*(fB[1]-1))
  fB<-c(fB,rep(0,xmax-ll+1))
  for (i in 1:xmax)
  {
    j<-i+1
    fX<-c(fX,(1/(1-aa*fB[1]))*sum(fB[2:j]*fX[i:1]*(bb*(1:i)/i+aa)))    
  }
  
  return(fX)
}



# Exemple numérique
# Loi de B : discrète 
# Support : {1,2,...,2000}
# fB(0) = 0
vk<-1:2000
fB1<-c(0,(40/(40+vk-1))^5-(40/(40+vk))^5)
sum(fB1)
vk0<-c(0,vk)
EB<-sum(vk0*fB1)
EB

lambda<-5
EM<-lambda

EX<-EM*EB
c(EM,EB,EX)


fX1<-panjer.poisson(lam=5,fB=f1,xmax=2500)
sum(fX1)

FX1<-cumsum(fX1)
plot(0:2500,FX1)

cbind(0:174,FX1[1:175])

# Calcul de la VaR
kappa<-0.99
VaRX1<-sum(FX1<kappa)
VaRX1

# Calcul de la TVaR avec espérance tronquée

long<-length(fX1)-1
vj<-0:long
EX1tronq<-sum(vj*fX1*(vj>VaRX1))
TVaRX1<-(EX1tronq+VaRX1*(FX1[VaRX1+1]-kappa))/(1-kappa)
TVaRX1

# Calcul de la TVaR avec prime SL

long<-length(fX1)-1
vj<-0:long
SLX1<-sum(fX1*pmax(vj-VaRX1,0))
TVaRX1b<-VaRX1+SLX1/(1-kappa)
TVaRX1b

c(kappa,VaRX1,TVaRX1,TVaRX1b)
