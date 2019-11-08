# Lyon A2018
# Date : mardi 2018-11-06
# Loi Poisson vs loi binomiale négative

lam<-5
EM<-lam

r1<-1
q1<-r1/(r1+EM)
r1*(1-q1)/q1


r2<-2
q2<-r2/(r2+EM)
r2*(1-q2)/q2

r3<-0.5
q3<-r3/(r3+EM)
r3*(1-q3)/q3

r4<-5
q4<-r4/(r4+EM)
r4*(1-q4)/q4


k<-0:20
pp<-dpois(k,lam)
p1<-dnbinom(k,r1,q1)
p2<-dnbinom(k,r2,q2)
p3<-dnbinom(k,r3,q3)
p4<-dnbinom(k,r4,q4)

matplot(k,cbind(pp,p1,p2,p3,p4),type="h",xlab="k",ylab="Pr(M=k)")

round(cbind(k,pp,p1,p2,p3,p4),4)

k<-0:20
Fp<-ppois(k,lam)
F1<-pnbinom(k,r1,q1)
F2<-pnbinom(k,r2,q2)
F3<-pnbinom(k,r3,q3)
F4<-pnbinom(k,r4,q4)

matplot(k,cbind(Fp,F1,F2,F3,F4),type="s",xlab="k",ylab="Pr(M=k)")

round(cbind(k,Fp,F1,F2,F3,F4),4)




