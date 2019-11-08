# Lyon A2018
# Date : mardi 2018-11-06
# Somme de v.a. i.i.d.
# 
# Loi de X : Gamma
#
aa<-0.5
bb<-aa/10
EX<-aa/bb
EX
VarX<-EX/bb
# 
x<-0:100
n1<-1
f1<-dgamma(x,aa*n1,bb*n1)

n2<-10
f2<-dgamma(x,aa*n2,bb*n2)


n3<-100
f3<-dgamma(x,aa*n3,bb*n3)

n4<-1000
f4<-dgamma(x,aa*n4,bb*n4)

n3<-100

n4<-1000

f4<-dgamma(x,aa*n4,bb*n4)

matplot(x,cbind(f1,f2,f3,f4),type="l",xlab="x",ylab="f(x)")
