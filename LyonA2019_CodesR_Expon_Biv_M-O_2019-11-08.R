# Cours : Lyon
# Semestre : A2019
#
# Loi exponentielle bivariée Marshall-Olkin
#
be1<-0.1
be2<-0.2
a0<-0.05

set.seed(20191108)

mm<-1000
Y0<-rexp(mm,a0)
Y1<-rexp(mm,be1-a0)
Y2<-rexp(mm,be2-a0)

X1<-pmin(Y1,Y0)
X2<-pmin(Y2,Y0)

plot(X1,X2, xlab="réalisations de X1",ylab="réalisations de X2")

plot(pexp(X1,be1),pexp(X2,be2),xlab="réalisations de U1",ylab="réalisations de U2")




