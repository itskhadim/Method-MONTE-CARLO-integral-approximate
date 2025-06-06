#PROJET
library("rgl")
library("pracma")
n=1000
phi= function(x,y) 2*x*y*sqrt(x^2+y^2)
x=seq(0,2,length=20)
y=seq(0,1,length=20)
z=outer(x,y,phi)
persp3d(x,y,z,col="lightblue")

#Valeur de l'integrale
I=quad2d(phi,xa=0,xb=2,ya=0,yb=1)
I
#1ER METHODE (FREQUENCE EMPIRIQUE)
set.seed(1)
U=runif(n,0,2)
V=runif(n,0,1)
K=max(phi(U,V))
Z=runif(n,0,K)
plot(U,V,main='nuage de points')
I_meth_1=K*2*(1/n)*sum(Z<phi(U,V))
I_meth_1

#Illustration de LFGN
meth_1_evol=K*1*2*cumsum(Z<phi(U,V))/(1:n)
plot(meth_1_evol,type='l',col="red",main='Illustration LFGN Meth1')
abline(h=I,col="blue")
legend("topright", legend=c("meth_1_evol", "h=I"),
       col=c("red", "blue"), lty=c(1,1), cex=0.8)

#Illustration du TCL 
set.seed(1)
N=200
U=runif(n*N,0,2)
V=runif(n*N,0,1)
calcul_meth_1 = function(X) 2*K*mean(X>Z)
M = matrix(U^2, ncol=N)
meth_1_vect =apply(M,MARGIN=1,FUN=calcul_meth_1)
hist(meth_1_vect,col="lightblue",main = "Illustration du TCL Meth1")
abline(v=I,col="red",lty=2,lwd=5)
legend("topright", legend="h=I",
       col="red", lty=c(1,1), cex=0.9)


#2e METHODE (2 UNIFORMES independantes)
set.seed(2)
U=runif(n,0,2)
V=runif(n,0,1)
I_meth_2=2*mean(phi(U,V))
I_meth_2

#Illustration de LFGN
meth_2_evol=2*cumsum(phi(U,V))/(1:n)
plot(meth_2_evol,type='l',col="red",main='Illustration LFGN Meth2')
abline(h=I,col="blue")
legend(450,1.5,c("meth_2_evol", "h=I"), col=c("red", "blue"), lty = c(1,1))

#Illustration du TCL
U=runif(n*N,0,2)
V=runif(n*N,0,1)
calcul_meth_2 = function(X,Y) 1*2*mean(phi(X,Y))
M = matrix(U, ncol=N)
L = matrix(V, ncol=N)
meth_2_vect=apply(M,L,MARGIN=2,FUN=calcul_meth_2)
hist(meth_2_vect,col="lightblue",main = "Illustration du TCL Meth2")
abline(v=I,col="red",lty=2,lwd=5)
legend("topright", legend="h=I",
       col="red", lty=c(1,1), cex=0.9)


#3e methode (p1(x) et p2(y))
set.seed(2)
U=runif(n,0,2)
V=runif(n,0,1)
P1inv=function(u) (728*u+1)^(1/6)-1
P2inv=function(u) u^(1/3)
X=P1inv(runif(n,0,1))
hist(X,col="orange") 
Y=P2inv(runif(n,0,1))
hist(Y,col="lightblue")
p1=function(x) (6/728)*(1+x)^5
p2=function(y) 3*y^2
p=function(x,y) p1(x)*p2(y)
I_meth_3=mean(phi(X,Y)/p(X,Y))
I_meth_3

#Illustration de LFGN
meth_3_evol=cumsum(phi(X,Y)/p(X,Y))/(1:n)
plot(meth_3_evol,type='l',col="red",main='Illustration LFGN Meth3')
abline(h=I,col="blue")
legend(440,2.6,c("meth_3_evol", "h=I"), col=c("red", "blue"), lty = c(1,1))

#Illustration du TCL
set.seed(2)
U=runif(n*N,0,2)
V=runif(n*N,0,1)
calcul_meth_3 = function(X) 2*mean(X)
M = matrix(phi(U,V), ncol=N)
meth_3_vect=apply(M,MARGIN=2,FUN=calcul_meth_3)
hist(meth_3_vect,col="lightblue",main = "Illustration du TCL Meth3")
abline(v=I,col="red",lty=2,lwd=5)
legend("topright", legend="h=I",
       col="red", lty=c(1,1), cex=0.9)


#Methode 4 (Avec les copules)
#Copule de claton
# Simualtion de la copule de Clayton 
set.seed(111)
cop_Clayton=function(x,y,t=1.143)(x^(-t)+y^(-t)-1)^(-1/t)
u=seq(0,1,length=20)
v=seq(0,1,length=20)
phi=function(x,y) 2*x*y*(sqrt(x**2+y**2))
z2=outer(u,v,phi)
z_cop=outer(u,v,cop_Clayton)
rgl::persp3d(u,v,z_cop,col="lightgreen")
rgl::persp3d(u,v,z2,col="lightblue")

#Simulation de la densit?  copule de Clayton avec les vecteurs alé???atoires
rcop_clayton_final=function(n,t=1.143){
  Fcond_inv=function(u,z,t=1.143){((z^(-t/(t+1))-1)*u^(-t)+1)^(-1/t)}
  U1=runif(n,0,1) 
  Z=runif(n,0,1)
  V1=1-Fcond_inv(U1,Z)
  return (c(1-U1,V1))
}

#Fct densite p(x,y) -- plus optimale
ph<-function(u,v){
  return(dunif(u)*dunif(v)*cop_Clayton(punif(u,0,2),punif(v,0,1),t=1.143))
} 

#Ici mettre la densit? de la copule avec le changement de variable 
d_clayton2 = function(u,v,t=1.143){
  y1 = (1+t)*((1-u)*(1-v))^(-t-1)
  y2 = (1-u)^(-t) + (1-v)^(-t) - 1
  y3 = y2^(-(1/t) - 2)
  return(y1*y3)
}

ph<-function(u,v){
  return(dunif(u,0,2)*dunif(v,0,1)*d_clayton2(punif(u,0,2),punif(v,0,1),t=1.143))
} #Ici mettre la densit? de copule avec le changement de variable 


x=seq(0,2,length=20)
y=seq(0,1,length=20)
z = outer(x,y,ph)
rgl::persp3d(x,y,z,col="lightblue")

toto=rcop_clayton_final(n)
U2=toto[1:n]
V2=toto[(n+1):(2*n)]
plot(U2,V2)

#pour avoir le nuage de point sur [0,2]x[0,1]
toto=rcop_clayton_final(n)
U2=qunif(toto[1:n], 0,2)
V2=qunif(toto[(n+1):(2*n)], 0,1)
plot(U2,V2,main='nuage de points (U2,V2)')


#Illustration par histogramme
hist(U2,col="orange") 
hist(V2,col="lightblue")

I_meth_4=mean(phi(U2,V2)/ph(U2,V2))
I_meth_4

#Illustration de LFGN
meth_4_evol=cumsum(phi(U2,V2)/ph(U2,V2))/(1:n)
plot(meth_4_evol,type='l',col="red",main='Illustration LFGN Meth4')
abline(h=I,col="blue")
legend(440,4.1,c("meth_4_evol", "h=I"), col=c("red", "blue"), lty = c(1,1))

#Illustration du TCL
U2=runif(n*N,0,2)
V2=runif(n*N,0,1)
calcul_meth_4 = function(X) 2*mean(X)
M = matrix(phi(U2,V2), ncol=N)
meth_4_vect=apply(M,MARGIN=2,FUN=calcul_meth_4)
hist(meth_4_vect,col="lightblue",main = "Illustration du TCL Meth4")
abline(v=I,col="red",lty=2,lwd=5)
legend("topright", legend="h=I",
       col="red", lty=c(1,1), cex=0.9)

#Approximation selon les diff?rentes valeurs de teta
t=c(1.143,2,3,4,5)
tab=matrix(nrow=length(t),ncol = n)
for (i in 1:length(t)){
  toto=rcop_clayton_final(n,t)
  U2=qunif(toto[1:n], 0,2)
  V2=qunif(toto[(n+1):(2*n)], 0,1)
  tab[i,]=cumsum(phi(U2,V2)/ph(U2,V2))/(1:n)
}
plot(tab[1,],type='l',col=1,ylim=c(-1,10),main = "Illustration MC4 LFGN selon diff?rents theta")
for (i in 1:length(t)-1){
  lines(tab[i+1,],type='l',col=i+1)
  abline(h=I,col="blue")
}
legend(620,11,
       legend=c('ta=1.143','tta=2','tta=3','tta=5','I'),
       col=c("black","red","green","lightblue","orange")
       ,lty=c(1,1,1))


#Evolutions des quatre methodes
plot(meth_1_evol,type='l',col="red",main='Evolution des 4 methodes')
lines(meth_2_evol,col="blue")
lines(meth_3_evol,col="black")
lines(meth_4_evol,col="orange")
abline(h=I, col="yellow")
legend(550,2.5, legend=c("meth_1_evol","meth_2_evol" ,"meth_3_evol","meth_4_evol","h=I"),
       col=c("red","blue","black","orange","yellow"), lty=c(1,1,1,1), cex=0.8)

#Moyennes et Variances
mean(meth_1_vect)
mean(meth_2_vect)
mean(meth_3_evol)
mean(meth_4_vect)
sd(meth_1_vect)
sd(meth_2_vect)
sd(meth_3_vect)
sd(meth_4_vect)
#FIN
