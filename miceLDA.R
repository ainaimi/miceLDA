rm(list=ls());cat("\014")
source(file="~/Dropbox/Documents/R/ipak.R")
packages <- c("VIM","mice","cmprsk","ggplot2")
ipak(packages)
##
### Data generation
##
set.seed(123)
mc<-1000 # number of MC samples
n<-500 #Number of subjects
K<-1 # Number of causes of death
N<-5 # N number of intervals per subject

## This is the matrix of parameters of interest, possibly different
## at each interval
psi.mat<-matrix(0,nrow=K,ncol=N+1)

##Here are the effect sizes for the K=2 causes
psi.mat[1,]<- log(2.5)
#psi.mat[2,]<- log(2)

##Here the (untreated) all-cause rate is set to lambda=0.01, with
##lambda/K per cause; muK=lambda is used in the algorithm.

lambda<-0.075
gamma.vec<-rep(log(lambda/K))
muK<-sum(exp(gamma.vec))
R<-A<-L<-ID<-Y<-Z<-Tv<-Int<-ALast<-LLast<-RLast<-LFirst<-numeric()
T0.vec<-T.vec<-Y.vec<-Z.vec<-rep(0,n)

##Here are the coefficients determining the
##mediation and treatment assignment mechanisms.

bevec<-c(log(3/7),log(4),log(0.5),log(1.5))
alvec<-c(log(2/7),0.5,0.5,log(4))
rlvec<-c(log(1/12),log(1/4),0.5,0.5,log(4))

##cval is used as in Young's algorithm to introduce the confounding
cval<-18.5

##Begin the data-generation loop
a<-1;b<-n;c<-N
sim_dat<-numeric()
simfunc<-function(a,b,c){
      for(i in 1:b){
        ##Generate the counterfactual (untreated) survival time
        T0<-rexp(1,lambda)
        Ival<-as.numeric(T0<cval) # 
        ##Begin the interval-by-interval simulation
        m<-0
        mu.tot<-0
        R.vec<-A.vec<-L.vec<-ALast.vec<-RLast.vec<-LLast.vec<-LFirst.vec<-rep(0,c+1)
        ##Implement Young's algorithm with multiple causes
        ##Generate the survival time, then the cause
        while(muK*T0 > mu.tot & m <= c){
          if(m == 0){
            ##First interval
            eta<-bevec[1]+bevec[2]*Ival+bevec[3]*0+bevec[4]*0
            pval<-1/(1+exp(-eta))
            L.vec[m+1]<-rbinom(1,1,pval)
            eta<-alvec[1]+alvec[2]*L.vec[m+1]+alvec[3]*0+alvec[4]*0
            pval<-1/(1+exp(-eta))
            A.vec[m+1]<-rbinom(1,1,pval)
            
            eta<-rlvec[1]+rlvec[2]*Ival+rlvec[3]*L.vec[m+1]+rlvec[4]*0+rlvec[5]*0
            pval<-1/(1+exp(-eta))
            R.vec[m+1]<-rbinom(1,1,pval)
            
            ALast.vec[m+1]<-0;LLast.vec[m+1]<-0;RLast.vec[m+1]<-0
            LFirst.vec<-rep(L.vec[m+1],c+1)
          }else{
            ##Subsequent intervals
            eta<-bevec[1]+bevec[2]*Ival+bevec[3]*A.vec[m]+
              bevec[4]*L.vec[m]
            pval<-1/(1+exp(-eta))
            L.vec[m+1]<-rbinom(1,1,pval)
            eta<-alvec[1]+alvec[2]*L.vec[m+1]+alvec[3]*L.vec[m]+
              alvec[4]*A.vec[m]
            pval<-1/(1+exp(-eta))
            A.vec[m+1]<-rbinom(1,1,pval)
            
            eta<-rlvec[1]+rlvec[2]*Ival+rlvec[3]*L.vec[m+1]+rlvec[4]*L.vec[m]+
              rlvec[5]*A.vec[m]
            pval<-1/(1+exp(-eta))
            R.vec[m+1]<-rbinom(1,1,pval)
            
            ALast.vec[m+1]<-A.vec[m];LLast.vec[m+1]<-L.vec[m];RLast.vec[m+1]<-R.vec[m]
          }
          muval<-sum(exp(gamma.vec+A.vec[m+1]*psi.mat[,m+1]))
          ##Tval is computed for each interval, but is overwritten
          ##until the final interval
          Tval<-m+(muK*T0-mu.tot)/muval
          mu.tot<-mu.tot+muval
          m<-m+1
        }
        ##After exiting the loop, the survival time has been generated as Tval
        ##Now need to generate the failure type.
        if(m > c){
          ##In the case of censoring at tenth interval, no failure.
          Tval<-m-1
          Z.vec[i]<-0
        }else{
          ##In the case of failure, use the ratio hazards to define the
          ##relevant multinomial distribution on the K causes.
          Z.vec[i]<-sample(c(1:K),1,prob=exp(gamma.vec+A.vec[m]*psi.mat[,m]))
        }
        ##Store the outcomes
        T0.vec[i]<-T0
        T.vec[i]<-Tval
        Y.vec[i]<-m-1
        ID<-c(ID,rep(i,m))
        Int<-c(Int,c(1:m))
        A<-c(A,A.vec[1:m])
        L<-c(L,L.vec[1:m])
        R<-c(R,R.vec[1:m])
        ALast<-c(ALast,ALast.vec[1:m])
        LLast<-c(LLast,LLast.vec[1:m])
        RLast<-c(RLast,RLast.vec[1:m])
        LFirst<-c(LFirst,LFirst.vec[1:m])
        Z<-c(Z,rep(0,m-1),Z.vec[i])
        tv<-c(1:m);tv[m]<-Tval
        Tv<-c(Tv,tv)
      }
      D<-data.frame(a,ID,Int,Tv,A,ALast,L,LLast,R,RLast,Z)
      ##Trim off the intervals beyond the Nth (loop goes one too far)
      D<-D[D$Int<=c,]
      names(D)<-c("mc","id","Int","t","x","x1","z","z1","r","r1","y")
      sim_dat<-rbind(sim_dat,D)
      return(sim_dat)
}

# mc number of MC samples
# n Number of subjects
# N number of intervals per subject
t2<-lapply(1:mc,function(x) simfunc(x,n,N))

r_mean<-unlist(lapply(1:mc,function(x) mean(t2[[x]]$r)))
mean(r_mean)
str(t2)
head(t2[[12]])
mFunc<-function(a){
  a$x_m<-a$x
  a$x_m1<-a$x1
  a$x<-ifelse(a$r==1,NA,a$x)
  a$x1<-ifelse(a$r1==1,NA,a$x1)
  a$r<-a$r1<-NULL
  return(a)
}

t2<-lapply(1:mc,function(x) mFunc(t2[[x]]))
str(t2)
head(t2[[12]],20)
aggr(t2[[14]])

g2<-rgb(169/255,169/255,169/255,alpha=0.05)
plot(NULL)
for(i in 1:mc){
  par(new=T)
  plot(cuminc(t2[[i]]$t,t2[[i]]$y,cencode=0),col=g2)  
}

# MSM
# exposure model
pFunc<-function(a,b,c){
  aa<-predict(glm(b,data=c),type="response")
  aa<-aa*a+(1-aa)*(1-a)
  return(aa)
}

wght<-function(a){
  f_num<-as.formula(x~Int)
  f_den<-as.formula(x~Int+x1+z+z1)
  
  a$x_num<-pFunc(a$x,f_num,a)
  a$x_den<-pFunc(a$x,f_den,a)
  
  a$x_num<-ave(a$x_num,a$id,FUN=cumprod)
  a$x_den<-ave(a$x_den,a$id,FUN=cumprod)
  
  a$sw<-a$x_num/a$x_den
  
  f_num<-as.formula(x_m~Int)
  f_den<-as.formula(x_m~Int+x_m1+z+z1)
  
  a$x_num<-pFunc(a$x_m,f_num,a)
  a$x_den<-pFunc(a$x_m,f_den,a)
  
  a$x_num<-ave(a$x_num,a$id,FUN=cumprod)
  a$x_den<-ave(a$x_den,a$id,FUN=cumprod)
  
  a$sw_m<-a$x_num/a$x_den
  a$x_num<-NULL
  a$x_den<-NULL

  return(a)
}

q<-lapply(1:mc,function(x) wght(t2[[x]]))

ggplot(q[[10]], aes(x=as.factor(Int), y=log(sw))) +
  geom_boxplot() +
  geom_point(position = position_jitter(width = 0.2))

ggplot(q[[10]], aes(x=as.factor(Int), y=log(sw_m))) +
  geom_boxplot() +
  geom_point(position = position_jitter(width = 0.2))

psi_a<-psi_b<-psi_c<-psi_d<-numeric()
for(i in 1:mc){
  #TRUE
  m1<-coxph(Surv(Int-1,t,y)~x_m+cluster(id),data=q[[i]])
  psi_a<-rbind(psi_a,c(exp(coef(m1)),summary(m1)$coefficients[4]))
  m2<-coxph(Surv(Int-1,t,y)~x_m+cluster(id),weight=sw_m,data=q[[i]])
  psi_b<-rbind(psi_b,c(exp(coef(m2)),summary(m2)$coefficients[4]))
  #COMPLETE CASE
  m1<-coxph(Surv(Int-1,t,y)~x+cluster(id),data=q[[i]])
  psi_c<-rbind(psi_c,c(exp(coef(m1)),summary(m1)$coefficients[4]))
  m2<-coxph(Surv(Int-1,t,y)~x+cluster(id),weight=sw,data=q[[i]])
  psi_d<-rbind(psi_d,c(exp(coef(m2)),summary(m2)$coefficients[4]))
}
c(mean(psi_a[,1]),mean(psi_a[,2]),sd(log(psi_a[,1])))
c(mean(psi_b[,1]),mean(psi_b[,2]),sd(log(psi_b[,1])))
c(mean(psi_c[,1]),mean(psi_c[,2]),sd(log(psi_c[,1])))
c(mean(psi_d[,1]),mean(psi_d[,2]),sd(log(psi_d[,1])))






