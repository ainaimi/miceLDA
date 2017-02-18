rm(list=ls());cat("\014")
source(file="~/Dropbox/Documents/R/ipak.R")
packages <- c("VIM","mice")
ipak(packages)
##
### Data generation
##
set.seed(123)
K<-2 # Number of causes of death

## This is the matrix of parameters of interest, possibly different
## at each interval
psi.mat<-matrix(0,nrow=K,ncol=N+1)

##Here are the effect sizes for the K=2 causes
psi.mat[1,]<- -log(2.5)
psi.mat[2,]<- log(4.5)

##Here the (untreated) all-cause rate is set to lambda=0.01, with
##lambda/K per cause; muK=lambda is used in the algorithm.

lambda<-0.075
gamma.vec<-rep(log(lambda/K))
muK<-sum(exp(gamma.vec))
A<-L<-ID<-Y<-Z<-Tv<-Int<-ALast<-LLast<-LFirst<-numeric()
T0.vec<-T.vec<-Y.vec<-Z.vec<-rep(0,n)

##Here are the coefficients determining the
##mediation and treatment assignment mechanisms.

bevec<-c(log(3/7),log(4),log(0.5),log(1.5))
alvec<-c(log(2/7),0.5,0.5,log(4))

##cval is used as in Young's algorithm to introduce the confounding

cval<-80

##Begin the data-generation loop
sim_dat<-numeric()
simfunc<-function(a,b,c){
      for(i in 1:b){
        ##Generate the counterfactual (untreated) survival time
        T0<-rexp(1,lambda)
        Ival<-as.numeric(T0 < cval)
        ##Begin the interval-by-interval simulation
        m<-0
        mu.tot<-0
        A.vec<-L.vec<-ALast.vec<-LLast.vec<-LFirst.vec<-rep(0,c+1)
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
            ALast.vec[m+1]<-0;LLast.vec[m+1]<-0
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
            ALast.vec[m+1]<-A.vec[m];LLast.vec[m+1]<-L.vec[m]
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
        ALast<-c(ALast,ALast.vec[1:m])
        LLast<-c(LLast,LLast.vec[1:m])
        LFirst<-c(LFirst,LFirst.vec[1:m])
        Z<-c(Z,rep(0,m-1),Z.vec[i])
        tv<-c(1:m);tv[m]<-Tval
        Tv<-c(Tv,tv)
      }
      D<-data.frame(a,ID,Int,Tv,A,ALast,L,LLast,Z)
      ##Trim off the intervals beyond the Nth (loop goes one too far)
      D<-D[D$Int<=c,]
      names(D)<-c("mc","id","Int","t","x","x1","z","z1","y")
      sim_dat<-rbind(sim_dat,D)
      return(sim_dat)
}

mc<-10 # number of MC samples
n<-500 #Number of subjects
N<-5 #number of intervals per subject
t2<-lapply(1:mc,function(x) simfunc(x,n,N))
str(t2)

head(sim_dat,50)
tail(sim_dat,50)






