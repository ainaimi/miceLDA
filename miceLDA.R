rm(list=ls());cat("\014")
source(file="~/Dropbox/Documents/R/ipak.R")
packages <- c("VIM","mice","cmprsk","ggplot2",
              "data.table","parallel","microbenchmark",
              "profvis")
ipak(packages)
##
### Data generation
##
time1<-proc.time()
set.seed(123)
mc<-500 # number of MC samples
n<-250 #Number of subjects
K<-1 # Number of causes of death
N<-5 # N number of intervals per subject
iterations<-50 # number of iterations for MICE
imputations<-10 # number of imputations for MICE
cores <- detectCores()-1 # for all parallel processing

## This is the matrix of parameters of interest, possibly different
## at each interval
psi.mat<-matrix(0,nrow=K,ncol=N+1)

##Effect size 
psi.mat[1,]<- log(2.5)

##Here are the coefficients determining the
##mediation and treatment assignment mechanisms.
cvec<-c(-1*.75,.25)
  b.Int<-(log((1/(1-.25))-1)-log(4)*0.75-log(0.5)*.5-log(1.5)*.5-log(2)*0)
bevec<-c(b.Int,log(4),log(0.5),log(1.5),log(2))
  a.Int<-(log((1/(1-.25))-1)-.5*.5-.5*.5-log(4)*.5-log(2)*0)
alvec<-c(a.Int,.5,.5,log(4),log(2))
  r.Int<-(log((1/(1-.15))-1)-log(1/4)*.75-.5*.5-.5*.5-log(4)*.5-log(2)*0)
rlvec<-c(r.Int,log(1/4),0.5,0.5,log(4),log(2))

##Here the (untreated) all-cause rate is set to lambda=0.01, with
##lambda/K per cause; muK=lambda is used in the algorithm.

lambda<-0.075
gamma.vec<-rep(log(lambda/K))
muK<-sum(exp(gamma.vec))
C<-R<-A<-L<-ID<-Y<-Z<-Tv<-Int<-ALast<-LLast<-RLast<-LFirst<-numeric()
T0.vec<-C.vec<-T.vec<-Y.vec<-Z.vec<-rep(0,n)

##cval is used as in Young's algorithm to introduce the confounding
cval<--log(1-.75)/lambda

quickdf <- function(l) {
  class(l) <- "data.frame"
  attr(l, "row.names") <- .set_row_names(length(l[[1]]))
  l
}

##Begin the data-generation loop
a<-1;b<-n;c<-N
sim_dat<-numeric()
simfunc<-function(a,b,c){
      #for(i in 1:b){
      loopr<-function(i,a,c){
        ##Generate the counterfactual (untreated) survival time
        T0<-rexp(1,lambda)
        Ival<-as.numeric(T0<cval)
        ##Continuous confounder
        eta<-cvec[1]+cvec[2]*Ival
        C<-rnorm(1,eta,1)
        ##Begin the interval-by-interval simulation
        m<-0
        mu.tot<-0
        R.vec<-A.vec<-L.vec<-ALast.vec<-RLast.vec<-LLast.vec<-LFirst.vec<-rep(0,c+1)
        ##Implement Young's algorithm with multiple causes
        ##Generate the survival time, then the cause
        while(muK*T0 > mu.tot & m <= c){
          if(m == 0){
            ##First interval
            eta<-bevec[1]+bevec[2]*Ival+bevec[3]*0+bevec[4]*0+bevec[5]*C
            pval<-1/(1+exp(-eta))
            L.vec[m+1]<-rbinom(1,1,pval)
            eta<-alvec[1]+alvec[2]*L.vec[m+1]+alvec[3]*0+alvec[4]*0+alvec[5]*C
            pval<-1/(1+exp(-eta))
            A.vec[m+1]<-rbinom(1,1,pval)
            eta<-rlvec[1]+rlvec[2]*Ival+rlvec[3]*L.vec[m+1]+rlvec[4]*0+rlvec[5]*0+rlvec[6]*C
            pval<-1/(1+exp(-eta))
            R.vec[m+1]<-rbinom(1,1,pval)
            ALast.vec[m+1]<-0;LLast.vec[m+1]<-0;RLast.vec[m+1]<-0
            LFirst.vec<-rep(L.vec[m+1],c+1)
          }else{
            ##Subsequent intervals
            eta<-bevec[1]+bevec[2]*Ival+bevec[3]*A.vec[m]+
              bevec[4]*L.vec[m]+bevec[5]*C
            pval<-1/(1+exp(-eta))
            L.vec[m+1]<-rbinom(1,1,pval)
            eta<-alvec[1]+alvec[2]*L.vec[m+1]+alvec[3]*L.vec[m]+
              alvec[4]*A.vec[m]+alvec[5]*C
            pval<-1/(1+exp(-eta))
            A.vec[m+1]<-rbinom(1,1,pval)
            eta<-rlvec[1]+rlvec[2]*Ival+rlvec[3]*L.vec[m+1]+rlvec[4]*L.vec[m]+
              rlvec[5]*A.vec[m]+rlvec[6]*C
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
        C.vec[i]<-C
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
        SD<-quickdf(list(rep(a,m),ID,Int,Tv,rep(C,m),A,ALast,L,LLast,R,RLast,Z)) # using data.table is faster than data.frame
        return(SD)
      }
      SD<-lapply(1:b,function(x) loopr(x,a,c))
      D<-rbindlist(SD)
      #D<-do.call(rbind,SD);D<-data.table(D)
      ##Trim off the intervals beyond the Nth (loop goes one too far)
      names(D)<-c("mc","id","Int","t","c","x","x1","z","z1","r","r1","y")
      D<-D[Int<=5]
      return(D)
}

# mc number of MC samples
# n Number of subjects
# N number of intervals per subject
t2<-mclapply(1:mc,function(x) simfunc(x,n,N),mc.cores=cores)

# Replace observed with NA for those with R=1
# r_mean<-unlist(lapply(1:mc,function(x) mean(t2[[x]]$r)))
# mean(r_mean)
# str(t2)
# head(t2[[12]])
mFunc<-function(a){
  a$x_m<-a$x
  a$x_m1<-a$x1
  a$x<-ifelse(a$r==1,NA,a$x)
  a$x1<-ifelse(a$r1==1,NA,a$x1)
  a$Int<-as.factor(a$Int)
  a$r<-a$r1<-NULL
  return(a)
}

t2<-mclapply(1:mc,function(x) mFunc(t2[[x]]),mc.cores=cores)

# g2<-rgb(169/255,169/255,169/255,alpha=0.05)
# plot(NULL)
# for(i in 1:mc){
#   par(new=T)
#   plot(cuminc(t2[[i]]$t,t2[[i]]$y,cencode=0),col=g2,xlim=c(0,5),ylim=c(0,.3))
# }

# MSM
# exposure model
pFunc<-function(a,b,c){
  aa<-predict(glm(b,data=c,family=binomial("logit")),type="response")
  aa<-aa*a+(1-aa)*(1-a)
  return(aa)
}

wght<-function(a){
  f_num<-as.formula(x~as.factor(Int))
  f_den<-as.formula(x~as.factor(Int)+x1+z+z1+c)

  a$x_num<-ifelse(is.na(a$x),1,pFunc(a$x,f_num,a))
  a$x_den<-ifelse(is.na(a$x),1,pFunc(a$x,f_den,a))

  a$x_num<-ave(a$x_num,a$id,FUN=cumprod)
  a$x_den<-ave(a$x_den,a$id,FUN=cumprod)

  a$sw<-a$x_num/a$x_den
  
  f_num<-as.formula(x_m~as.factor(Int))
  f_den<-as.formula(x_m~as.factor(Int)+x_m1+z+z1+c)
  
  a$x_num<-pFunc(a$x_m,f_num,a)
  a$x_den<-pFunc(a$x_m,f_den,a)
  
  a$x_num<-ave(a$x_num,a$id,FUN=cumprod)
  a$x_den<-ave(a$x_den,a$id,FUN=cumprod)
  
  a$sw_m<-a$x_num/a$x_den
  a$x_num<-NULL
  a$x_den<-NULL
  return(a)
}

q<-mclapply(1:mc,function(x) wght(t2[[x]]),mc.cores=cores)
head(q[[2]])

# ggplot(q[[10]], aes(x=as.factor(Int), y=log(sw))) +
#   geom_boxplot() +
#   geom_point(position = position_jitter(width = 0.2))
# 
# ggplot(do.call(rbind,q), aes(x=as.factor(Int), y=log(sw_m))) +
#   geom_boxplot() +
#   geom_point(position = position_jitter(width = 0.2),color=g2)

psi_a<-psi_b<-psi_c<-psi_d<-numeric()
eFunc<-function(a,b,c){
  m1<-coxph(Surv(as.numeric(Int)-1,t,y)~a+cluster(id),weights=b,data=c)
  p<-c(exp(coef(m1)),summary(m1)$coefficients[4])
  return(p)
}

psi_a<-do.call(rbind,(mclapply(1:mc,function(x) eFunc(q[[x]]$x_m,NULL,q[[x]]),mc.cores=cores)))
psi_b<-do.call(rbind,(mclapply(1:mc,function(x) eFunc(q[[x]]$x_m,q[[x]]$sw_m,q[[x]]),mc.cores=cores)))
psi_c<-do.call(rbind,(mclapply(1:mc,function(x) eFunc(q[[x]]$x,NULL,q[[x]]),mc.cores=cores)))
psi_d<-do.call(rbind,(mclapply(1:mc,function(x) eFunc(q[[x]]$x,q[[x]]$sw,q[[x]]),mc.cores=cores)))

c(mean(psi_a[,1]),mean(psi_a[,2]),sd(log(psi_a[,1])))
c(mean(psi_b[,1]),mean(psi_b[,2]),sd(log(psi_b[,1])))
c(mean(psi_c[,1]),mean(psi_c[,2]),sd(log(psi_c[,1])))
c(mean(psi_d[,1]),mean(psi_d[,2]),sd(log(psi_d[,1])))

##multiple imputation
head(q[[10]])
fFunc<-function(a){
  a$x1<-a$x_m<-a$x_m1<-a$sw<-a$sw_m<-NULL
  a$x<-as.factor(a$x)
  a$Int<-as.factor(a$Int)
  return(a)
}
q<-mclapply(1:mc,function(x) fFunc(q[[x]]),mc.cores=cores)
str(q[[10]])

ini <- mice(q[[10]],seed=123,maxit=0)
pMatrix<-ini$predictorMatrix
pMatrix[6,]<-c(0,0,1,0,1,0,1,1,1)
pMatrix

methd<-c("","","","","","logreg","","","")

# cores <- detectCores()-1
# cl <- makeCluster(cores)
# clusterCall(cl, function() library("mice"))
# clusterExport(cl, c("methd","pMatrix","iterations","imputations"))
imp<-mclapply(q,function(x) mice(x,seed=123,maxit=iterations,m=imputations,
                               method=methd,
                               predictorMatrix=pMatrix,
                               diagnostics=F),mc.cores=cores)

# alph<-.02
# rr<-rgb(1,0,0,alpha=alph)
# gg<-rgb(0,205/255,0,alpha=alph)
# bb<-rgb(0,0,1,alpha=alph)
# cc<-rgb(0,1,1,alpha=alph)
# mm<-rgb(1,0,1,alpha=alph)
# pal=c(rr,gg,bb,cc,mm)
# plot(1:iterations,imp[[1]]$chainMean[1,,1]-1,col=rgb(1,0,0,alpha=0),
#      type="l",ylim=c(0,.4),xlab="Iteration",ylab="Mean of X",
#      las=1,tcl=-.1,yaxs="i",xaxs="i")
# for(i in 1:mc){
#   par(new=T)
#   plot(1:iterations,imp[[i]]$chainMean[1,,1]-1,type="l",ylim=c(0,.4),col=pal[1],ylab="",xlab="",xaxt="n",yaxt="n",yaxs="i",xaxs="i")
#   lapply(2:imputations,function(x) lines(imp[[i]]$chainMean[1,,x]-1,col=pal[x]))
# }

imp_dat<-mclapply(1:mc,function(x) complete(imp[[x]],"long"),mc.cores = cores)
nu<-function(a){
  # convert from factor to numeric
  a$x<-as.numeric(a$x)-1
  a$Int<-as.numeric(a$Int)-1
  # lag the imputed exposure
  a<-data.table(a)
  a[, x1:=c(0, x[-.N]), by=id]
  #a$i<-a$i2<-a$i3<-NULL
  return(a)
  }
impD<-mclapply(1:mc,function(a) nu(imp_dat[[a]]),mc.cores = cores)
mean(t2[[1]]$x,na.rm=T)
mean(impD[[1]]$x)

pFunc<-function(a){
  for(i in 1:imputations){
    a$num<-predict(glm(x~as.factor(Int)+cluster(id),
                       family=binomial("logit"),
                       data=a,subset=.imp==i),type="response")
    a$num<-a$num*a$x+(1-a$num)*(1-a$x)

    a$den<-predict(glm(x~as.factor(Int)+x1+z+z1+c+cluster(id),
                       family=binomial("logit"),
                       data=a,subset=.imp==i),type="response")
    a$den<-a$den*a$x+(1-a$den)*(1-a$x)
    
    a$sw<-ave(a$num/a$den,a$id,FUN=cumprod)
    
    a$sw<-ifelse(a$sw>quantile(a$sw,.9),quantile(a$sw,.9),a$sw)
    
    a$num<-a$den<-NULL
    return(a)
  }
}
t<-mclapply(1:mc,function(x) pFunc(impD[[x]]),mc.cores = cores)

hist(do.call(rbind,lapply(1:mc,function(x) mean(t[[x]]$sw))))

psi_e<-p<-numeric()
eFunc<-function(a,b){
  for(i in 1:imputations){
  m1<-coxph(Surv(Int,t,y)~x+cluster(id),weights=sw,data=a,subset=.imp==i)
  p<-rbind(p,cbind(b,i,coef(m1),summary(m1)$coefficients[4]))
  }
  return(p)
}

psi_e<-mclapply(1:mc,function(x,y) eFunc(t[[x]],x),mc.cores=cores)
# p<-data.frame(do.call(rbind,psi_e))
# head(psi_e,50)
# names(p)<-c("mc","imp","psi","se")

RubinsRules<-function(estimate,se){
  q<-exp(mean(estimate))
  u<-mean(se)
  b<-var(estimate)
  m<-nrow(as.matrix(estimate))
  t<-u+(1+(1/m))*b
  se<-sqrt(t)
  results<-c(q,se,u,b)
  names(results)<-c("psi","se","se_bar","sd_psi")
  return(results)
}
res<-mclapply(1:mc,function(a) RubinsRules(psi_e[[a]][,3],psi_e[[a]][,4]),mc.cores=cores)
psi_e<-do.call(rbind,res)

do.call(rbind,lapply(1:mc,function(x) mean(t[[x]]$sw)))
abline(v=1,col="red",lwd=2)

res<-data.frame(rbind(
  c("true unadjusted",mean(psi_a[,1]),mean(psi_a[,2]),sd(log(psi_a[,1]))),
  c("true weighted",mean(psi_b[,1]),mean(psi_b[,2]),sd(log(psi_b[,1]))),
  c("cc unadjusted",mean(psi_c[,1]),mean(psi_c[,2]),sd(log(psi_c[,1]))),
  c("cc weighted",mean(psi_d[,1]),mean(psi_d[,2]),sd(log(psi_d[,1]))),
  c("mice weighted",mean(psi_e[,1]),mean(psi_e[,2]),sd(log(psi_e[,1])))
))
names(res)<-c("type","psi","se","se_mc")
res
(proc.time()-time1)/60

## to do: 
##        
##        3) modify code to get results under more 
##          missing data mechanisms
##        4) modify code to increase MC sample, 
##          imputatino and iteration number
##        5) look into splitting up the data/imputation in each MC 
##          sample into pieces

## https://www.ncbi.nlm.nih.gov/pubmed/28034175

