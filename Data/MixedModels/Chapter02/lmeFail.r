lmeFail <-
function(Nexp=5000,method="REML")
{
 dump("lmeFail","c:\\MixedModels\\Chapter02\\lmeFail.r")
 library(nlme);library(lme4)
 k=m=2
 alpha=-1;beta=.2;sigma2eps=.1;sigma2a=sigma2b=.1;ro=.5
 Dstar=matrix(c(sigma2a,sqrt(sigma2a*sigma2b)*ro,sqrt(sigma2a*sigma2b)*ro,sigma2b),k,k)
 TD=t(chol(Dstar))
 NS=seq(from=10,to=50,by=5)
 nNS=length(NS)
 lmeF=lmeF4=be=rep(0,nNS)

 for(iN in 1:nNS)
 {
   N=NS[iN]
   ni=round(runif(N,min=3,max=7))	
   NT=sum(ni)
   d=matrix(ncol=3,nrow=NT)
   for(iexp in 1:Nexp)
   {
     j <- 1
     for(i in 1:N) {
        n <- ni[i]
        d[j:(j + n - 1), 1] <- i
        d[j:(j + n - 1), 2] <- 1:n
        b2=TD%*%rnorm(2)
        d[j:(j+n-1),3] <- (alpha+b2[1])+(beta+b2[2])*(1:n)+rnorm(n,mean=0,sd=sqrt(sigma2eps))
        j <- j + n
    }
    dL <- as.data.frame(d)
    names(dL)<-c("id","x","y")
    out<-try(lme(fixed=y~x,random=~x|id,data=dL,method=method))
    if(attr(out,"class")=="try-error") lmeF[iN]=lmeF[iN]+1/Nexp

    out<-try(lmer(y~x+(x|id),data=dL))
    if(attr(out,"class")=="try-error") lmeF4[iN]=lmeF4[iN]+1/Nexp
	print(summary(out))
print(NT)
dmy=cbind(dL$id,dL$y,rep(1,nrow(d)),dL$x,rep(1,nrow(d)),dL$x)
print(dmy)
o.lme=lmeFS(m=2,k=2,dmy,MLRML="RML",pr=1)
  Dstar.lme <- o.lme$D*o.lme$s2
print(Dstar)


 
  }
 }
cbind(NS,lmeF,lmeF4)
}
