ssLOG <-
function(Prx=.4,Prz=.3,sigmaINT=.25,alpha=-1,beta=.2,gamma=.1,delta=log(2),nExp=300)
{
dump("ssLOG","c:\\MixedModels\\Chapter07\\ssLOG.r")
library(MASS)
NS=seq(from=50,to=100,by=10);LNS=length(NS)
ns=seq(from=100,to=200,by=20);Lns=length(ns)
power=as.data.frame(matrix(0,nrow=Lns,ncol=LNS))
names(power)=paste("N =",NS)
row.names(power)=paste("n =",ns)
print(date())
for(iN in 1:LNS)
for(inn in 1:Lns)
{
	N=NS[iN];n=ns[inn]
	NT=N*n
	d=as.data.frame(matrix(nrow=NT,ncol=5))
	names(d)=c("id","y","x","z","xz")
	pv=rep(NA,nExp)
    for(ix in 1:nExp)
    {
		k=1
		for(i in 1:N) 
	    {
			 xi<-rep(0,n);xi[runif(n)<Prx]=1
			 zi<-rep(0,n);zi[runif(n)<Prz]=1
			 d[k:(k + n - 1), 1] <- i
			 d[k:(k + n - 1), 3] <- xi
			 d[k:(k + n - 1), 4] <- zi
			 d[k:(k + n - 1), 5] <- xi*zi
		     yi=rep(0,n)
		     LEX=alpha+beta*xi+gamma*zi+delta*(xi*zi)+rnorm(1,mean=0,sd=sigmaINT)
			 yi[rnorm(n)<1/(1+exp(-LEX))]=1 
		     d[k:(k + n - 1), 2] <- yi
			 k <- k + n
		}
		o <- glmmPQL(y~x+z+xz, random=~1|id,family=binomial,data=d,verbose=F)
		if((attr(o,"class"))[1]!="try-error")  pv[ix]=(summary(o)[[21]])[4,5]	
    } #ix loop is over
    pv=pv[!is.na(pv)]
    power[iN,inn]=mean(pv<.05)
print(c(iN,inn))
}
return(power)
}
