lmesim <-
function(N=6,minn=5,maxn=8,betas=c(1,.1),s2=1,D=matrix(c(1,.09,.09,.1),2,2),nExp=100,sr=2)
{
library(nlme)
set.seed(sr)
dump("lmesim","c:\\MixedModels\\Chapter03\\lmesim.r")
methods=c("ML","RML","MINQUE","MM","UVLS")
nmethods=length(methods)
ni=round(runif(N,min=minn,max=maxn))
NT=sum(ni)
TD=t(chol(D))
Dstar=s2*D
d=matrix(nrow=NT,ncol=1+1+2+2)
MSE=rep(0,nmethods)
names(MSE)=methods
for(iexp in 1:nExp)
{
  j <- 1
  for(i in 1:N)
  {
    nii <- ni[i]
    d[j:(j + nii - 1), 1] <- rep(i, nii)
    Xi <- cbind(rep(1, nii), 1:nii)
    d[j:(j + nii - 1), 3:4] <- Xi
    Zi <- Xi
    d[j:(j + nii - 1), 5:6] <- Zi
    bi=sqrt(s2)*TD%*%rnorm(2)
    yi=Xi%*%betas +Zi%*%bi+sqrt(s2)*rnorm(nii)
    d[j:(j + nii - 1), 2] <- yi
    j <- j + nii
  }
  
  o.lme=lmeFS(m=2,k=2,d,MLRML="ML")
  Dstar.lme <- o.lme$D*o.lme$s2
  MSE[1]=MSE[1]+sum((Dstar.lme-D)^2)/nExp

  o.lme=lmeFS(m=2,k=2,d,MLRML="RML")
  Dstar.lme <- o.lme$D*o.lme$s2
  MSE[2]=MSE[2]+sum((Dstar.lme-D)^2)/nExp
   
  o.MINQUE=lmevarMINQUE(m=2,k=2,d)
  Dstar.MINQUE=o.MINQUE[[2]]
  MSE[3]=MSE[3]+sum((Dstar.MINQUE-D)^2)/nExp

  Dstar.MM=lmevarMM(m=2,k=2,d,s2=o.MINQUE[[1]])[[2]]   
  MSE[4]=MSE[4]+sum((Dstar.MM-D)^2)/nExp

  Dstar.UVLS=lmevarUVLS(m=2,k=2,d)[[2]]   
  MSE[5]=MSE[5]+sum((Dstar.UVLS-D)^2)/nExp
}
return(cbind(MSE))
}
