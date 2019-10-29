heightlog.nlme <-
function(sex=0)
{
 dump("heightlog.nlme","c:\\MixedModels\\Chapter08\\heightlog.nlme.r")
 library(nlme)
 gender=c("Boy","Girl")
 if(sex==1) a0=c(165,-2.2,-0.42,0.039) else a0=c(182,-1,-0.16,0.018) #start
 da=height.dat[height.dat$sex==sex,] #extract gender set
 cat("\nQLogist model for height of",gender[sex+1],"\n\n")

 cat("\n\n============== All four parameters are random and do not correlate\n")
 outnlme<-nlme(height~QLogist(a1,a2,a3,a4,x=year),fixed=a1+a2+a3+a4~1,random=pdDiag(a1+a2+a3+a4~1),groups=~id, data=da,start=c(a1=a0[1],a2=a0[2],a3=a0[3],a4=a0[4]))
 print(summary(outnlme))
 a=outnlme$coefficients$fixed
 cat("\nAcceleration time =",round((sqrt(2*a[4])-a[3])/(2*a[4]),1))

 cat("\n\n============== Parameters a1,a2, and a4 are random with unrestricted cov matrix\n")
 outnlme<-nlme(height~QLogist(a1,a2,a3,a4,x=year),fixed=a1+a2+a3+a4~1,random=a1+a2+a4~1,groups=~id, data=da,start=c(a1=a0[1],a2=a0[2],a3=a0[3],a4=a0[4]))
 print(summary(outnlme))
 a=outnlme$coefficients$fixed
 cat("\nAcceleration time =",round((sqrt(2*a[4])-a[3])/(2*a[4]),1))
 k=3 #the number of random effects
 os <- matrix(as.numeric(VarCorr(outnlme)),ncol=k+1)
 sdb=diag(os[1:k,2],k,k)
 LR=os[2:k,3:(k+1)]
 Db=matrix(ncol=k,nrow=k)
 R=diag(rep(1,k),ncol=k,nrow=k)
 for(i in 2:k)
   R[i,1:(i-1)]=R[1:(i-1),i]=LR[(i-1),1:(i-1)]
 Db=sdb%*%R%*%sdb
 cat("\nMatrix D-star=cov(a1,a2,a4):\n")
 print(Db)

 cat("\n\n============== Only adult height is subject-specific, a1=random\n")
 outnlme<-nlme(height~QLogist(a1,a2,a3,a4,x=year),fixed=a1+a2+a3+a4~1,random=a1~1,groups=~id, data=da,start=c(a1=a0[1],a2=a0[2],a3=a0[3],a4=a0[4]))
 print(summary(outnlme))
 a=outnlme$coefficients$fixed
 cat("\nAcceleration time =",round((sqrt(2*a[4])-a[3])/(2*a[4]),1))
 
 par(mfrow=c(7,10),mar=c(2,2,1,1))
 a1.id=outnlme$coefficients$random$id
 xx=seq(from=7,to=18,by=.1)
 id=da$id;uid=unique(id);nud=length(uid)
 for(i in 1:nud)
 {
     x=da$year[id==uid[i]];y=da$height[id==uid[i]]
     y=y[order(x)];x=x[order(x)] # order to make lines
     plot(x,y,type="o",xlab="",ylab="",main=paste(gender[sex+1],uid[i]))
     yp=QLogist(a1=a[1]+a1.id[i],a2=a[2],a3=a[3],a4=a[4],x=xx)
     lines(xx,yp,lwd=3,col=22)
 }
}
