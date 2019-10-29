heightlog2S <-
function(sex=0)
{
 dump("heightlog2S","c:\\MixedModels\\Chapter08\\heightlog2S.r")
 gender=c("Boy","Girl")
 if(sex==1) a0=c(165,-2.2,-0.42,0.039) else a0=c(182,-1,-0.16,0.018) #start
 da=height.dat[height.dat$sex==sex,] #extract gender set
 nkids=nrow(da)
 id=da$id
 xx=seq(from=7,to=18,by=.1) #values to plot the curve
 uid=unique(id);nud=length(uid)
 par4cov=matrix(nrow=nud,ncol=4+4^2) #matrix of all fits
 par(mfrow=c(7,10),mar=c(2,2,1,1))
 Scov=matrix(0,4,4)
 Sa=rep(0,4)
 SSres=0;dfpool=0
 for(i in 1:nud)
 {
  x=da$year[id==uid[i]];y=da$height[id==uid[i]]
  y=y[order(x)];x=x[order(x)] # order to make lines
  plot(x,y,type="o",xlab="",ylab="",main=paste(gender[sex+1],uid[i]))
  oi=try(nls(y~a1/(1+exp(a2-a3*x-a4*x^2)),start=list(a1=a0[1],a2=a0[2],a3=a0[3],a4=a0[4]),nls.control(maxiter=200))) #prevent from stop
  if(attr(oi,"class")!="try-error") 
  {
     a=coef(oi)
     par4cov[i,1:4]=a
     yp=a[1]/(1+exp(a[2]-a[3]*xx-a[4]*xx^2))
     lines(xx,yp,col=2,lwd=3)
     covpar=summary(oi)$cov.unscaled #cov matrix/sigma2i
     par4cov[i,5:20]=as.vector(covpar)
     SSres=SSres+sum((y-predict(oi))^2)
     dfpool=dfpool+length(y)-4
     i.covpar=solve(covpar)
     Sa=Sa+i.covpar%*%a
     Scov=Scov+i.covpar
  }  
 }
 s2pool=SSres/dfpool # pooled variance
 print(gender[sex+1])
 meanpar=as.data.frame(matrix(ncol=2,nrow=4))
 names(meanpar)=c("Mean","Weighted Mean")
 row.names(meanpar)=c("a1","a2","a3","a4")
 for(i in 1:4) meanpar[i,1]=mean(par4cov[,i],na.rm=T)
 meanpar[,2]=solve(Scov)%*%Sa #weighted mean
 print(meanpar)
 par4cov[,5:20]=par4cov[,5:20]*s2pool # matrices Ti
 # par4cov.height0=heightlog2S(sex=0)
 return(par4cov) # return all ind fits for further analysis
}
