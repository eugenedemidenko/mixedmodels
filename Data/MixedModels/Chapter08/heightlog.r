heightlog <-
function()
{
 dump("heightlog","c:\\MixedModels\\Chapter08\\heightlog.r")
 gender=c("Boy","Girl")
 par(mfrow=c(1,2))
 ag=seq(from=7,to=18,by=.1)
 for(sex in 0:1)
 {
  print(gender[sex+1])
  da=height.dat[height.dat$sex==sex,]
  matplot(da$year,da$height,type="n",xlab="",ylab="",ylim=c(100,200))
  title(gender[sex+1])
  uid=unique(da$id);nud=length(uid)
  for(i in 1:nud)
  {
    x=da$year[da$id==uid[i]];y=da$height[da$id==uid[i]]
    y=y[order(x)];x=x[order(x)]
    lines(x,y);points(x,y)
  }
  mtext(side=2,"Height, cm",line=2.5,cex=1.5)
  o=nls(height~a1/(1+exp(a2-a3*year)),data=da,start=c(a1=180,a2=1,a3=.2))
  print(summary(o))
  a=coef(o)
  lines(ag,a[1]/(1+exp(a[2]-a[3]*ag)),col=2,lwd=3)
  
  o=nls(height~a1/(1+exp(a2-a3*year-a4*year^2)),data=da,start=list(a1=a[1],a2=a[2],a3=a[3],a4=0),nls.control(maxiter=200))
  print(summary(o))
  a=coef(o)
  lines(ag,a[1]/(1+exp(a[2]-a[3]*ag-a[4]*ag^2)),col=3,lwd=3)
  acc=(sqrt(2*a[4])-a[3])/(2*a[4])
  segments(acc,0,acc,300,lty=2,col=3)
  text(14,100,paste("Acceleration at",round(acc,1),"years old"))
  legend(7,200,c("Logistic","QLogistic"),lwd=3,lty=1,col=2:3,cex=1.25)
 }
 mtext(side=1,"Age, years",cex=1.5,line=-1.5,outer=T)
}
