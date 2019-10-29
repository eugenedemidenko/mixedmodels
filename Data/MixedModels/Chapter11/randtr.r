randtr <-
function(roeps=.2,distr="run",m=3)
{
    par(mfrow=c(4,4),mar=c(1,1,1,1))
    ang=(2*pi)/m*(0:(m-1)) # angles for vertices
    Prot=matrix(ncol=2,nrow=2)
    for(itr in 1:16)
    {
        plot(c(-1.5,1.5),c(-1.5,1.5),type="n",xlab="",ylab="",axes=F)
        if(distr=="run") leng=runif(n=m,min=1-roeps,max=1+roeps)
        else leng=rnorm(n=m,mean=1,sd=roeps) #rotation distribution
        x=leng*cos(ang);y=leng*sin(ang)
        theta=runif(n=1,min=0,max=2*pi) #rotation angle
        Prot[1,1]=cos(theta);Prot[1,2]=-sin(theta) #rotation matrix
        Prot[2,1]=sin(theta);Prot[2,2]=cos(theta)
        xtr=Prot%*%t(cbind(x,y))
        xtr=cbind(xtr,xtr[,1]) # connectivity
        polygon(xtr[1,],xtr[2,],col=2)
        lines(xtr[1,],xtr[2,],lwd=4)
    }
}
