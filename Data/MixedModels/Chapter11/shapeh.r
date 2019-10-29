shapeh <-
function(sigma=.01,N=15)
{
    dump("shapeh","c:\\MixedModels\\Chapter11\\shapeh.r")
    angtrue=(2*pi)*(0:5)/6
    mutrue=cbind(cos(angtrue),sin(angtrue))
    par(mfrow=c(4,4),mar=c(0,0,0,0))
    plot(c(-1.25,1.25),c(-1.25,1.25),type="n",axes=F,xlab="",ylab="")
    lines(c(mutrue[,1],mutrue[1,1]),c(mutrue[,2],mutrue[1,2]),lwd=4)
    text(0,0,"True shape",cex=1.5)
    u6=rep(1,6);R=matrix(ncol=2,nrow=2)    
    for(i in 1:15)
    {
        mur=mutrue+matrix(sigma*rnorm(2*6),ncol=2,nrow=6)
        ti=rnorm(2,mean=0,sd=.1)
        ri=runif(1,min=.9,max=1.1)
        theta=runif(1,min=0,max=2*pi)
        R[1,1]=cos(theta);R[1,2]=-sin(theta)
        R[2,1]=sin(theta);R[2,2]=cos(theta)
        ui=u6%*%t(ti)+mur%*%t(solve(R))/ri
        plot(c(-1.25,1.25),c(-1.25,1.25),type="n",axes=F,xlab="",ylab="")
        polygon(c(ui[,1],ui[1,1]),c(ui[,2],ui[1,2]),col=3)
        lines(c(ui[,1],ui[1,1]),c(ui[,2],ui[1,2]),lwd=4)
    }
}
