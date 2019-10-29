nlsnif <-
function()
{
    dump("nlsnif","c:\\MixedModels\\Chapter09\\nlsnif.r")
    dat=read.table("c:\\MixedModels\\Chapter09\\NLSNIF.dat")
    n=nrow(dat)
    thetaDEL=matrix(ncol=4,nrow=n)
    outnls=nls(y~theta1+theta2/(1+exp(theta4*(x-theta3))),
               start=c(theta1=2000,theta2=3200,theta3=-8.3,theta4=1.3)
               ,data=dat)
    thetaALL=coef(outnls)
    for(i in 1:n)
    {
        w=rep(1,n);w[i]=0
        outnlsi=nls(y~theta1+theta2/(1+exp(theta4*(x-theta3))),
                    weights=w,start=list(theta1=thetaALL[1],
                                         theta2=thetaALL[2],theta3=thetaALL[3],
                                         theta4=thetaALL[4]),data=dat,
                    control=list(maxiter=500))
        thetaDEL[i,]=coef(outnlsi)
    }
    return(thetaDEL)
}
