wbf <-
function()
{
    dump("wbf","c:\\MixedModels\\Chapter09\\wbf.r") # save the code
    #Influence analysis for Chapter 9
    #Example 1: Women's body fat, data in WomenBF.dat
    dat=read.table("c:\\MixedModels\\Chapter09\\WomenBF.dat")
    n=nrow(dat)
    outlm=lm(Fat~Triceps+Thigh,data=dat)
    betaALL=coef(outlm)[2:3]
    ri=outlm$residuals
    betai=matrix(nrow=n,ncol=2)
    dbdw1=dbdw0=matrix(nrow=n,ncol=2)
    X=as.matrix(cbind(rep(1,n),dat[,2:3]))
    iXtX=solve(t(X)%*%X)
    dbdw1=iXtX%*%t(ri*X)
    for(i in 1:n)
    {
        w=rep(1,n); w[i]=0
        outi=lm(Fat~Triceps+Thigh,data=dat,weights=w)
        betai[i,]=(coef(outi)[2:3]-betaALL)/betaALL*100 
    }
    matplot(1:n,betai,lty=1,type="b",xlab="Case deleted",
            ylab="% beta change")
    lines(1:n,-100*dbdw1[2,],lty=2)
    lines(1:n,-100*dbdw1[3,],lty=2,col=2)
    legend(12,2,c("1=Triceps case deletion","2=Thigh case deletion",
                  "Triceps db/dw=1","Thigh db/dw=1"),
           lty=c(1,1,2,2),col=c(1,2,1,2))
}
