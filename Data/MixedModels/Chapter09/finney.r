finney <-
function()
{
    dump("finney","c:\\MixedModels\\Chapter09\\finney.r")
    dat=read.table("c:\\MixedModels\\Chapter09\\Finney.dat")
    n=nrow(dat)
    lnv=log(dat$Volume)
    lnr=log(dat$Rate)
    dat=cbind(dat,lnv,lnr)
    outL=glm(y~lnv+lnr,data=dat,family=binomial)
    beta1OUT=matrix(nrow=n,ncol=2)
    for(i in 1:n)
    {
        w=rep(1,n);w[i]=0
        outLi=glm(y~lnv+lnr,weights=w,data=dat,family=binomial)
        beta1OUT[i,]=coef(outLi)[2:3]
    }
    beta1OUT
}
