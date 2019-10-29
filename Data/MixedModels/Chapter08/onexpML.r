onexpML <-
function(beta=.1,sigma=.1,omega=.2,n=10,N=100,ss=2)
{
    dump("onexpML","c:\\MixedModels\\Chapter08\\onexpML.r")
    set.seed(ss)
    LInt=function(a,y,beta,sigma,omega)
    {
        Sy2=sum(y^2);Sy=sum(y)
        Sya=Sy2-2*Sy*exp(a)+length(y)*exp(2*a)
        t1=-Sya/2/sigma^2
        t2=-(a-beta)^2/2/omega^2/sigma^2
        f=exp(t1+t2)
        return(f)
    }
    
    
    ai=rnorm(N,mean=beta,sd=sigma*omega)
    eps=matrix(sigma*rnorm(N*n),nrow=N,ncol=n)
    y=exp(ai)%*%t(rep(1,n))+eps
    bs=L=seq(from=beta-1,to=beta+1,length=100)
    for(ib in 1:100)
    {
        L[ib]=0
        for(i in 1:N)
        {
            yi=y[i,]
            Li=log(integrate(LInt,y=yi,sigma=sigma,beta=bs[ib],omega=omega,lower=-Inf,upper=Inf)$value)
            L[ib]=L[ib]+Li
        }
    }
    plot(bs,L,type="o")
    bs[L==max(L)]
}
