ups2 <-
function(y,s2i,alpha=0.05,maxit=10,eps=0.0001)
{
    dump("ups2","c:\\MixedModels\\Chapter05\\ups2.r")
    m0=sum(y/s2i)/sum(1/s2i)
    cff0=sum((y-m0)^2/s2i)
    n=length(s2i)
    qc=qchisq(alpha,df=n-1)
    if(cff0<=qc) return(0)
    TL=0
    for(iter in 1:maxit)
    {
        mut=sum(y/(TL+s2i))/sum(1/(TL+s2i))
        dbd1=-1/sum(1/(TL+s2i))*sum(y/(TL+s2i)^2)
        dbd2=sum(y/(TL+s2i))*sum(1/(TL+s2i)^2)/sum(1/(TL+s2i))^2
        dbd=dbd1+dbd2
        p1=-2*dbd*sum((y-mut)/(TL+s2i))-sum((y-mut)^2/(TL+s2i)^2)
        PS2=sum((y-mut)^2/(TL+s2i))
        TLL=TL-(PS2-qc)/p1 # Newton's iteration
        if(abs(TLL-TL)<eps) break
        TL=TLL
    }
    return(TL)
}
