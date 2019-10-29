LNGHint <-
function(K=13)
{
dump("LNGHint","c:\\MixedModels\\Chapter07\\LNGHint.r")

LNint=function(x,s,sigma)
{
    ex1=exp(-(x-s)^2/2/sigma^2)
    ex2=1/sigma/sqrt(2*pi)*exp(x)/(1+exp(x))
    return(ex1*ex2)
}
xw=gauher(K)
out=as.data.frame(matrix(ncol=4,nrow=6))
names(out)=c("s","sigma","GH quadrature","integrate")
j=0
for(s in c(-1,0,1))
    for(sigma in c(1,2))
    {
        j=j+1     
        ex=exp(s+sqrt(2)*sigma*xw[,1])
        out[j,1]=s;out[j,2]=sigma
        out[j,3]=sum(ex/(1+ex)*xw[,2])/sqrt(pi)
        out[j,4]=integrate(f=LNint,s=s,sigma=sigma,low=s-6*sigma,up=s+6*sigma)$value
    }
print(paste("K =",K))
return(out)
}
