logtree <-
function()
{
    dump("logtree","c:\\MixedModels\\Chapter06\\logtree.r")
    library(nlme);library(lattice)
    dat=read.csv("c:\\MixedModels\\Chapter06\\trunktree.txt")    
    xyplot(trunk~day|id,data=dat,type="b")
    out.nlme<-nlme(trunk~a1/(1+exp(a2-a3*day)),fixed=a1+a2+a3~1,random=a1~1|id, data=dat,start=c(200,2,0.003))
    summary(out.nlme)
}
