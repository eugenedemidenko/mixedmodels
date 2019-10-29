logG <-
function()
{
    dump("logG","c:\\MixedModels\\Chapter06\\logG.r")
    library(nlme);library(lattice)
    dat=read.csv("c:\\MixedModels\\Chapter06\\TUMspher.txt")    
    xyplot(lntumvol~day|id,type="b",data=dat)
    out.nlme<-nlme(lntumvol~a1-a2*exp(-a3*day),fixed=a1+a2+a3~1,random=a1~1|id, data=dat,start=c(6,5,0.07))
    summary(out.nlme)
}
