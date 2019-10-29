phototum <-
function()
{
    library(lattice);library(nlme)
    dump("phototum","c:\\MixedModels\\Chapter06\\phototum.r")
    datpt=read.csv("c:\\MixedModels\\Chapter06\\phototumdat.csv")
    print(datpt)
    xyplot(lntv~day|id,type="b",data=datpt)

    outnlme<-nlme(lntv~alpha+gamma*day+beta*(exp(-delta*day)-1), fixed=alpha+gamma+beta+delta~1,
                  random=pdDiag(alpha+beta~1),groups=~id,data=datpt,start=c(5,0.5,2.7,0.23))
    summary(outnlme)
    
}
