alckid <-
function()
{
    dump("alckid","c:\\MixedModels\\Chapter09\\alckid.r")
    dat=read.csv("c:\\MixedModels\\Chapter09\\alckid.dat")
    outL=glm(yn~agen+sensn+frn+alcmov,family=binomial,data=dat)
    summary(outL)
}
