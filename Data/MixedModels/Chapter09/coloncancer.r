coloncancer <-
function()
{
    library(nlme)
    dump("coloncancer","c:\\MixedModels\\Chapter09\\coloncancer.r")
    dat=read.table("c:\\MixedModels\\Chapter09\\coloncancer.dat")
    out.lme=lme(LNtotexpplus1~female+black+age+stage2+stage3+xchrlson+t1,random=~1|id,data=dat)
    summary(out.lme)
    
}
