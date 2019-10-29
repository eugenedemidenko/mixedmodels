familyI <-
function()
{
    dump("familyI","c:\\MixedModels\\Chapter03\\familyI.r")
    library(nlme)
    famdat=read.table("c:\\MixedModels\\Chapter02\\Family.txt",header=T,stringsAsFactors=F)
    lmout=lme(Weight~Height,random=~1|FamilyID,data=famdat,method="ML")
    print(summary(lmout))
    fam.uniq=unique(famdat$FamilyID)
    nfam=length(fam.uniq)
    plot(famdat$Height,famdat$Weight,type="n",xlab="Height, inches",ylab="Weight, pounds",main=paste("Weight versus Height for",nfam,"families"))
    text(famdat$Height,famdat$W,famdat$FamilyID)
    af=lmout$coefficients$fixed
    ar=lmout$coefficients$random$FamilyID
    lines(famdat$Height,af[1]+af[2]*famdat$Height,col=2,lwd=5)
    maxH=max(famdat$Height)+.1
    x=range(famdat$Height)
    for(i in 1:nfam)
    {
       lines(x,af[1]+ar[i]+af[2]*x,col=3,lwd=1)
       text(x[2],af[1]+ar[i]+af[2]*x[2],i,adj=0)
    }
    legend(x[1],max(famdat$Weight),c("Fixed effect","Random effect"),lty=1,col=2:3,lwd=c(5,1))
 }
