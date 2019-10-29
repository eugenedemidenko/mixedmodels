untrlme <-
function()
{
    dump("untrlme","C:\\MixedModels\\Chapter10\\untrlme.r")
    library(nlme)
    dat = read.table("c:\\MixedModels\\Chapter10\\DEregrowth.dat",stringsAsFactor=F)
    names(dat)=c("TreatmentGroup","MiceID","TumorVolume","Day")
    utrgr=unique(dat[,1])
    plot(1,1,xlim=c(-10,30),ylim=c(-5,1),main=utrgr[1],type="n")
    xd=seq(from=-10,to=30,by=.1)
    daUNTR=dat[dat$TreatmentGroup==utrgr[1] | dat$Day<=0,]
    y=log(daUNTR$TumorVolume); day=daUNTR$Day
    id=daUNTR$MiceID; uid=unique(id); nid=length(uid)
    for(j in 1:nid)
    {
        yi=y[id==uid[j]];xi=day[id==uid[j]]
        lines(xi,yi);points(xi,yi,pch=16)
    }
    o=lme(fixed=y~day,random=~1|id)
    print(summary(o))
    aUNTR=as.vector(o$coefficients$fixed)
    lines(xd,aUNTR[1]+aUNTR[2]*xd,col=2,lwd=3)
    segments(0,-6,0,2,col=3) 
    covpar <- o$varFix
    DT <- log(2)/aUNTR[2]
    lines(c(-12,DT,DT),c(aUNTR[1]+log(2),aUNTR[1]+log(2),-6),col=3)
    SE.TD <- DT/aUNTR[2] * sqrt(covpar[2, 2])
    text(10,-4,paste("Doubling Time =", round(DT,1), "\nSE =",round(SE.TD,2)),adj=0,cex=1.25)    
    
    
}
