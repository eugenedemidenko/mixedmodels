trnlme <-
function()
{
    library(nlme)
    dat = read.table("c:\\MixedModels\\Chapter10\\DEregrowth.dat",stringsAsFactor=F)
    dump("trnlme","C:\\MixedModels\\Chapter10\\trnlme.r")
    names(dat)=c("TreatmentGroup","MiceID","TumorVolume","Day")
    xd=seq(from=1,to=30,by=.1)
    utrgr=unique(dat[,1])
    par(mfrow=c(1,3))
    for(ig in 2:4)
    {
        print(utrgr[ig])
        plot(1,1,xlim=c(0,30),ylim=c(-5,1),main=utrgr[ig],type="n",xlab="",ylab="")
        daUNTR=dat[dat$TreatmentGroup==utrgr[ig] & dat$Day>0,]
        y=log(daUNTR$TumorVolume); day=daUNTR$Day
        id=daUNTR$MiceID; uid=unique(id); nid=length(uid)
        for(j in 1:nid)
        {
            yi=y[id==uid[j]];xi=day[id==uid[j]]
            lines(xi,yi);points(xi,yi,pch=16)
        }
        o <- nlme(model=y~log(exp(a1+a2*day)+exp(a3-a4*day)),fixed=list(a1~1,a2~1,a3~1,a4~1), random=a1+a3~1|id,start=c(-7.,0.2,-0.8,0.2))
        print(summary(o))
        a <- as.vector(o$coefficients[[1]])
        yfit=log(exp(a[1]+a[2]*xd)+exp(a[3]-a[4]*xd))
        lines(xd,yfit,col=2,lwd=3)
        TR <- (log(a[4]/a[2])-a[1]+a[3])/(a[4]+a[2])
        ytr <-log(exp(a[1]+a[2]*TR)+exp(a[3]-a[4]*TR))
        segments(-1, 0, 30, 0, col=3,lty=2)
        segments(TR, ytr, TR, -7, col=3)
        d1 <- -1/(a[2] + a[4])
        d2 <- - d1
        al <- log(a[4]/a[2])
        d3 <- (-a[2]*al+a[2]*(a[1]-a[3])-a[2]-a[4])/(a[2]+a[4])^2/a[2]
        d4 <- (-a[4]*al+a[4]*(a[1]-a[3])+a[2]+a[4])/(a[2]+a[4])^2/a[4]
        dd <- c(d1, d3, d2, d4)
        SETR <- sqrt(t(dd) %*% o$varFix %*% dd)
        text(20,-4,paste("TR=",round(TR,1),"\nSE=",round(SETR,1),sep=""),adj=0,cex=2)
    }
    mtext(side=1,"Days Post Treatment",outer=T,cex=1.5,line=-1.5)
    mtext(side=2,"LOG Tumor Volume",outer=T,cex=1.5,line=-1.5)
    
}
