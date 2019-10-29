histgr <-
function()
{
dump("histgr","c:\\MixedModels\\Chapter12\\histgr.r")
par(mfrow=c(1,2),mar=c(1,1,1,1))
d=scan("c:\\MixedModels\\Chapter12\\IMAGE038.pgm",what="",quiet=T)
nc=as.numeric(d[2]);nr=as.numeric(d[3])
d=dr=matrix(as.numeric(d[5:length(d)]),ncol=nr,nrow=nc)
for(i in 1:nr) d[,i]=dr[,nr-i+1]
image(1:nc,1:nr,d,col=gray(0:255/255),main="Original image",xlab="",ylab="",axes=F)
d=scan("c:\\MixedModels\\Chapter12\\IMAGE038.eq.pgm",what="",quiet=T)
nc=as.numeric(d[2]);nr=as.numeric(d[3])
d=dr=matrix(as.numeric(d[5:length(d)]),ncol=nr,nrow=nc)
for(i in 1:nr) d[,i]=dr[,nr-i+1]
image(1:nc,1:nr,d,col=gray(0:255/255),main="After histogram equalization",xlab="",ylab="",axes=F)
}
