potato <-
function()
{
dump("potato","c:\\MixedModels\\Chapter11\\potato.r")
par(mfrow=c(2,3),mar=c(1,1,3,1))
for(i in 1:6)
{
	potdat=scan(paste("c:\\MixedModels\\Chapter11\\pot",i,".pgm",sep=""),what="",quiet=T) # reading the image files
	nr=as.numeric(potdat[3]);nc=as.numeric(potdat[2]) 
	potdat=matrix(as.numeric(potdat[5:length(potdat)]),byrow=T,ncol=nc,nrow=nr) # the matrix image
	image(1:nr,1:nc,potdat,xlab="",ylab="",axes=F,col=grey(0:255/255))
	mtext(side=3, paste("Potato #",i,sep=""),cex=2,line=.5) 
}
}
