hypoxiaRAT <-
function()
{
dump("hypoxiaRAT", "c:\\MixedModels\\Chapter12\\hypoxiaRAT.r")
n <- 1024
x <- 1:n
# we will save the graph in the file 
bmp(file="c:\\MixedModels\\Chapter12\\Hypoxia\\hypoxiaRAT.bmp",width=500,height=1500)
par(mfrow=c(8,2),mar=c(1,1,1,1),omi=c(0,0.25,0.25,0))
for(i in 1:8)
for(igr in 1:2)
{
   if(igr==1) cc="c" else cc=""
   fn=paste("c:\\MixedModels\\Chapter12\\Hypoxia\\Group",igr,"\\_",i,cc,"_1a_p.pgm",sep="")
   d <- scan(fn,what="")
   d <- matrix(as.numeric(d[12:length(d)]), n, n)
   image(x, x, d, xlab = "", ylab = "", axes = F,col=gray(0:255/255))
   if(igr==1) mtext(side=2,paste("Rat #",i,sep=""),line=0.25,cex=1.25)
   if(i==1 & igr==1) mtext(side=3,"Control group",line=1,cex=1.5)
   if(i==1 & igr==2) mtext(side=3,"Treatment group",line=1,cex=1.5)
}
dev.off() # saving the graph
}
