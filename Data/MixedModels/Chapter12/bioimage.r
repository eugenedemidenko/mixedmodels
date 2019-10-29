bioimage <-
function()
{
dump("bioimage", "c:\\MixedModels\\Chapter12\\bioimage.r")
im <- c(1, 2, 3, 5)
nm <- c("Control", "Radiation", "Drug", "Radiation+Drug")
nr <- 512
nc <- 384
cc <- "c:\\MixedModels\\Chapter12\\grp"
jpeg("c:\\MixedModels\\Chapter12\\bioimage.jpg",height=1500,width=800)
par(mfrow = c(7, 4), mar = c(1, 1, 1, 1), omi = c(0, 0, .2, 0))
for(j in 1:7)
for(i in 1:4)
{
	d <- scan(paste(cc, as.character(im[i]), as.character(j), ".pgm", sep = ""),what="")
	d <- matrix(as.numeric(d[12:length(d)]), ncol = nc, nrow = nr)
	image(1:nr, 1:nc, d, axes = F, xlab = "", ylab = "",col=gray(0:255/255))
    if(j==1) mtext(side = 3, nm[i], line = 0.1,cex=1.5)	
}
dev.off()
}
