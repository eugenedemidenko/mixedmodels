schiz <-
function () 
{
    dump("schiz", "c:\\MixedModels\\Chapter12\\schiz.r")
    cc = "c:\\MixedModels\\Chapter12\\schiz\\case"
    for (i in 1:30) {
        par(mfrow = c(3, 7), mar = c(1, 1, 1, 1),omi=c(0,0,.25,0))
        for (j in 1:21) {
            d = scan(paste(cc, i, "\\case", i, ".0", j + 50,".pgm", sep = ""), what = "", quiet = T)
            dm = matrix(as.numeric(d[12:length(d)]), nrow = 256, ncol = 256)
            image(1:256, 1:256, dm, col = gray(0:255/255),xlab="",ylab="",axes=F)
			mtext(side=3,paste("Frame",j+50),line=.2,cex=.75)
        }
		mtext(side=3,paste("Case",i),outer=T,cex=1.25,line=.25)
    }
}
