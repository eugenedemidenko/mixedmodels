KSimage <-
function(job=1)
{
dump("KSimage","c:\\MixedModels\\Chapter12\\KSimage.r")
if(job==1)
{
	par(mfrow = c(1, 2), mar = c(1, 1, 3, 1), omi = c(0, 0, 0, 0))
	d <- scan("c:\\MixedModels\\Chapter12\\grp11.pgm",what="")
	d <- as.numeric(d[9:length(d)])
	nr <- d[1];nc <- d[2]
	d <- matrix(d[4:length(d)], nrow = nr, ncol = nc)
	image(1:nr, 1:nc, d, xlab = "", ylab = "", axes = F,col=gray(0:255/255))
	mtext(side = 3, "Control", line = 0.25, cex = 2)
	d <- scan("c:\\MixedModels\\Chapter12\\grp51.pgm",what="")
	d <- as.numeric(d[9:length(d)])
	nr <- d[1];nc <- d[2]
	d <- matrix(d[4:length(d)], nrow = nr, ncol = nc)
	image(1:nr, 1:nc, d, xlab = "", ylab = "", axes = F,col=gray(0:255/255))
	mtext(side = 3, "Drug+Radiation", line = 0.25, cex = 2)
}
if(job==2)
{
	par(mfrow = c(1, 2), mar = c(4, 4, 3, 1))
	d <- scan("c:\\kluwer\\image\\sujatha\\grp11.pgm",what="")
	d <- as.numeric(d[9:length(d)])
	nr <- d[1]
	nc <- d[2]
	J1 <- nr * nc
	d1 <- d[4:length(d)]
	d <- scan("c:\\kluwer\\image\\sujatha\\grp51.pgm",what="")
	d <- as.numeric(d[9:length(d)])
	nr <- d[1]
	nc <- d[2]
	print(c(nr, nc))
	J2 <- nr * nc
	d2 <- d[4:length(d)]
	h1 <- f1 <- h2 <- f2 <- rep(0, 256)
	for(i in 0:255) {
		h1[i + 1] <- length(d1[d1 == i])/length(d1)
		h2[i + 1] <- length(d2[d2 == i])/length(d2)
	}
	for(i in 2:256) {
		f1[i] <- f1[i - 1] + h1[i]
		f2[i] <- f2[i - 1] + h2[i]
	}
	f1[256] <- f2[256] <- 1
	matplot(cbind(0:255, 0:255), cbind(h1, h2), type = "l", col = 1,xlab="",ylab="")
	mtext(side = 2, "Probability", line = 2.5, cex = 1.75)
	mtext(side = 3, "Histogram", line = 1, cex = 1.75)
	matplot(cbind(0:255, 0:255), cbind(f1, f2), type = "l", col = 1,xlab="",ylab="")
	lines(0:255, f1, lwd = 3)
	lines(0:255, f2, lwd = 3, lty = 2)
	mf <- max(abs(f1 - f2))
	jm <- (0:255)[abs(f1 - f2) == mf]
	segments(jm, -0.25, jm, 0.6)
	text(jm+5, 0.63, paste("max=", round(mf, 3)),adj=0)
	mtext(side = 3, "Distribution function", line = 1, cex = 1.75)
	legend(0, 1, c("Control", "Treatment"), lty = 1:2, cex = 1.25,lwd = 3)
	mtext(side = 1, "Grayscale level, byte", line = -1, outer = T,cex = 1.5)
	J <- (J1 * J2)/(J1 + J2)
	lambda <- mf * (sqrt(J) + 0.11/sqrt(J) + 0.12)
	j <- 1:10000
	js <- rep(1, 10000)
	js[seq(from = 2, to = 10000, by = 2)] <- -1
	Q <- 2 * sum(js * exp(-2 * j^2 * lambda^2))
	cat("\nbyte.max =",jm," max.cdf.distance =",round(mf,3)," lambda =",round(lambda,3)," QKS =", round(Q,4),"\n")
}
}
