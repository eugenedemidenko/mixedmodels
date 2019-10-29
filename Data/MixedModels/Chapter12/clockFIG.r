clockFIG <-
function()
{
	dump("clockFIG", "c:\\MixedModels\\Chapter12\\clockFIG.r")
	matR=function(M)	#matrix reflection about y-axis
	{
		nr=nrow(M);nc=ncol(M)
		MR=M
		for(i in 1:nc) MR[,nc-i+1]=M[,i]
		return(MR)
	}
	
	c1 <- scan("c:\\MixedModels\\Chapter12\\clock1.pgm",what="")
	nr1 <- as.numeric(c1[2]); nc1 <- as.numeric(c1[3])
	M1 <- matrix(as.numeric(c1[5:length(c1)]), nrow = nr1, ncol = nc1)
	c2 <- scan("c:\\MixedModels\\Chapter12\\clock2.pgm",what="")
	nr2 <- as.numeric(c2[2]); nc2 <- as.numeric(c2[3])
	M2 <- matrix(as.numeric(c2[5:length(c2)]), nrow = nr2, ncol = nc2)
	cL1 <- matrix(c(327, 148, 376, 256, 269, 344, 128, 280, 191, 101, 252, 111), nrow = 2)
	cL2 <- matrix(c(302, 154, 349, 258, 245, 343, 108, 283, 169, 108, 228, 116), nrow = 2)
	cL1 <- cL1[, 1:5];	cL2 <- cL2[, 1:5]
	nLc1 <- ncol(cL1)
	X <- matrix(0, ncol = 6, nrow = 2 * nLc1)
	for(i in 1:nLc1) {
		i1 <- 1 + 2 * (i - 1)
		i2 <- 2 * i
		X[i1, 1] <- 1
		X[i1, 2:3] <- cL1[1:2, i]
		X[i2, 4] <- 1
		X[i2, 5:6] <- cL1[1:2, i]
	}
	xx <- t(X) %*% X
	xy <- t(X) %*% as.vector(cL2)
	beta <- b <- solve(xx) %*% xy
	par(mfrow = c(2, 2), mar = c(1, 1, 3, 1))
	image(1:nr1, 1:nc1, matR(M1), xlab = "", ylab = "", axes = F,col=gray(0:255/255))
	mtext(side = 3, "Source image", line = 0.25, cex = 1.25)
	for(i in 1:nLc1)
		points(cL1[1 + 2 * (i - 1)], nc1 - cL1[2 * i] + 1,	pch = 16, cex = 1.25)
	image(1:nr2, 1:nc2, matR(M2), xlab = "", ylab = "", axes = F,col=gray(0:255/255))
	mtext(side = 3, "Target image", line = 0.25, cex = 1.25)
	for(i in 1:nLc1)
		points(cL2[1 + 2 * (i - 1)], nc1 - cL2[2 * i] + 1, pch = 16, cex = 1.25)
	del <- matrix(nrow = nr1, ncol = nc1)
	r1 <- min(nr1, nr2)
	r2 <- min(nc1, nc2)
	image(1:r1, 1:nc1, matR(M1[1:r1, 1:r2] - M2[1:r1, 1:r2]), xlab = "", ylab = "", axes = F,col=gray(0:255/255))
	mtext(side = 3, "Difference before alignment", line = 0.25,	cex = 1.25)
	for(i in 1:nr1)
		for(j in 1:nc1) {
			p1 <- round(b[1] + b[2] * i + b[3] * j)
			q1 <- round(b[4] + b[5] * i + b[6] * j)
			if(q1 > 0 & q1 <= nc2 & p1 > 0 & p1 <= nr2)	del[i, j] <- M1[i, j] - M2[p1, q1]
		}
	image(1:nr1, 1:nc1, matR(del), xlab = "",	ylab = "", axes = F,col=gray(0:255/255))
	mtext(side = 3, "Difference after alignment", line = 0.25, cex = 1.25)
}
