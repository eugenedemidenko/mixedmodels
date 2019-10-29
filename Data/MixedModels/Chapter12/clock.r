clock <-
function(job = 1, st = 2)
{
	dump("clock", "c:\\MixedModels\\Chapter12\\clock.r")
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
	print(as.vector(b))
	if(job == 1) {
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
		nn <- length(!is.na(as.vector(del)))
		rr <- del[!is.na(del)]
		print(mean(rr)/sqrt(var(rr)))
		return(sum(rr^2)/nn)
	}
	if(job == 2) {
		print("DUD algorithm")
		print(date())
		P1 <- nr1
		Q1 <- nc1
		P2 <- nr2
		Q2 <- nc2
		print(paste("P1=", P1, " Q1=", Q1, "P2=", P2, "Q2=", Q2))
		p <- 1:P1
		q <- 1:Q1
		set.seed(st)
		par(mfrow = c(2, 2), mar = c(0, 0, 2, 0))
		TT <- 7
		print(paste("TT =", TT))
		ee <- rep(1, TT)
		A <- cbind(beta, beta, beta, beta, beta, beta, beta)
		n <- P1 * Q1
		y <- as.vector(M1)
		h <- 5*c(2, 0.01, 0.01, 2, 0.01, 0.01)
		for(i in 1:6)
			A[i, 1 + i] <- A[i, 1 + i] + h[i]
		FF <- matrix(nrow = n, ncol = TT)
		SS0 <- 1e+100
		eps <- 0.001
		p <- 1:P1
		q <- 1:Q1
		for(j in 1:TT) {
			M2M <- M1
			for(ii in p)
				for(ji in q) {
					pp <- round(A[1, j] + A[2, j] * ii + A[3, 1] * ji)
					qq <- round(A[4, j] + A[5, j] * ji +A[6, j] * ii)
					if(pp <= P2 & pp > 0 & qq <= Q2 & qq >	0)	M2M[ii, ji] <- M2[pp, qq]
				}
			FF[, j] <- as.vector(M2M)
		}
		RES <- y %*% t(ee) - FF
		ssj <- diag(t(RES) %*% RES)
		os <- order(ssj)
		A <- A[, os]
		FF <- FF[, os]
		ssj <- ssj[os]
		SS <- ssj[1]
		for(iter in 0:10) {
			print(date())
			print(c(iter, SS, A[, 1]))
			FFd <- FF[, 2:TT] - FF[, 1]
			Ad <- A[, 2:TT] - A[, 1]
			FtA <- FFd %*% t(Ad)
			print(t(FtA) %*% FtA)
			iq <- solve(t(FtA) %*% FtA)
			print(diad(iq))
			res <- y - FF[, 1]
			dd <- (Ad %*% t(Ad)) %*% iq %*% t(FtA) %*% res
			a <- A[, 1] + dd
			M2M <- M1
			for(ii in p)
				for(ji in q) {
					pp <- round(a[1] + a[2] * ii + a[3] * ji)
					qq <- round(a[4] + a[5] * ji + a[6] * ii)
					if(pp <= P2 & pp > 0 & qq <= Q2 & qq >0)	M2M[ii, ji] <- M2[pp, qq]
				}
			FF1 <- as.vector(M2M)
			res <- y - FF1
			SS1 <- sqrt(sum(res^2/length(res)))
			#			if(SS1 >= SS)				break
			image(p, seq(from = Q1, to = 1, by = -1), M2M - M1, axes = F, xlab = "", ylab = "", zlim = c(-200, 200))
			mtext(side = 3, paste("Iter =", iter + 1, ", SS =",
				round(SS1, 1)), cex = 1.25)
			FF[, 2:TT] <- FF[, 1:(TT - 1)]
			FF[, 1] <- FF1
			A[, 2:TT] <- A[, 1:(TT - 1)]
			A[, 1] <- a
			ssj[2:TT] <- ssj[1:(TT - 1)]
			ssj[1] <- SS1
			SS <- SS1
		}
	}
	if(job == 3) {
		print(paste("GN algorithm, date:", date()))
		n <- nr1 * nc1
		npar <- 6
		h <- 5*c(1, 10/nr1, 10/nc1, 1, 10/nr1, 10/nc1)
		FF <- matrix(nrow = n, ncol = 6)
		for(iter in 0:10) {
			del <- matrix(ncol = nc1, nrow = nr1)
			for(i in 1:nr1)
				for(j in nc1) {
					p1 <- b[1] + b[2] * i + b[3] * j
					q1 <- b[4] + b[5] * i + b[6] * j
					if(q1 > 0 & q1 <= nc2 & p1 > 0 & p1 <= nr2) 
					{
						F0 <- M2
						del[i, j] <- M1[i, j] - M2[p1, q1]
					}
				}
			nn <- length(!is.na(as.vector(del)))
			RMSE <- sqrt(sum(del[!is.na(del)]^2)/nn)
			print(paste("Iter =", iter, " RMSE =", RMSE))
			for(de in 1:npar) {
				bb <- b
				bb[de] <- bb[de] + h[de]
				del <- matrix(ncol = nc1, nrow = nr1)
				for(i in 1:nr1)
					for(j in nc1) {
						p1 <- bb[1] + bb[2] * i + bb[
							3] * j
						q1 <- bb[4] + bb[5] * i + bb[
							6] * j
						if(q1 > 0 & q1 <= nc2 & p1 >
							0 & p1 <= nr2)
							del[i, j] <- M2[p1,
								q1]
					}
				FF[, de] <- as.vector(del - F0)/h[de]
			}
			FtA <- FFd %*% t(Ad)
			iq <- solve(t(FtA) %*% FtA)
			res <- y - FF[, 1]
			dd <- (Ad %*% t(Ad)) %*% iq %*% t(FtA) %*% res
			a <- A[, 1] + dd
			M2M <- M1
			for(ii in p)
				for(ji in q) {
					pp <- round(a[1] + a[2] * ii + a[3] *
						ji)
					qq <- round(a[4] + a[5] * ji + a[6] *
						ii)
					if(pp <= P2 & pp > 0 & qq <= Q2 & qq >
						0)
						M2M[ii, ji] <- M2[pp, qq]
				}
			FF1 <- as.vector(M2M)
			res <- y - FF1
			SS1 <- sum(res^2)
			#			if(SS1 >= SS)				break
			image(p, seq(from = Q1, to = 1, by = -1), M2M - M1,
				axes = F, xlab = "", ylab = "", zlim = c(-200
				, 200))
			mtext(side = 3, paste("Iter =", iter + 1, ", SS =",
				round(SS1)), cex = 1.25)
			FF[, 2:TT] <- FF[, 1:(TT - 1)]
			FF[, 1] <- FF1
			A[, 2:TT] <- A[, 1:(TT - 1)]
			A[, 1] <- a
			ssj[2:TT] <- ssj[1:(TT - 1)]
			ssj[1] <- SS1
			SS <- SS1
		}
	}
	return()
}
