logFSL1 <-
function(dat, xwFSL=gauher(13), bsigma, maxiter = 20, eps = 0.001, silent = 1)
{
	#Fixed GH likelihood method with lambda=1 for logistic regression with random intercept
	nSim <- nrow(xwFSL)
	u <- xwFSL[, 1] * sqrt(2)
	wFSL <- xwFSL[, 2]/sqrt(atan(1) * 4)
	y <- dat[, 2]
	es <- rep(1, nSim)
	X <- matrix(dat[, 3:ncol(dat)], ncol = (ncol(dat) - 2))
	di <- dat[, 1]
	N <- length(unique(di))
	m <- ncol(X) + 1
	m1 <- m + 1
	theta <- bsigma
	b <- theta[1:m]
	sigma <- theta[m1]
	if(silent == 0) print("FSL, Fixed Sample Likelihood: iter   loglik   grad   beta    sigma^2")
	loglik.old <-  - 10^10
	for(iter in 1:maxiter) {
		der <- rep(0, m1)
		H1 <- H3 <- H2 <- matrix(0, m1, m1)
		loglik <- 0
		for(id in unique(di)) {
			yi <- y[di == id]
			ni <- length(yi)
			ee <- rep(1, ni)
			Xi <- cbind(ee, X[di == id,  ])
			ri <- t(Xi) %*% yi
			ki <- sum(yi)
			Xib <- Xi %*% b
			bri <- sum(ri * b)
			sl <- es %*% t(Xib) + (sigma * u) %*% t(ee)
			esl <- exp(sl)
			esl2 <-  - esl/(1 + esl)^2
			ei <- (as.vector(exp(bri + sigma * ki * u - log(1 + esl) %*% ee))) * wFSL
			db <- es %*% t(ri) - (esl/(1 + esl)) %*% Xi
			dbs <- ki * u - u * (esl/(1 + esl)) %*% ee
			dbsc <- cbind(db, dbs)
			gi <- sum(ei)
			qq <- cbind(db, dbs)
			Ei <- diag(ei, nSim, nSim)
			H2 <- H2 + (t(qq) %*% Ei %*% qq)/gi
			H3 <- H3 + (t(qq) %*% Ei %*% qq)/gi
			deri <- c((t(db) %*% ei), sum(dbs * ei))/gi
			der <- der + deri
			loglik <- loglik + log(gi)
			H1 <- H1 + (deri %*% t(deri))
		}
		iH <- solve(H1)
		delta <- iH %*% der
        if(max(abs(delta))<eps) break
		theta <- theta + delta
		b <- theta[1:m]
		sigma <- theta[m1]
		if(silent == 0)
			print(c(iter, loglik, sqrt(sum(der^2)), theta[1:m], theta[m1]^2))
		if(sum(abs(delta)) < eps)
			return(list(b, sigma^2, iH, iter))
	}
	return(list(b, sigma^2, iH, maxiter))
}
