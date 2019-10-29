logFS1 <-
function(dat, nSim = 100, bsigma, maxiter = 20, eps = 0.001, silent = 1)
{
	#Fixed Simulation approach with lambda=1 for logistic regression with random intercept
	u <- rnorm(nSim)
	es <- rep(1, nSim)
	y <- dat[, 2]
	X <- matrix(dat[, 3:ncol(dat)], ncol = (ncol(dat) - 2))
	di <- dat[, 1]
	N <- length(unique(di))
	m <- ncol(X) + 1
	m1 <- m + 1
	theta <- bsigma
	b <- theta[1:m]
	sigma <- theta[m1]
	if(silent == 0) print("Fixed Simulation: iter   lambda   loglik   grad   beta    sigma^2")
	for(iter in 1:maxiter) {
		der <- rep(0, m1)
		H <- matrix(0, m1, m1)
		loglik <-  - N * log(nSim)
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
			ei <- as.vector(exp(bri + sigma * ki * u - log(1 + esl) %*% ee))
			db <- es %*% t(ri) - (esl/(1 + esl)) %*% Xi
			dbs <- ki * u - u * (esl/(1 + esl)) %*% ee
			dbsc <- cbind(db, dbs)
			gi <- sum(ei)
			deri <- c((t(db) %*% ei), sum(dbs * ei))/gi
			der <- der + deri
			loglik <- loglik + log(gi)
			H <- H + (deri %*% t(deri))
		}
		iH <- solve(H)
		delta <- iH %*% der
		if(sum(abs(delta)) < eps)
		    return(list(b, sigma^2, iH, iter))		
		theta <- theta+delta
		b <- theta[1:m]
		sigma <- theta[m1]
		if(silent == 0)
			print(c(iter, lambda, loglik, sqrt(sum(der^2)), theta[1:m], theta[m1]^2))
	}
	return(list(b, sigma^2, iH, maxiter))
}
