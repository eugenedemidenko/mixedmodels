logMLE1 <-
function(dat, b, s2, maxiter = 20, eps = 1e-010, silent = 1)
{
	hi <- function(u, s2, B, ki)
	{
		n <- length(u)
		r <- rep(0, n)
		for(i in 1:n)
			r[i] <- ki * u[i] - u[i]^2/(2 * s2) - sum(log(1 + B * exp(u[i])))
		return(exp(r))
	}
	u2hi <- function(u, s2, B, ki)
	{
		n <- length(u)
		r <- rep(0, n)
		for(i in 1:n)
			r[i] <- ki * u[i] - u[i]^2/(2 * s2) - sum(log(1 + B * exp(u[i])))
		return(exp(r) * u^2)
	}
	ehi <- function(u, s2, B, j, ki)
	{
		n <- length(u)
		r <- rep(0, n)
		for(i in 1:n) {
			yy <- B * exp(u[i])
			p <- exp(ki * u[i] - u[i]^2/(2 * s2) - sum(log(1 + yy)))
			r[i] <- (p * yy[j])/(1 + yy[j])
		}
		return(r)
	}
	y <- dat[, 2]
	X <- cbind(rep(1, nrow(dat)), dat[, 3:ncol(dat)])
	r <- t(X) %*% y
	m <- length(b)
	di <- dat[, 1]
	N <- length(unique(di))
	bs2 <- c(b, s2)
	PI <- atan(1) * 4
	if(silent == 0)
		print("ML EFS: iter   grad    lambda     loglik      beta     s2")
	for(iter in 1:maxiter) {
		grad <- rep(0, m + 1)
		Hess <- matrix(0, m + 1, m + 1)
		loglik <-  - N/2 * log(s2) + sum(b * r) - N/2 * log(2 * PI)
		A <- -3.5 * sqrt(s2)
		B <- 3.5 * sqrt(s2)
		for(id in unique(di)) {
			yi <- dat[di == id, 2]
			sy <- sum(yi)
			ni <- length(yi)
			Xi <- cbind(rep(1, ni), dat[di == id, 3:ncol(dat)])
			Bi <- exp(Xi %*% b)
			Ihi <- integrate(f = hi, lower = A, upper = B, s2 = s2, B = Bi, ki = sy)$value
			loglik <- loglik + log(Ihi)
			I4 <- I3 <- rep(0, ni)
			for(j in 1:ni)
				I3[j] <- integrate(f = ehi, lower = A, upper = B, s2 = s2, B = Bi, j = j, ki = sy)$value
			I3 <- t(Xi) %*% I3
			is21 <- integrate(f = u2hi, lower = A, upper = 0, s2 = s2, B = Bi, ki = sy)$value
			is22 <- integrate(f = u2hi, lower = 0, upper = B, s2 = s2, B = Bi, ki = sy)$value
			is2 <- is21 + is22
			d1 <- t(Xi) %*% yi - I3/Ihi
			d2 <- -0.5/s2 + (0.5/s2^2 * is2)/Ihi
			grad[1:m] <- grad[1:m] + d1
			grad[m + 1] <- grad[m + 1] + d2
			Hess[1:m, 1:m] <- Hess[1:m, 1:m] + d1 %*% t(d1)
			Hess[m + 1, m + 1] <- Hess[m + 1, m + 1] + d2^2
			Hess[1:m, m + 1] <- Hess[1:m, m + 1] + d2 * d1
		}
		if((iter - 1 + silent) == 0)
			print(c(0, loglik, bs2))
		Hess[m + 1, 1:m] <- Hess[1:m, m + 1]
		iH <- solve(Hess)
		delta <- iH %*% grad
		for(lambda in 0:10) {
			lambdacur <- lambda
			bs2new <- bs2 + delta/2^lambda
			bnew <- bs2new[1:m]
			s2new <- bs2new[m + 1]
			A <- -3.5 * sqrt(s2new)
			B <- 3.5 * sqrt(s2new)
			if(s2new > 0) {
				logliknew <-  - N/2 * log(s2new) + sum(bnew * r) - N/2 * log(2 * PI)
				for(id in unique(di)) {
					yi <- dat[di == id, 2]
					sy <- sum(yi)
					ni <- length(yi)
					Xi <- cbind(rep(1, ni), dat[di == id, 3:ncol(dat)])
					Bi <- exp(Xi %*% bnew)
					Ihi <- integrate(f = hi, lower = A, upper = B, s2 = s2new, B = Bi, ki = sy)$value
					logliknew <- logliknew + log(Ihi)
				}
				if(logliknew > loglik)
					break
			}
		}
		loglik <- logliknew
		if(sum(abs(bs2 - bs2new)) < eps)
			return(list(b, s2, iH, iter))
		bs2 <- bs2new
		b <- bs2new[1:m]
		s2 <- bs2new[m + 1]
		if(silent == 0)
			print(c(iter, 1/2^lambdacur, sum(abs(grad)), logliknew, bs2new))
	}
	return(list(b, s2, iH, maxiter))
}
