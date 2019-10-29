poissMLE <-
function(dat, s2, maxiter = 20, eps = 1e-010, silent = 1)
{
	hi <- function(u, s2, B, ki)
	{
		n <- length(u)
		r <- rep(0, n)
		for(i in 1:n)
			r[i] <- ki * u[i] - (0.5 * u[i]^2)/s2 - exp(u[i]) * B
		return(exp(r))
	}
	u2hi <- function(u, s2, B, ki)
	{
		n <- length(u)
		r <- rep(0, n)
		for(i in 1:n)
			r[i] <- ki * u[i] - (0.5 * u[i]^2)/s2 - exp(u[i]) * B
		return(exp(r) * u^2)
	}
	y <- dat[, 2]
	X <- dat[, 3:ncol(dat)]
	b <- glm(y ~ X, family = poisson())$coefficients
	X <- cbind(rep(1, nrow(dat)), X)
	q <- t(X) %*% y
	m <- length(b)
	di <- dat[, 1]
	N <- length(unique(di))
	dd <- rep(0, m + 1)
	for(iter in 1:maxiter) {
		grad <- rep(0, m + 1)
		Hess <- matrix(0, m + 1, m + 1)
		loglik <-  - N/2 * log(s2) + sum(b * q)
		s2s <- 0
		A <- -10 * sqrt(s2)
		B <- 10 * sqrt(s2)
		for(id in unique(di)) {
			yi <- dat[di == id, 2]
			sy <- sum(yi)
			ni <- length(yi)
			Xi <- cbind(rep(1, ni), dat[di == id, 3:ncol(dat)])
			exi <- exp(Xi %*% b)
			Bi <- sum(exi)
			Ihi <- integrate(f = hi, lower = A, upper = B, s2 = s2, B = Bi, ki = sy)$value
			Ji <- integrate(f = hi, lower = A, upper = B, s2 = s2, B = Bi, ki = sy + 1)$value
			is2 <- integrate(f = u2hi, lower = A, upper = 0, s2 = s2, B = Bi, ki = sy)$value + integrate(
				f = u2hi, lower = 0, upper = B, s2 = s2, B = Bi, ki = sy)$value
			loglik <- loglik + log(Ihi)
			dd[1:m] <- t(Xi) %*% yi - Ji/Ihi * (t(Xi) %*% exi)
			dd[m + 1] <- -0.5/s2 + (0.5/s2^2 * is2)/Ihi
			grad <- grad + dd
			Hess <- Hess + dd %*% t(dd)
		}
		iH <- solve(Hess)
		delta <- iH %*% grad
		lambda <- 1
		for(ii in 1:10) {
			bnew <- b + lambda * delta[1:m]
			s2new <- s2 + lambda * delta[m + 1]
			if(s2new > 0) {
				logliknew <-  - N/2 * log(s2new) + sum(bnew * q)
				A <- -10 * sqrt(s2)
				B <- 10 * sqrt(s2)
				for(id in unique(di)) {
					yi <- dat[di == id, 2]
					sy <- sum(yi)
					ni <- length(yi)
					Xi <- cbind(rep(1, ni), dat[di == id, 3:ncol(dat)])
					exi <- exp(Xi %*% bnew)
					Bi <- sum(exi)
					Ihi <- integrate(f = hi, lower = A, upper = B, s2 = s2new, B = Bi, ki = sy)$value
					logliknew <- logliknew + log(Ihi)
				}
				if(logliknew > loglik)
					break
			}
			lambda <- lambda/2
		}
		if(sum(abs(bnew - b)/(1 + abs(b))) < eps)
			return(list(b, s2, iH, iter))
		b <- bnew
		s2 <- s2new
		if(silent != 1)
			print(c(iter, sum(abs(grad)), lambda, loglik, b, s2))
	}
	return(list(b, s2, iH, maxiter))
}
