logVARLINK1 <-
function(dat, b, s2, maxiter = 20, eps = 0.001, silent = 1)
{
	li1 <- function(Bj, k, u = 0, eps = 1e-008, MaxIter = 100)
	{
		for(iter in 1:MaxIter) {
			it <- iter
			eu <- exp(u)
			d1 <-  - k + eu * sum(Bj/(1 + Bj * eu))
			d2 <- eu * sum(Bj/(1 + Bj * eu)^2)
			unew <- u - d1/d2
			if(abs(unew)>6) return(c(6,it))

			if(abs(unew - u) < eps)
				break
			u <- unew
		}
		return(c(u, it))
	}
	twoprob <- function(s, sig2)
	{
		s1 <- s/sqrt(5.2448 + sig2)
		s2 <- s/sqrt(1.6944 + sig2)
		p1 <- pnorm(s1)
		p2 <- pnorm(s2)
		II <- 0.4353 * p1 + 0.5647 * p2
		derII <- (0.4353 * dnorm(s1))/sqrt(5.2448 + sig2) + (0.5647 * dnorm(s2))/sqrt(1.6944 + sig2)
		return(cbind(as.vector(II), as.vector(derII)))
	}
	di <- dat[, 1]
	m <- ncol(dat) - 1
	N <- nrow(dat)
	X <- cbind(rep(1, N), dat[, 3:ncol(dat)])
	y <- dat[, 2]
	eem <- rep(1, m)
	if(silent == 0)
		print("VARLINK1: iter   gradient     betas          s2")
	for(iter in 1:maxiter) {
		u <- matrix(ncol = 2, nrow = length(unique(di)))
		i <- 1
		s2s <- rep(NA, length(unique(di)))
		u0=0
		for(id in unique(di)) {
			yi <- dat[di == id, 2]
			ni <- length(yi)
			sy <- sum(yi)
			Xi <- cbind(rep(1, ni), dat[di == id, 3:ncol(dat)])
			Bj <- exp(Xi %*% b)
			u <- li1(Bj = Bj, k = sy, u = u0, eps = 0.0001, MaxIter = 6)
			ei <- Bj * exp(u[1])
			d2i <- sum(ei/(1 + ei)^2)
			if(u[2] < 6)
				s2s[i] <- u[1]^2 - 1/d2i
			u0=u[1]
			i <- i + 1
		}
		s2 <- max(mean(s2s, na.rm = T), 0)
		b.old <- b
		for(i in 1:5) {
			Xb <- X %*% b
			ID <- twoprob(Xb, s2)
			d1 <- ((y - ID[, 1]) * ID[, 2])/ID[, 1]/(1 - ID[, 1])
			grad <- t(X) %*% d1
			d2 <- ID[, 2]^2/ID[, 1]/(1 - ID[, 1])
			XX <- X * (d2 %*% t(eem))
			Hess <- t(X) %*% XX
			iH <- solve(Hess)
			b <- b + iH %*% grad
		}
		delta <- b - b.old
		if(sum(abs(delta)) < eps)
			return(list(b, s2, iH, iter))
		if(silent == 0)
			print(c(iter, sum(sqrt(grad^2)), b, s2))
	}
	return(list(b, s2, iH, maxiter))
}
