poissGEE <-
function(d, pr = 1, eps = 1.0000000000000001e-005, maxiter = 100)
{
	# It is assumed that the 1st col is id, 2nd y, from 3 to ncol(d) is X, no ones.
	# This is poissGEE function with a estimated by VLS
	di <- d[, 1]
	uid <- unique(di)
	N <- length(uid)
	y <- d[, 2]
	NT <- length(y)
	X <- d[, 3:ncol(d)]
	X1 <- cbind(rep(1, NT), X)
	m <- ncol(X1)
	o <- glm(y ~ X, family = poisson())
	b <- o$coefficients
	a <- 0
	if(pr == 1)
		print("GEE")
	for(iter in 1:maxiter) {
		M1 <- matrix(0, m, m)
		m1 <- rep(0, m)
		s1 <- s2 <- s3 <- 0
		for(id in unique(di)) {
			yi <- d[di == id, 2]
			ni <- length(yi)
			Xi <- cbind(rep(1, ni), d[di == id, 3:ncol(d)])
			exi <- exp(Xi %*% b)
			ri <- yi - exi
			Ei <- diag(as.vector(exi), ni, ni)
			Vi <- Ei + a * (exi %*% t(exi))
			iVi <- solve(Vi)
			Fi <- t(Xi) %*% Ei
			M1 <- M1 + Fi %*% iVi %*% t(Fi)
			m1 <- m1 + Fi %*% iVi %*% ri
			s1 <- s1 + sum(ri * exi)^2
			s2 <- s2 + sum(exi^3)
			s3 <- s3 + sum(exi^2)^2
		}
		delta <- solve(M1) %*% m1
		b <- b + delta
		anew <- min(max((s1 - s2)/s3, 0), 0.98999999999999999)
		s2 <-  - log(1 - anew)
		if(pr == 1)
			print(c(iter, b, s2))
		if(sum(abs(a - anew)) < eps) {
			b[1] <- b[1] - s2/2
			return(list(b, s2, solve(M1), iter))
		}
		a <- anew
	}
	return(list(b, s2, solve(M1), maxiter))
}
