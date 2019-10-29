poissGEE1 <-
function(d, pr = 1, eps = 1.0000000000000001e-005, maxiter = 100)
{
	# It is assumed that the 1st col is id, 2nd y, from 3 to ncol(d) is X, no ones.
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
		print(paste("GEE1"))
	for(iter in 1:maxiter) {
		M1 <- matrix(0, m, m)
		m1 <- rep(0, m)
		s1 <- s2 <- 0
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
			s1 <- s1 + sum(ri * (iVi %*% ri))
			s2 <- s2 + sum(ri * (iVi %*% exi))^2
		}
		anew <- a + (s1 - (NT - m))/s2
		anew <- min(max(0, anew), 0.98999999999999999)
		delta <- solve(M1) %*% m1
		b <- b + delta
		s2 <-  - log(1 - a)
		if(pr == 1)
			print(c(iter, b, a))
		if(sum(abs(a - anew)) < eps) {
			b[1] <- b[1] - s2/2
			return(list(b, s2, solve(M1), iter))
		}
		a <- anew
	}
	return(list(b, s2, solve(M1), maxiter))
}
