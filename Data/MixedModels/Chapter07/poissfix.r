poissfix <-
function(d, pr = 1, eps = 1.0000000000000001e-005, maxiter = 100)
{
	# It is assumed that the 1st col is id, 2nd y, from 3 to ncol(d) is X, no ones.
	id <- d[, 1]
	uid <- unique(id)
	y <- d[, 2]
	m <- ncol(d) - 3 + 1
	X <- d[, 3:ncol(d)]
	q <- t(X) %*% y
	o <- glm(y ~ X, family = poisson())
	b <- o$coefficients
	b <- b[2:(m + 1)]
	if(pr == 1)
		print("Poisson fixed intercepts:")
	for(iter in 1:maxiter) {
		M1 <- matrix(0, m, m)
		sm <- rep(0, m)
		for(ii in uid) {
			yi <- y[id == ii]
			Xi <- X[id == ii,  ]
			ki <- sum(yi)
			bx <- Xi %*% b
			ei <- exp(bx)
			se <- sum(ei)
			ni <- length(yi)
			xe <- t(Xi) %*% ei
			M1 <- M1 + ((t(Xi) %*% diag(as.vector(ei), ni, ni) %*% Xi)/se - (xe %*% t(xe))/se^2) * ki
			sm <- sm + xe/se * ki
		}
		iM1 <- solve(M1)
		bnew <- b + iM1 %*% (q - sm)
		if(sum(abs(bnew - b)) < eps)
			return(list(b, 0, iM1, iter))
		if(pr)
			print(c(iter, bnew))
		b <- bnew
	}
	return(list(b, 0, iM1, maxiter))
}
