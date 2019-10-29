poissHeck <-
function(dat, s2, maxiter = 20, eps = 0.0001, silent = 1)
{
	hi <- function(u, s2, B, ki)
	{
		r <- ki * u - (0.5 * u^2)/s2 - exp(u) * B
		return(exp(r))
	}
	u2hi <- function(u, s2, B, ki)
	{
		r <- ki * u - (0.5 * u^2)/s2 - exp(u) * B
		return(exp(r) * u^2)
	}
	y <- dat[, 2]
	X <- dat[, 3:ncol(dat)]
	b <- glm(y ~ X, family = poisson())$coefficients
	if(silent != 1) {
		print(paste("Heckman method, varb0=", s2, "betaPoisson:"))
		print(b)
	}
	X <- cbind(rep(1, nrow(dat)), X)
	q <- t(X) %*% y
	m <- length(b)
	di <- dat[, 1]
	N <- length(unique(di))
	dd <- rep(0, m + 1)
	h <- 0.001
	NN <- 20
	LL <- s2s <- seq(from = s2/2, to = s2 * 2, length = NN)
	for(iter in 1:maxiter) {
		B <- 0
		A <-  - sum(y)/N
		us <- seq(from = -3 * sqrt(s2), to = 3 * sqrt(s2), by = h)
		cc <- matrix(0, m, m)
		for(id in unique(di)) {
			yi <- dat[di == id, 2]
			sy <- sum(yi)
			ni <- length(yi)
			Xi <- cbind(rep(1, ni), dat[di == id, 3:ncol(dat)])
			Bi <- sum(exp(Xi %*% b)) * exp( - s2/2)
			Ii <- sum(hi(u = us, s2 = s2, B = Bi, ki = sy)) * h
			Ki <- sum(u2hi(u = us, s2 = s2, B = Bi, ki = sy)) * h
			Mi <- sum(hi(u = us, s2 = s2, B = Bi, ki = sy + 1)) * h
			B <- B + Ki/Ii/N
			A <- A + (Mi * Bi)/N/Ii
			eei <- as.vector(exp(Xi %*% b - s2/2))
			Vi <- diag(as.vector(eei), ni, ni) + (1 - exp( - s2)) * eei %*% t(eei)
			cc <- cc + t(Xi) %*% solve(Vi) %*% Xi
		}
		s2new <- (0.5 * (1 - sqrt(1 - 4 * A * B)))/A
		if(silent != 1)
			print(paste("iter=", iter, "varb=", s2new))
		if(abs(s2 - s2new) < eps) {
			b[1] <- b[1] - s2new/2
			return(list(b, s2new, solve(cc), iter))
		}
		s2 <- s2new
	}
	return(list(b, s2, maxiter))
}
