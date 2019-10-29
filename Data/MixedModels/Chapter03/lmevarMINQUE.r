lmevarMINQUE <-
function(m, k, d)
{
# This function returns s2MINQQUE and DstarMINQUE for lme model
	dump("lmevarMINQUE", "c:\\MixedModels\\Chapter03\\lmevarMINQUE.r")
	options(width = 180)
	id <- d[, 1]
	idu <- unique(id)
	Mk <- matrix(0, k^2, k^2)
	yZ <- y <- d[, 2]
	Xz <- X <- d[, 3:(2 + m)]
	Nmat <- solve(t(X) %*% X)
	Xy <- t(X) %*% y
	bOLS <- Nmat %*% Xy
	kXy <- kronecker(Xy, Xy)
	yp <- X %*% bOLS
	N <- length(idu)
	Nt <- length(y)
	s2y <- sum((y - yp)^2)/(Nt - m)
	Nt2m <- Nt^2 - m
	D <- matrix(0, k, k)
	Qi <- matrix(0, k^2, k^2)
	qi <- rep(0, k^2)
	G <- matrix(0, ncol = k^2, nrow = m^2)
	Cv <- rep(0, k^2)
	for(i in idu) {
		Xi <- d[id == i, 3:(2 + m)]
		Zi <- d[id == i, (3 + m):(2 + m + k)]
		ZtX <- t(Xi) %*% Zi
		G <- G + kronecker(ZtX, ZtX)
	}
	for(i in idu) {
		yi <- d[id == i, 2]
		ni <- length(yi)
		Xi <- d[id == i, 3:(2 + m)]
		Zi <- d[id == i, (3 + m):(2 + m + k)]
		ZtZ <- t(Zi) %*% Zi
		ZtX <- t(Zi) %*% Xi
		ZtXN <- ZtX %*% Nmat
		ypi <- yp[id == i]
		pp <- as.vector(ZtZ - ZtX %*% Nmat %*% t(ZtX))
		Cv <- Cv + pp
		Qi <- Qi + kronecker(ZtZ, ZtZ) - kronecker(ZtXN, ZtXN) %*% G
		Zy <- t(Zi) %*% yi
		Zypi <- t(Zi) %*% ypi
		Xyi <- t(Xi) %*% yi
		qi <- qi + kronecker(Zy, Zy) - kronecker(Zypi, Zypi) - s2y * pp
		Pi <- Zi %*% ginverse.sym(ZtZ) %*% t(Zi)
		yZ[id == i] <- yi - Pi %*% yi
		Xz[id == i,  ] <- Xi - Pi %*% Xi
	}
	Qi <- Qi - Cv %*% t(Cv)/(Nt - m)
	DMINQUE <- solve(Qi) %*% qi
	DMINQUE <- matrix(DMINQUE, ncol = k, nrow = k)
	XzX <- t(Xz) %*% Xz
	if(sum(abs(XzX)) < 1e-010)
		bb <- matrix(0, m, m)
	else bb <- ginverse.sym(t(Xz) %*% Xz)
	rr <- t(Xz) %*% yZ
	yp <- Xz %*% bb %*% rr
	s2MINQUE <- sum((yZ - yp)^2)/(Nt - N * k)
	return(list(s2MINQUE, DMINQUE))
}
