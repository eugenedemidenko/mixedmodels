logric <-
function(dat, maxiter = 20, eps = 1e-005, silent = 1)
{
    
    permniki=function(ni,ki)
    {
        # it returns all ni-dimensional 0/1 vectors
        # with the sum of elements=ki
        mall=expand.grid(rep(list(0:1),ni))   
        ms=rowSums(mall)
        mall=mall[ms==ki,]
        return(as.matrix(mall,ncol=ni))
    }
    
	y <- dat[, 2]
	X <- dat[, 3:ncol(dat)]
    r=t(X)%*%y    
	m <- ncol(dat) - 2
	if(m > 1)
		b <- glm(y ~ X)$coefficients[2:(ncol(X) + 1)]
	else b <- glm(y ~ X)$coefficients[2]
	di <- dat[, 1]
	b0 <- s2 <- NA
	for(iter in 1:maxiter) {
		grad <- rep(0, m)
		Hess <- matrix(0, m, m)
		loglik <- sum(b*r)
		for(id in unique(di)) {
			y <- dat[di == id, 2]
			ni <- length(y)
			ki <- sum(y)
			if(m == 1)
				X <- matrix(dat[di == id, 3], ncol = 1)
			else X <- dat[di == id, 3:ncol(dat)]
			tXy <- t(X) %*% y
			xb <- X %*% b
			Si <- 0
			gr <- rep(0, m)
			he <- matrix(0, m, m)
            perm=permniki(ni,ki)
            ri=as.vector(exp(perm%*%xb))
            Si=sum(ri)
            xz=perm%*%X
            gr=t(xz)%*%ri
            he=t(xz)%*%(ri*xz)
			if(Si > 0) {
				grad <- grad + tXy - gr/Si
				Hess <- Hess + he/Si - gr %*% t(gr)/Si^2
				loglik <- loglik - log(Si)
			}
		}
		iH <- solve(Hess)
		delta <- iH %*% grad
		if(sum(abs(delta)) < eps)
			return(list(b, iH, iter))
		b <- b + delta
		if(silent == 0)
			print(c(iter, loglik, sum(abs(grad)), b))
	}
    return(list(b, iH, maxiter))
}
