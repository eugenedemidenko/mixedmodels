poissGEEM <-
function(da, m, D, pr = 1, eps = 1.0000000000000001e-005, maxiter = 100)
{
	# It is assumed that the 1st col is id, 2nd y, from 3 to ncol(d) is X, no ones.
	# This is poissGEE function with a estimated by VLS
	di <- da[, 1]
	uid <- unique(di)
	N <- length(uid)
	y <- da[, 2]
	NT <- length(y)
	X <- da[, 3:(2 + m)]
	Z <- da[, (3 + m):ncol(da)]
	kz <- ncol(Z)
	DP <- dupp(kz)
	DPP <- dupplus(kz)
	X1 <- cbind(rep(1, NT), X)
	m <- ncol(X1)
	o <- glm(y ~ X, family = poisson())
	b <- o$coefficients
	k2 <- (kz * (kz + 1))/2
	if(pr == 1)
		print("GEEM")
	for(iter in 1:maxiter) {
		d <- DPP %*% as.vector(D)
		M1 <- matrix(0, k2, k2)
		m1 <- rep(0, k2)
		for(id in unique(di)) {
			yi <- y[di == id]
			ni <- length(yi)
			Xi <- X1[di == id,  ]
			Zi <- Z[di == id,  ]
			ex <- exp(Xi %*% b)
			V <- diag(as.vector(ex), ni, ni)
			ee <- yi - ex
			Ki <- matrix(ncol = k2, nrow = ni^2)
			kk <- 1
			for(j in 1:ni)
				for(k in 1:ni) {
					v1 <- exp(sum(b * (Xi[j,  ] + Xi[k,])))
					v2 <- exp(0.5 * t(Zi[j,  ] + Zi[k,]) %*% D %*% (Zi[j,  ] + Zi[k,  ]))
					v2 <- as.numeric(v2)
					v3 <- exp(0.5 * t(Z[j,  ]) %*% D %*% Z[	j,  ] + 0.5 * t(Z[k,  ]) %*% D %*%	Z[k,  ])
					v3 <- as.numeric(v3)
					V[k, j] <- V[k, j] + v1 * v2 - v1 *v3
					dV <- 0.5 * v1 * v2 * kronecker(Zi[j,  ] + Zi[k,  ], Zi[j,  ]+Zi[k,  ]) - 0.5 * v1 * v3 *(kronecker(Zi[j,  ], Zi[j,]) - kronecker(Zi[k,  ], Zi[	k,  ]))
					Ki[kk,  ] <- DPP %*% dV
					kk <- kk + 1
				}
			rr <- as.vector(ee %*% t(ee) - V)
			M1 <- M1 + t(Ki) %*% Ki
			m1 <- m1 + t(Ki) %*% rr
		}
		d <- d + solve(M1) %*% m1
		D <- matrix(DP %*% d, kz, kz)
		o <- eigen(D)
		PD <- o$vectors
		LD <- o$values
		LD[LD < 0] <- 0
		D <- PD %*% diag(LD, kz, kz) %*% t(PD)
		M1 <- matrix(0, m, m)
		m1 <- rep(0, m)
		for(id in unique(di)) {
			yi <- y[di == id]
			ni <- length(yi)
			Xi <- X1[di == id,  ]
			Zi <- Z[di == id,  ]
			ee <- exp(Xi %*% b)
			V <- diag(as.vector(ee), ni, ni)
			r <- yi - ee
			Ei <- diag(as.vector(ee), ni, ni)
			for(j in 1:ni)
				for(k in 1:ni) {
					v1 <- exp(sum(b * (Xi[j,  ] + Xi[k,])))
					v2 <- exp(0.5 * t(Zi[j,  ] + Zi[k,	]) %*% D %*% (Zi[j,  ] + Zi[k,  ]))
					v3 <- exp(0.5 * t(Z[j,  ]) %*% D %*% Z[	j,  ] + 0.5 * t(Z[k,  ]) %*% D %*%	Z[k,  ])
					V[k, j] <- V[k, j] + v1 * v2 - v1 *v3
				}
			iV <- solve(V)
			M1 <- M1 + t(Xi) %*% Ei %*% iV %*% Ei %*% Xi
			m1 <- m1 + t(Xi) %*% Ei %*% iV %*% r
		}
		b <- b + solve(M1) %*% m1
		print(b)
	}
	return(list(b, s2, solve(M1), maxiter))
}
