lmevarMM <-
function(m, k, d, s2)
{
	# This function estimated Dstar by MM
    # s2 is the MINQUE, s2=lmeMINQUE(m,k,d)[[1]]
	dump("lmevarMM", "c:\\MixedModels\\Chapter03\\lmevarMM.r")
	id <- d[, 1]
	idu <- unique(id)
	N <- length(idu)
	X <- d[, 3:(2 + m)]
	Z <- d[, (3 + m):(2 + m + k)]
	y <- d[, 2]
	iN <- solve(t(X) %*% X)
	DbMMs <- Lmat <- Ebb <- matrix(0, k, k)
	Fmat <- matrix(0, k^2, k^2)
	M1 <- M2 <- matrix(0, k^2, m^2)
	bOLS <- iN %*% t(X) %*% y
	iXX <- matrix(0, m, m)
	for(i in idu) {
		yi <- y[id == i]
		nii <- length(yi)
		Ini <- diag(rep(1, nii), ncol = nii, nrow = nii)
		Xi <- X[id == i,  ]
		Zi <- Z[id == i,  ]
		ZtZ <- t(Z) %*% Z
		iZtZ <- ginverse.sym(ZtZ)
		Zp <- iZtZ %*% t(Zi)
		Ji <- Zp %*% Zi
		XtZ <- t(Xi) %*% Zi
		iXX <- iXX + ginverse.sym(t(Xi) %*% Xi)
		ZpX <- Zp %*% Xi
		M1 <- M1 + kronecker(ZpX, ZpX)
		M2 <- M2 + kronecker(XtZ, XtZ)
		XnX <- Xi %*% iN %*% t(Xi)
		Rii <- Zp %*% XnX %*% Zi
		tPi <- Ji %*% Rii
		Fmat <- Fmat + kronecker(Ji, Ji) - kronecker(Ji, tPi) - kronecker(tPi, Ji)
		ei <- yi - Xi %*% bOLS
		bi <- Zp %*% ei
		bib <- bi %*% t(bi)
		Ebb <- Ebb + bib
		Lmat <- Lmat + bib - s2 * (Zp %*% (Ini - XnX) %*% t(Zp))
		DbMMs <- DbMMs + (bib - s2 * iZtZ)/N
	}
	Fmat <- Fmat + M1 %*% kronecker(iN, iN) %*% M2
	DstarMM <- matrix(solve(Fmat) %*% as.vector(Lmat), ncol = k, nrow = k)
	return(DstarMM)
}
