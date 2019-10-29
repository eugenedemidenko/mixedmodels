lmevarUVLS <-
function(m, k, d)
{
# This function computes VLSU estimates for s2 and Dstar
id <- d[, 1]
idu <- unique(id)
X <- d[, 3:(2 + m)]
Z <- d[, (3 + m):(2 + m + k)]
y <- d[, 2]
k2 <- k^2
k21 <- k2 + 1
Fmat <- matrix(0, k21, k21)
vf <- rep(0, k21)
tXy <- t(X) %*% y
iN <- solve(t(X) %*% X)
bOLS <- iN %*% tXy
Nt <- length(y)
Gmat <- matrix(0, k * m, k * m)
Hmat <- matrix(0, k^2, k^2)
for(i in idu) 
{
   yi <- y[id == i]
   nii <- length(yi)
   Xi <- X[id == i,  ]
   Zi <- Z[id == i,  ]
   ei <- yi - Xi %*% bOLS
   tZie <- t(Zi) %*% ei
   ZtZ <- t(Zi) %*% Zi
   tZX <- t(Zi) %*% Xi
   Gmat <- Gmat + kronecker(tZX, tZX)
   Pi <- Xi %*% iN %*% t(Xi)
   tZPZ <- t(Zi) %*% Pi %*% Zi
   Hmat <-Hmat+kronecker(ZtZ,ZtZ)-kronecker(ZtZ,tZPZ)-kronecker(tZPZ,ZtZ)
   Fmat[2:k21, 1] <- Fmat[2:k21, 1] + as.vector(ZtZ - tZPZ)
   vf[2:k21] <- vf[2:k21] + kronecker(tZie, tZie)
   vf[1] <- vf[1] + sum(ei^2)
}
Fmat[1, 1] <- Nt - m
Fmat[2:k21, 2:k21] <- Hmat + Gmat %*% kronecker(iN, iN) %*% t(Gmat)
Fmat[1, 2:k21] <- Fmat[2:k21, 1]
ds2 <- solve(Fmat) %*% vf
s2VLSU <- ds2[1]
DVLSU <- matrix(ds2[2:k21], nrow = k, ncol = k)
return(list(s2VLSU, DVLSU))
}
