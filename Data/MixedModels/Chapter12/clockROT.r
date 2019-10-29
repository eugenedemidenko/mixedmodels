clockROT <-
function () 
{
    dump("clockROT", "c:\\MixedModels\\Chapter12\\clockROT.r")
    matR = function(M) # image reflection
	{
        nr = nrow(M)
        nc = ncol(M)
        MR = M
        for (i in 1:nc) MR[, nc - i + 1] = M[, i]
        return(MR)
    }
    matFILL = function(M) # filling NA breaks
	{
		MF=M
        nr = nrow(M)
        nc = ncol(M)
        for (i in 2:(nr - 1))
		for (j in 2:(nc - 1))
		if (is.na(M[i, j])) {
            mv = as.vector(M[(i - 1):(i + 1), (j - 1):(j + 1)])
            MF[i, j] = mean(mv, na.rm = T)
        }
        return(MF)
    }
    c1 <- scan("c:\\MixedModels\\Chapter12\\clock1.pgm", what = "")
    P <- as.numeric(c1[2])
    Q <- as.numeric(c1[3])
    M <- matrix(as.numeric(c1[5:length(c1)]), nrow = P, ncol = Q)
    par(mfrow = c(2, 2), mar = c(0, 0, 0, 0))
    thetaALL = c(0, pi/4, pi/2, -pi/6)
    thetaL = c("0", "p/4", "p/2", "-p/6")
    PQ = sqrt(2) * max(P, Q)
    for (irot in 1:4) {
        MROT = matrix(ncol = PQ, nrow = PQ)
        theta = thetaALL[irot]
        for (i in 1:P) for (j in 1:Q) {
            pp = (PQ/2) + (i - P/2) * cos(theta) + (j - Q/2) * 
                sin(theta)
            pp = round(pp)
            qq = (PQ/2) - (i - P/2) * sin(theta) + (j - Q/2) * 
                cos(theta)
            qq = round(qq)
            MROT[pp, qq] = M[i, j]
        }
        image(1:PQ, 1:PQ, matR(matFILL(MROT)), xlab = "", ylab = "", axes = F, col = gray(0:255/255))
        mtext(side = 3, paste("q =", thetaL[irot]), font = 5, cex = 1.5, line = -2)
    }
}
