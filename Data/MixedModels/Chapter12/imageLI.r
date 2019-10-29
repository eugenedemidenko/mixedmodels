imageLI <-
function() 
{
#Linear image interpolation
    dump("imageLI", "c:\\MixedModels\\Chapter12\\imageLI.r")
    
    par(mfrow = c(1, 2), mar = c(3, 2, 3, 1))
    Mor <- 4
    q <- p <- 1:Mor
    ee <- rep(1, Mor)
    M <- abs(p - Mor/2) %*% t(ee) + ee %*% t(abs(q - Mor/2))
	gr=0:255/255
    image(p, q, M, axes = F, xlab = "", ylab = "",col=gray(gr))
    mtext(side = 3, "Original 4 x 4 image", line = .5, cex = 1)
    axis(side = 1, at = 1:Mor, cex = 1.25)
    axis(side = 2, at = 1:Mor, srt = 90, cex = 1.25)
    N <- Mor * 4
    x <- y <- (1:N)/N
    MM <- matrix(0, N, N)
    for (i in 1:N) for (j in 1:N) {
        p <- min(max(floor(i/N * Mor), 1), Mor - 1)
        q <- min(max(floor(j/N * Mor), 1), Mor - 1)
        dx <- i/N * Mor - p;dy <- j/N * Mor - q
        if (dy > dx) 
            A <- (M[p + 1, q + 1] - M[p, q + 1]) * dx + (M[p,q + 1] - M[p, q]) * dy
            else A <- (M[p + 1, q] - M[p, q]) * dx + (M[p + 1,q + 1] - M[p + 1, q]) * dy
			
            MM[i, j] <- M[p, q] + A
        }
        image(1:N, 1:N, MM, axes = F, xlab = "", ylab = "",col=gray(gr))
        mtext(side = 3, "Linearly interpolated 16 x 16 image", line = .5, cex = 1)
        axis(side = 1, at = seq(from = 2, to = N, by = 2), cex = 1.25)
        axis(side = 2, at = seq(from = 2, to = N, by = 2), srt = 90,cex = 1.25)
}
