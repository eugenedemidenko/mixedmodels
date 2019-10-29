randmat <-
function (ss = 3) 
{
    dump("randmat", "c:\\MixedModels\\Chapter12\\randmat.r")
    gr = 0:255/255
    set.seed(ss)
    par(mfrow = c(6, 6), mar = c(0.1, 0.1, 0.1, 0.1), omi = c(0, 0, 0.2, 0))
    P <- 1
    M <- 100
    m <- 1:M
    em <- rep(1, M)
    OM <- matrix(c(1, 0.9, 0.9, 1), 2, 2)
    for (i in 1:36) {
        lambda <- runif(2)/3
        eta <- rnorm(2)/3
        A <- 2 + OM %*% rnorm(2)
        X1 <- lambda[1] * m %*% t(em) + eta[1] * em %*% t(m)
        X2 <- lambda[2] * m %*% t(em) + eta[2] * em %*% t(m)
        y <- A[1] * cos(X1) + A[2] * cos(X2) + 0.1 * matrix(rnorm(M^2), M, M)
        image(1:M, 1:M, y, axes = F, xlab = "", ylab = "", col = gray(gr))
        text(1, 0.9 * M, paste("l1=", round(lambda[1], 2),", l2=", round(lambda[2], 2), sep = ""), font = 5, col = 2, cex = 1, adj = 0)
        text(1, 0.8 * M, paste("h1=", round(eta[1], 2), ", h2=", round(eta[2], 2), sep = ""), font = 5, col = 2, cex = 1, adj = 0)
    }
}
