decE <-
function(A, nu, B, la, Cmat, maxiter=10, eps=10^-8)
{
    # Cmat is the 4x4 covariance matrix for A, nu, B, la
    t0 <- - A/nu
    for(iter in 1:maxiter) {
        exp1 <- exp(A + nu * t0)
        exp2 <- exp(B - la * t0)
        t1 <- t0 - (exp1 + exp2 - 1)/(nu * exp1 - la * exp2)
        if(abs(exp1 + exp2 - 1) < eps) break
        t0 <- t1
    }
    dFdA <- exp1
    dFdB <- exp2
    dFdnu <- t1 * exp1
    dFdla <- exp2 * t1
    dFdt <- nu * exp1 - exp2 * la
    h <- - c(dFdA, dFdnu, dFdB, dFdla)/dFdt #derivative vector
    var1 <- as.numeric(t(h) %*% Cmat %*% h) # delta-method
    return(c(t1, var1))
}
