GLSest <-
function(d,m,k,D)
#economical GLS estimate of beta
{
    id <- d[, 1]
    undi <- unique(id)
    y <- d[, 2]
    X <- as.matrix(d[, 3:(2 + m)], ncol = m)
    Z <- as.matrix(d[, (3 + m):(2 + m + k)], ncol = k)
    Xty=t(X)%*%y
    XtX=t(X)%*%X
    Ik=diag(rep(1,k),k,k)
    for(i in undi)
    {
        yi <- y[id == i]
        Xi <- X[id == i, ]
        Zi <- Z[id == i, ]
        ZtZ <- t(Zi)%*%Zi
        XtZ <- t(Xi)%*%Zi
        Zty <- t(Zi)%*%yi
        iM=solve(Ik+D%*%ZtZ)%*%D
        Xty=Xty-XtZ%*%iM%*%Zty
        XtX=XtX-XtZ%*%iM%*%t(XtZ)
    }
    beta.GLS=solve(XtX)%*%Xty
    return(beta.GLS)
}
