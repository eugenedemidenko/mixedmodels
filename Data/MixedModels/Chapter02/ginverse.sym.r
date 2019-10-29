ginverse.sym <-
function(A,eps=10^(-8))
{
    #Generalized inverse of a symmetric matrix A
    PV=eigen(A,symmetric=T)
    V0=IV=PV$values
    IV[abs(V0)>eps]=1/V0[abs(V0)>eps]
    IV[abs(V0)<=eps]=0
    Ainv=PV$vectors%*%(IV*(t(PV$vectors)))
    return(Ainv)
}
