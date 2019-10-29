twoint <-
function (K=13)
    {
    dump("twoint","c:\\MixedModels\\Chapter07\\twoint.r")
    xw=gauher(K)
    x=rep(xw[, 1],times=K)
    y=rep(xw[, 1],each=K)
    wx=rep(xw[, 2],times=K)
    wy=rep(xw[, 2],each=K)
    INT=sum((x^2+y^2)*wx*wy)
    INT 
    }
