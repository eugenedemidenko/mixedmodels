lmeD <-
function(N=100,m=3,rs=434,s2=.02)
{
    dump("lmeD","c:\\MixedModels\\Chapter03\\lmeD.r")
    library(nlme)
    s=sqrt(s2)
    set.seed(rs)
    k=m+1
    b=rnorm(m)
    ni=round(runif(n=N,min=80,max=100))
    NT=sum(ni)
    d=matrix(nrow=NT,ncol=2+m+k)
    D=matrix(rnorm(k^2),k,k)
    D=t(D)%*%D
    TD=t(chol(D))
    j <- 1
    for(i in 1:N) {
        n <- ni[i]
        d[j:(j + n - 1), 1] <- i
        Xi <- matrix(rnorm(n*m),ncol=m)
        d[j:(j + n - 1), 3:(2+m)] <- Xi
        Zi <- cbind(rep(1,n),Xi)        
        d[j:(j + n - 1), (3+m):(2+m+k)] <- Zi
        d[j:(j + n - 1), 2] <- Xi %*% b + Zi %*% TD %*% rnorm(k,0,s) + rnorm(n,0,s)
        j <- j + n
    }
    dL <- as.data.frame(d)
    names(dL)<-c("id","y","X1","X2","X3","Zint","Z1","Z2","Z3")
    print(date())
    o <- lme(fixed=y~X1+X2+X3,random=~Z1+Z2+Z3|id,data=dL,method="ML")
    print(date())
    print(summary(o))
    os <- matrix(as.numeric(VarCorr(o)),ncol=k+1)
    sdb=diag(os[1:k,2],k,k)
    LR=os[2:k,3:(k+1)]
    Db=matrix(ncol=k,nrow=k)
    R=diag(rep(1,k),ncol=k,nrow=k)
    for(i in 2:k)
        R[i,1:(i-1)]=R[1:(i-1),i]=LR[(i-1),1:(i-1)]
    Db=sdb%*%R%*%sdb
    print("Estimate of matrix D from lme:")
    print(Db/o$sigma^2)
    dFS=cbind(dL$id,dL$y,rep(1,NT),dL$X1,dL$X2,dL$X3,dL$Zint,dL$Z1,dL$Z2,dL$Z3)
    
    ofs <- lmeFS(d=dFS,m=4,k=4,D=matrix(0,4,4),MLRML="ML",MaxIter=25,epspar=0.0001,pr=F)
    print(ofs)
    print(date())
    ofM=lmevarMINQUE(d=dFS,m=4,k=4)
    print(ofM[[2]]/ofM[[1]])
    oD=lmevarMM(d=dFS,m=4,k=4,s2=ofM[[1]])
    print(oD/ofM[[1]])
    oD=lmevarUVLS(d=dFS,m=4,k=4)
    print("D UVLS estimate:")
    print(oD[[2]]/oD[[1]])
    print("Estimate of matrix D from lmeFS:")
    return(ofs[[9]])
}
