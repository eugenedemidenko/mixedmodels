lmeFS <-
function(d,m,k,D=matrix(0,k,k),MLRML="ML",MaxIter=25,epspar=0.0001,pr=F)
{
# This R function estimates LME model by ML/RML using FS algorithm
# Matrix D is free (can be NOT positive definite) 
# d is the data set
# The first column, d[, 1] must be id
# The second column, d[, 2] must be the response variable (y)
# m next columns, d[, 3:(2+m)] must be fixed effects variables (matrix X)
# k next columns, d[, (3+m):(2+m+k)] must be random effects variables (matrix Z)
# solcode=0: normal solution
# solcode=1: Iterations reached MaxIter
# solcode=2: min(eigen(D))<=0 at the final iteration

solcode=0 # default normal solution
id <- d[, 1]
undi <- unique(id)
N <- length(undi)
y <- d[, 2]
X <- as.matrix(d[, 3:(2 + m)], ncol = m)
Z <- as.matrix(d[, (3 + m):(2 + m + k)], ncol = k)
npar <- m + 1 + k * k # total number of parameters
Nt <- nrow(d)
cterm=-0.5*Nt*log(2*pi)
if(MLRML!="ML") cterm=cterm+0.5*m*log(2*pi)
b=GLSest(d,m,k,D) # compute GLS beta given D
ev=y-X%*%b
s2=sum(ev^2)
Ik=diag(rep(1,k),k,k)
# compute estimate sigma^2
for(i in undi)
{
   ei=ev[id == i]
   Zi=Z[id == i,  ]
   Zei=t(Zi)%*%ei
   ZtZ=t(Zi)%*%Zi
   IkZ=solve(Ik + D%*%ZtZ)%*%D
   s2=s2-sum(Zei*(IkZ%*%Zei))
}
s2=s2/Nt # s2 estimate
for(iter in 1:MaxIter)
{
  iterdone <- iter
  der <- rep(0, npar)
  H <- matrix(0, npar, npar)
  loglik <- cterm
  for(i in undi)
  {
     yi <- y[id == i]
     Xi <- X[id == i,  ]
     Zi <- Z[id == i,  ]
     ni <- length(yi)
     Ini <- diag(rep(1, ni), ni, ni)
     Vi <- Ini + Zi %*% D %*% t(Zi)
     iVi <- solve(Vi)
     ei <- yi - Xi %*% b
     Ri <- t(Zi) %*% iVi %*% Zi
     Pi <- t(Xi) %*% iVi %*% Zi
     ri <- t(Zi) %*% iVi %*% ei
     si <- t(Xi) %*% iVi %*% ei
     a <- ri %*% t(ri)
     a <- (kronecker(a, Ri) + kronecker(Ri, a))/s2
     der[1:m] <- der[1:m] + si
     ssi <- t(ei) %*% iVi %*% ei
     der[m + 1] <- der[m + 1] + ssi
     der[(m + 2):npar] <- der[(m + 2):npar] + as.vector(Ri - ri %*% t(ri)/s2)
     H[1:m, 1:m] <- H[1:m, 1:m] + t(Xi) %*% iVi %*% Xi/s2
     H[m + 1, m + 1] <- H[m + 1, m + 1] + (0.5 * ni)/s2^2
     H[m + 1, (m + 2):npar] <- H[m + 1, (m + 2):npar] + (0.5 * as.vector(Ri))/s2
     H[(m + 2):npar, (m + 2):npar] <- H[(m + 2):npar, (m + 2):npar] + 0.5 * kronecker(Ri, Ri)
     detVi=det(Vi)
     loglik <- loglik - 0.5 * log(abs(detVi)) - (0.5 * ssi)/s2 - 0.5 * ni * log(abs(s2))
  }
        
  der[1:m] <- der[1:m]/s2
  der[m + 1] <- (-0.5 * Nt)/s2 + (0.5 * der[m + 1])/s2^2
  der[(m + 2):npar] <- -0.5 * der[(m + 2):npar]
      
  if(MLRML!="ML")
  {
    delDer=matrix(0,ncol=k,nrow=k)
    iXX=ginverse.sym(H[1:m, 1:m]*s2)
    for(i in undi) 
    {
       Xi <- X[id == i,  ]
       Zi <- Z[id == i,  ]
       ni <- nrow(Xi)
       Ini <- diag(rep(1, ni), ni, ni)
       Vi <- Ini + Zi %*% D %*% t(Zi)
       iVi <- solve(Vi)
       ZVX=t(Zi)%*%iVi%*%Xi
       delDer=delDer+ZVX%*%iXX%*%t(ZVX)
    }
    der[(m + 2):npar]=der[(m + 2):npar]+0.5*as.vector(delDer)
    der[m + 1]=der[m + 1] + 0.5*m/s2   
    loglik=loglik-.5*log(abs(det(H[1:m,1:m])))
  }
        
  H[m + 1, 1:m] <- H[1:m, m + 1]
  H[(m + 2):npar, m + 1] <- H[m + 1, (m + 2):npar]
  H[(m + 2):npar, 1:m] <- H[1:m, (m + 2):npar]
        
  iH <- ginverse.sym(H)
  covbeta <- iH[1:m, 1:m]
  cov.s2vechD=iH[(m + 1):npar, (m + 1):npar]
        
  delta <- iH %*% der
        
  b=b+delta[1:m]
  s2=s2+delta[m+1]
  D<-D+matrix(delta[(m + 2):npar], ncol = k, nrow = k)
     
  if(pr)
  {
    o1 <- paste("================== ITERATION =", iter, " loglik =", loglik, ", grad =", sqrt(sum(der^2)), ", s2 =", s2)
    print(o1)
    print("beta-coefficients:")
    print(b)
    print("D-matrix:");print(D)
  }
        
  if(max(abs(delta))<epspar) break
} # maximum iterations
    
if(iterdone>=MaxIter) solcode=1	
if(min(eigen(D,only.values=T,symmetric=T)$values)<=0) solcode=2
grad=sqrt(sum(der^2))
retlist=list(solcode,iterdone,as.numeric(loglik),b,s2,covbeta,cov.s2vechD,grad,D)
names(retlist)=c("solcode","iterdone","loglik","b","s2","covbeta","cov.s2vechD","grad","D")
return(retlist)    
}
