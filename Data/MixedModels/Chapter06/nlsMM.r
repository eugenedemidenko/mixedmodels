nlsMM <-
function(b1=1,b2=2,n=10,sigma=.5,theta=.3,nExp=1000)
{
 dump("nlsMM","c:\\MixedModels\\Chapter06\\nlsMM.r")
 MMw=function(x,a1,a2,chW) # weighted Michaelis-Menten model
 {
    fx=a1*x/(a2+x)
    return(chW%*%fx)
 }
 x=1:n; un=rep(1,n)
 W=diag(rep(1,n),n,n)+theta*un%*%t(un)
 iW=solve(W); chW=chol(W) # Cholesky decomposition
 apar=matrix(ncol=2,nrow=nExp)
 for(iexp in 1:nExp)
 {
   y=b1*x/(b2+x)+rnorm(n,mean=0,sd=sigma)
   z=-y/x
   olm=lm(y~z) # starting values
   yc=chW%*%y # transform y vector
   onls <- try(nls(yc~MMw(x=x,b1,b2,chW=chW),
            start=list(b1=coef(olm)[1],b2=coef(olm)[2]),
            control = list(maxiter = 500)))
  if(attr(onls,"class")!="try-error") 
        apar[iexp,]=summary(onls)$coefficients[,1]
 }
 return(apar)
}
