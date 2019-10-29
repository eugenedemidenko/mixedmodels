callmetaMLFS <-
function(n=10,beta=1,sigma2=1,dfc=3,nExp=100)
{
    dump(c("callmetaMLFS","metaMLFS"),"c:\\Kluwer\\metaMLFS.r")
    sigma2i=rchisq(n,df=dfc)
    for(iexp in 1:nExp)
    {
        y=beta+rnorm(n,mean=0,sd=sqrt(sigma2+sigma2i))
      
        betas2=metaMLFS(y,sigma2i)
        print(betas2)
        
        
        
    }
    
}
metaMLFS <-
function(y, sigma2i,maxiter=10,eps=10^-7)
{
    #Fisher scoring algorithm for simple meta-anaysis model
    n = length(y)
    w <- 1/sigma2i
    sw <- sum(w)
    beta0 <- sum(y * w)/sw
    sigma2ML <- 0
    for(iter in 1:maxiter)
    {
        w <- 1/(sigma2i + sigma2ML)
        sw <- sum(w)
        sw2 <- sum(w^2)
        betaML <- sum(y * w)/sw
        sigma2ML.new <- sum(((y - betaML)^2 - sigma2ML) * w^2)/sw2
        if(abs(sigma2ML - sigma2ML.new) < eps) break
        sigma2ML=sigma2ML.new
    }
    return(c(betaML, sigma2ML))
}
