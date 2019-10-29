onexpSIM <-
function(nExp=100,beta=.5,n=10,N=100,method="gauher",sigma=.2,omega=2,K=13,ss=2)
{
dump(c("onexpSIM","onexpFSL"),"c:\\MixedModels\\Chapter08\\onexpSIMFSL.r")
set.seed(ss)
library(nlme)
d=as.data.frame(matrix(nrow=N*n,ncol=2))
names(d)=c("id","y")
bnlme=bFSL=matrix(NA,nrow=nExp,ncol=4)
print(date())
for(ix in 1:nExp)
{
	k=1
	for(i in 1:N)
	{
		d[k:(k+n-1),1]=i
		d[k:(k+n-1),2]=exp(rnorm(n=1,mean=beta,sd=sigma*omega))+rnorm(n,mean=0,sd=sigma)
		k=k+n
	}
	out=try(nlme(y~exp(b),fixed=b~1,random=b~1,groups=~id,data=d,start=c(b=beta)))
	if(attr(out,"class")[1]!="try-error")
	{
	    bnlme[ix,1]=out$coefficients$fixed
	    bnlme[ix,2:3]=as.numeric((as.matrix(VarCorr(out)))[,2])
	    bnlme[ix,4]=out$logLik
	}
    
	out=onexpFSL(d=d,N=N,n=n,beta=beta,sigma2=sigma^2,omega=omega*sigma,method=method,maxit=100,eps=1e-03,K=K,pr=F)
	if(out[5]<100)  bFSL[ix,]=c(out[2],out[3],out[1],out[4])
}
print(date())
MSE=as.data.frame(matrix(ncol=3,nrow=2))
names(MSE)=c("beta","sigma","sigma*omega")
row.names(MSE)=c("nlme","FSL")
partrue=c(beta,sigma,sigma*omega)
for(i in 1:3)
{
    MSE[1,i]=mean((bnlme[,i]-partrue[i])^2,na.rm=T)
    MSE[2,i]=mean((bFSL[,i]-partrue[i])^2,na.rm=T)
}
return(MSE)
}
onexpFSL <-
function (d, N, n, beta, sigma2, omega, method= "gauher", maxit = 100, eps = 1e-03, K = 131, pr=F) 
{
    dump("onexpFSL", "c:\\MixedModels\\Chapter08\\onexpFSL.r")
    if(method=="gauher")
    {
	    xw = gauher(K)
	    U = sqrt(2) * xw[, 1]
	    w = xw[, 2]
    }
    else 
    {
	w=rep(1/K,K)
	U=rnorm(K,mean=0,sd=sqrt(2))
    }		
    id = d[, 1]
    y = d[, 2]
    one = rep(1, n)
    oneN = rep(1, N)
    par = c(sigma2, beta, omega)
    for (it in 1:maxit) {
        iterdone=it
        loglik = -N * n*log(sqrt(2 * pi)) - N * n/2 * log(sigma2) + N*log(2)/2
        ders2 = derbeta = deromega = rep(0, N)
        for (i in 1:N) {
            yi = y[id == i]
            si2 = var(yi) * (n - 1)
            ybar = mean(yi)
            hk = exp(beta + U * omega)
            nyh = si2 + n * (ybar - hk)^2
            ek = exp(-nyh/2/sigma2)
            Swe = sum(w * ek)
            loglik = loglik + log(Swe)
            NUM = sum(w * ek * nyh)
            ders2[i] = NUM/Swe/2/sigma2^2 - n/2/sigma2
            derbeta[i] = sum(w * ek * hk * (ybar - hk))/Swe * n/sigma2
            deromega[i] = sum(w * ek * hk * (ybar - hk) * U)/Swe * n/sigma2
        }
        DER = cbind(ders2, derbeta, deromega)
        iH = solve(t(DER) %*% DER)
        grad = t(DER) %*% oneN
        parnew = par + iH %*% grad
        sigma2 = parnew[1]
        beta = parnew[2]
        omega = parnew[3]
        if (max(abs(parnew-par)) < eps) break
        if(pr) print(c(it,sqrt(par[1]),par[2:3], loglik, max(abs(grad))))
        par=parnew
    }
    # sigma,beta,omega,loglik
    return(c(sqrt(par[1]),par[2:3], loglik,iterdone))
}
