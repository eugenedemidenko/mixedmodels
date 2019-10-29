INT.ch08 <-
function(beta=.5,sigma=.3,omega=2,K=13,nExp=500000)
{
	dump("INT.ch08","c:\\MixedModels\\Chapter08\\INT.ch08.r")
	par(mfrow=c(1,1),mar=c(4,4,1,1))
	tau=function(beta,yb,n,omega,eps=0.00001,maxit=100)
	{
		ny=length(y)
		if(yb<=0) return(rep(NA,ny))
		tau=ln(yb)-beta
		for(it in 1:maxit)
		{
			ex1=exp(beta+tau);ex2=ex1^2
			delta=(ex2-ex1*yb+tau/(n*omega^2))/(2*ex2-ex1*yb+1/(n*omega^2))
			tau=tau-delta
			if(max(abs(delta))<eps) return(tau)
		}		
		return(tau)
	}

	tau.solve=function(beta,yb,n,omega,eps,maxit)
	{
		ny=length(y)
		tau=log(yb)-beta
		tau[yb<=0]=0
		for(it in 1:maxit)
		{
			ex1=exp(beta+tau);ex2=ex1^2
			num=ex2-ex1*yb+tau/(n*omega^2)
			den=2*ex2-ex1*yb+1/(n*omega^2)
			delta=num/den
			tau=tau-delta
			if(max(abs(delta))<eps) break
		}
		return(tau)
	}

	Sbar=function(beta,beta0,n,sigma,omega,K=13)
	{
		wxk=gauher(K)
		x=rep(wxk[,1],times=K)
		y=rep(wxk[,1],each=K)
		wx=rep(wxk[,2],times=K)
		wy=rep(wxk[,2],each=K)
		yb=exp(beta0+sqrt(2)*sigma*omega*x)+sqrt(2/n)*sigma*y
		taus=tau.solve(beta=beta,yb=yb,n,omega,eps=0.00001,maxit=10)
		Sbar.val=2*sigma^2*omega/sqrt(n)*sum(taus*wx*wy,na.rm=T)
		ex1=exp(beta+taus);ex2=ex1^2
		t1=(ex1-2*ex2)/(2*ex2-ex1+1/n/omega^2)
		Sbar.der=2*sigma^2*omega/sqrt(n)*sum(t1*wx*wy,na.rm=T)
		return(c(Sbar.val,Sbar.der))
	}

	wxk=gauher(K)
	x=rep(wxk[,1],times=K)
	y=rep(wxk[,1],each=K)
	wx=rep(wxk[,2],times=K)
	wy=rep(wxk[,2],each=K)
	INT.qn=wx*wy
	ns=1:15
	exb=bs=rbsLB=rep(NA,15)
	for(n in ns)
	{
		INT.dom=exp(beta+sigma*omega*sqrt(2)*x)+sigma*sqrt(2)/sqrt(n)*y
		qn=sum(INT.qn[!is.na(INT.dom)])/pi
		INT.G=log(INT.dom)*wx*wy
		G=sum(INT.G[INT.dom>0])/pi
		exb[n]=100*(G/qn/beta-1)
		ai=rnorm(nExp,mean=beta,sd=sigma*omega)
		eps=matrix(rnorm(nExp*n,mean=0,sd=sigma),ncol=n,nrow=nExp)
		Y=exp(ai)%*%t(rep(1,n))+eps
		my=Y%*%rep(1/n,n)
		bs[n]=mean(log(my),na.rm=T)
		bas=beta
		for(it in 1:10)
		{
			Svd=Sbar(beta=bas,beta0=beta,n=n,sigma=sigma,omega=omega)
			delta=Svd[1]/Svd[2]
			bas=bas-delta
			if(abs(delta)<0.00001) break
		}
		rbsLB[n]=bas

	}
	matplot(ns,cbind(exb,-sigma^2/2/ns*exp(-2*beta)/beta*100,(bs-beta)/beta*100,(rbsLB-beta)/beta*100),type="b",xlab="Sample Size, n",ylab="% Relative Asymptotic Bias",col=1,lty=1)
	text(6,8,"1=Two-step exact bias\n2=Two-step approximate bias\n3=Simulation-based (nExp=1,000,000)\n4=LB exact bias",adj=0,cex=1.5)
	cbind(ns,exb,-sigma^2/2/ns*exp(-2*beta)/beta*100,(bs-beta)/beta*100,(rbsLB-beta)/beta*100)
}
