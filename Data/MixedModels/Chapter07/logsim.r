logsim <-
function(N=100,n.min=10,n.max=15,sigma=2,nExp=10,ss=3)
{
	# Simulations for bias slope as function of SD intercept
    
	dump("logsim","c:\\MixedModels\\Chapter07\\logsim.r")
	print(date())  # it may take long time
	set.seed(ss)
	library(MASS) # for glmmPQL
	# besides it needss: logMLEgh, logMLE1, logFSL, logFS, logric (may be slow when n.max>15)
	ni <- round(runif(N,min=n.min,max=n.max))
	NT=sum(ni)
	beta1=1;beta2=1;beta3=-.05
	eps=0.001
	bt <- c(beta1, beta2,beta3)
	dat <- matrix(ncol = 4, nrow = sum(ni))
	X=matrix(0,nrow=NT,ncol=N)
	k=1
	for(i in 1:N) {
		xi <- rnorm(ni[i])
		dat[k:(k + ni[i] - 1), 1] <- i
		dat[k:(k + ni[i] - 1), 3] <- xi
		dat[k:(k + ni[i] - 1), 4] <- xi^2
        X[k:(k + ni[i] - 1),i]=1
		k <- k + ni[i]
	}
	id=dat[,1];x1=dat[,3];x2=dat[,4]
	b0=rep(0,3)            
	parML<-s2ML <- as.data.frame(matrix(nrow = nExp, ncol = 8))
	names(parML)=names(s2ML)=namMET=c("MLEgh13","MLEintegrate","FSgh13","FS500","logric","glmmPQL","FS1gh13","FS1500")

	for(isim in 1:nExp) {
		k <- 1
		ai <- sigma * rnorm(N)
		for(i in 1:N) {
			xi=dat[k:(k + ni[i] - 1), 3]
			lin <- ai[i]+beta1+beta2*xi+beta3*xi^2
			prob <- exp(lin)/(1 + exp(lin))
			yi <- rep(0, ni[i])
			yi[runif(ni[i]) < prob] <- 1
			dat[k:(k + ni[i] - 1), 2] <- yi
			k <- k + ni[i]
		}
        
        y=dat[,2]
        o=glm(y~x1+x2+X-1,family=binomial)
        a=coef(o)
        b0[1]=mean(a[3:(N+2)]);b0[2:3]=a[1:2]
        sigma20=var(a[3:(N+2)])
        
        o <- logMLEgh(dat=dat,b=b0,s2=sigma20,ngh=13,eps=eps)
	parML[isim, 1] <- (o[[1]])[2]
	s2ML[isim, 1] <- o[[2]]
		
	
	o <- logMLE1(dat = dat, b = bt, s2 = sigma20,eps=eps)
	parML[isim, 2] <- (o[[1]])[2]
	s2ML[isim, 2] <- o[[2]]
		
	o <- logFSL(dat = dat, bsigma = c(bt, sqrt(sigma20)), xwFSL = gauher(K = 13),eps=eps)
	parML[isim, 3] <- (o[[1]])[2]
	s2ML[isim, 3] <- o[[2]]
        
	o <- logFS(dat = dat, bsigma = c(bt, sqrt(sigma20)), nSim = 500,eps=eps)
	parML[isim, 4] <- (o[[1]])[2]
	s2ML[isim, 4] <- o[[2]]
        
	o <- logric(dat = dat, eps=eps)
	parML[isim, 5] <- (o[[1]])[1]
		        
	o <- glmmPQL(y~x1+x2, random=~1|id,family=binomial,verbose=F)
	if((attr(o,"class"))[1]!="try-error") 
	{
	        parML[isim, 6]=((o[[4]])$fixed)[2]
	       	s2ML[isim, 6]=as.numeric(as.matrix(VarCorr(o))[1,1])
	}
		
	o <- logFSL1(dat = dat, bsigma = c(bt, sqrt(sigma20)), xwFSL = gauher(K = 13),eps=eps)
	parML[isim, 7] <- (o[[1]])[2]
	s2ML[isim, 7] <- o[[2]]
		
	o <- logFS1(dat = dat, bsigma = c(bt, sqrt(sigma20)), nSim = 500,eps=eps)
	parML[isim, 8] <- (o[[1]])[2]
	s2ML[isim, 8] <- o[[2]]

	}

print("beta estimates:")
print(parML)
print("s2 estimates:")
print(s2ML)
    
mse=as.data.frame(matrix(ncol=2,nrow=8))
names(mse)=c("slope","s2")
row.names(mse)=namMET
for(i in 1:8)
{
    mse[i,1]=mean((parML[,i]-beta2)^2,na.rm=T)
    mse[i,2]=mean((s2ML[,i]-sigma^2)^2,na.rm=T)
}
print("MSE:")
print(mse)    
date()
}
