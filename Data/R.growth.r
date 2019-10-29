`R.growth` <-
function(job = 1, regrtime=c(0,1,3,4), groupNames=c("Control","3 Gy Radiation","5Gy Radiation","10 Gy Radiation") )
{

	dat=scan("c:\\Projects\\Rgrowth\\tumdat.csv",sep=",",what="")  # read csv example data file

	datnumb=as.data.frame(matrix(as.numeric(dat[5:length(dat)]),ncol=4,byrow=T)) 
	names(datnumb)=c("Group","id","Time","TumVol")
	LNTumVol=log(datnumb$TumVol)
	datnumb=data.frame(datnumb,LNTumVol)
	nGroups=length(unique(datnumb$Group))
	TVrange=range(datnumb$LNTumVol)

	wid=.1 # defines the space between points for better visualization 


if(job==1) #Plots individual tumor volume data on the log scale for each group as a 2x2 panel
{
		par(mfrow=c(2,2)) # 2x2 panel assuming 4 groups
		for(ig in 1:nGroups)
		{
			dig=datnumb[datnumb$Group==(ig-1),]
			plot(dig$Time,dig$LNTumVol,xlab="Time after treatment",ylab="LN Tumor Volume",main=groupNames[ig],ylim=TVrange)
			if(ig>1) 
			{
				segments(regrtime[ig], TVrange[1],regrtime[ig], TVrange[2],col=3)
				text(regrtime[ig], TVrange[2],"Regrowth", adj=0,col=3)
			}

			idig=dig$id
			uidig=unique(idig)
			for(animal in uidig)
			lines(dig$Time[idig==animal],dig$LNTumVol[idig==animal])
			utim=unique(dig$Time)
			utim=utim[order(utim)]
			nutim=length(utim)
			mTV=seTV=rep(0,nutim)
			for(i in 1:nutim)
			{
				x=dig$LNTumVol[dig$Time==utim[i]]
				mTV[i]=mean(x)
				seTV[i]=sd(x)/sqrt(length(x))
				points(utim[i], mTV[i],pch=16,col=2,cex=1.5)
				segments(utim[i], mTV[i]-seTV[i],utim[i], mTV[i]+seTV[i],col=2)
				segments(utim[i]-wid, mTV[i]-seTV[i],utim[i]+wid, mTV[i]-seTV[i],col=2)
				segments(utim[i]-wid, mTV[i]+seTV[i],utim[i]+wid, mTV[i]+seTV[i],col=2)
			}
		lines(utim,mTV,col=2,lwd=3)
	}
}

if(job==2) #Plots mean group tumor volume data with SE
{
	par(mfrow=c(1,1))
	utim=datnumb$Time
	plot(utim,utim,ylim=TVrange,xlab="Time after treatment",ylab="Mean LN Tumor Volume +/-SE",type="n")
	for(ig in 1:nGroups)
	{
		dig=datnumb[datnumb$Group==(ig-1),]
		utim=unique(dig$Time)
		utim=utim[order(utim)]
		nutim=length(utim)
		mTV=seTV=rep(0,nutim)

		for(i in 1:nutim)
		{
			x=dig$LNTumVol[dig$Time==utim[i]]
			mTV[i]=mean(x)
			seTV[i]=sd(x)/sqrt(length(x))
		}
		yr=range(c(mTV-seTV,mTV+seTV))
		for(i in 1:nutim)
		{
			points(utim[i], mTV[i],pch=15+ig,col=1+ig,cex=1.5)
			segments(utim[i], mTV[i]-seTV[i],utim[i], mTV[i]+seTV[i],col=1+ig)
			segments(utim[i]-wid, mTV[i]-seTV[i],utim[i]+wid, mTV[i]-seTV[i],col=1+ig)
			segments(utim[i]-wid, mTV[i]+seTV[i],utim[i]+wid, mTV[i]+seTV[i],col=1+ig)
		}
		lines(utim,mTV,col=1+ig,lwd=3)
		if (ig>1)	segments(regrtime[ig], TVrange[1]-10, regrtime[ig], mTV[regrtime[ig]], col=1+ig)

		legend(min(utim),TVrange[2],groupNames,pch=16:(15+nGroups),lty=1,lwd=3,col=2:(1+nGroups),cex=1.5)

	}

}

if(job==3)
{
	datan=datnumb[datnumb$Group==0,]
	for(ig in 1:(nGroups-1))
		datan=rbind(datan,datnumb[datnumb$Group==ig & datnumb$Time>=regrtime[ig+1],])
	datan=rbind(datnumb[datnumb$Group>0 & datnumb$Time==0, ],datan)
	y=datan$LNTumVol
	ny=length(y)
	X=matrix(0,ncol=nGroups,nrow=ny)
	rate=datan$Time
	ind=datan$id
	k=1
	for(ig in 1:nGroups)
		X[datan$Group==(ig-1),ig]=1

	utim=unique(datnumb$Time)
	out=lme(y~X+rate-1,random=~1|ind)
	print(summary(out))
	par(mfrow=c(1,1))
	plot(utim,utim,ylim=TVrange,xlab="Time after treatment",ylab="Mixed model LN Tumor Volume",type="n")
	xt=c(min(utim),max(utim))
	a=out$coefficients$fixed
	for(ig in 1:nGroups)
	{
		lines(xt,a[ig]+a[nGroups+1]*xt,col=1+ig,lty=2)
		xtr=c(regrtime[ig],max(utim))
		lines(xtr,a[ig]+a[nGroups+1]*xtr,col=1+ig,lty=1,lwd=3)
	
	}


	outTable=as.data.frame(matrix(ncol=5,nrow=2*nGroups))
	names(outTable)=c("DT","TGD","LN SF","%SF","%KF")
	row.names(outTable)[seq(from=1,to=2*nGroups,by=2)]=groupNames
	row.names(outTable)[seq(from=2,to=2*nGroups,by=2)]=paste("SE ",groupNames,sep="")
	beta=a[nGroups+1]
	cov=out$varFix
	leg=rep("",nGroups)
	for(ig in 1:nGroups)
	{
		DT=(a[1]-a[ig]+log(2))/beta
		outTable[1+2*(ig-1),1]=DT
		df=rep(0,nGroups+1)
		df[1]=1/beta
		df[ig]=-1/beta
		df[nGroups+1]=-DT/beta
		se=sqrt(t(df)%*%cov%*%df)
		outTable[2*ig,1]=se

		TGD=(a[1]-a[ig])/beta
		outTable[1+2*(ig-1),2]=TGD
	
		df[nGroups+1]=-TGD/beta
		se=sqrt(t(df)%*%cov%*%df)
		outTable[2*ig,2]=se

		LNSF=a[ig]-a[1]
		SF=exp(LNSF)
		outTable[1+2*(ig-1),3]=LNSF


		df=rep(0,nGroups+1)
		df[1]=-1
		df[ig]=1
		se=sqrt(t(df)%*%cov%*%df)
		outTable[2*ig,3]=se


		outTable[1+2*(ig-1),4]=SF*100
		outTable[2*ig,4]=se*100*SF

		outTable[1+2*(ig-1),5]=(1-SF)*100
		outTable[2*ig,5]=se*100*SF
		leg[ig]=paste(groupNames[ig], " (%cell kill=",round((1-SF)*100),")",sep="")

	}
	legend(min(utim),TVrange[2],leg,lty=1,lwd=3,col=2:(1+nGroups))
	return(round(outTable,4))
}

}

