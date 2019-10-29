asN <-
function(sigma=1,n=10,nExp=100)
{
dump("asN","c:\\MixedModels\\Chapter07\\asN.r")
# asymptotic properties of conditional 
# logistic regression
Ns=seq(from=100,to=1000,by=100)
NsN=length(Ns)
am=rep(NA,NsN)
a=rep(NA,nExp)    
beta1=1;beta2=-1;beta3=.5
    
for(isn in 1:NsN)
{
 NT=n*Ns[isn]
 ni=rep(n,Ns[isn])
 dat=matrix(nrow=NT,ncol=4)
 k=1
 for(i in 1:Ns[isn]) {
     xi <- rnorm(ni[i])
     dat[k:(k + ni[i] - 1), 1] <- i
     dat[k:(k + ni[i] - 1), 3] <- xi
     dat[k:(k + ni[i] - 1), 4] <- xi^2
     k <- k + ni[i]
 }
     
 for(isim in 1:nExp) {
     k <- 1
     ai <- sigma * rnorm(Ns[isn])
     for(i in 1:Ns[isn]) {
       xi=dat[k:(k + ni[i] - 1), 3]
       lin <- ai[i]+beta1+beta2*xi+beta3*xi^2
       prob <- exp(lin)/(1 + exp(lin))
       yi <- rep(0, ni[i])
       yi[runif(ni[i]) < prob] <- 1
       dat[k:(k + ni[i] - 1), 2] <- yi
       k <- k + ni[i]
 }
     
 o <- logric(dat = dat,silent=1)
 a[isim]=(o[[1]])[2] #beta3 estimate
 }
 am[isn]=mean(a) 
}
cbind(Ns,am)
}
