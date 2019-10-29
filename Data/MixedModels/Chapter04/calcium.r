calcium <-
function()
{
dump("calcium","c:\\MixedModels\\Chapter04\\calcium.r")
da=read.csv("c:\\MixedModels\\Chapter04\\calcium.txt")
par(mfrow=c(1,2))
grnam=c("Placebo Group","Calcium Group")
for(tre in c(1,0))
{
  dat=da[da$treat==tre,]
  plot(dat$time,dat$y,xlab="Weeks",ylab="Bone density, g/cm^2",xlim=c(0,120),ylim=c(.7,1.1))
  title(grnam[tre+1])
  uid=unique(dat$id)
  nuid=length(uid)
  for(j in 1:nuid)
    lines(dat$time[dat$id==uid[j]],dat$y[dat$id==uid[j]])
}
}
