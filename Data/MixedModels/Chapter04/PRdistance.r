PRdistance <-
function()
{
 dump("PRdistance", "c:\\MixedModels\\Chapter04\\PRdistance.r")
 library(nlme)
 da=read.csv(file="c:\\MixedModels\\Chapter04\\PRdistance.txt")
 oAR<-lme(fixed=y~ti+sex,random=~1|id,correlation=corAR1(),data=da,method="ML")
 summary(oAR)
}
