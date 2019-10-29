KSimageR <-
function(nu=1.4)
{
dump("KSimageR","c:\\MixedModels\\Chapter12\\KSimageR.r")
par(mfrow = c(1, 2), mar = c(1, 1, 3, 1), omi = c(0, 0, 0, 0))
d <- scan("c:\\MixedModels\\Chapter12\\grp11.pgm",what="")
d <- as.numeric(d[9:length(d)])
nr <- d[1];nc <- d[2]
d <- matrix(d[4:length(d)], nrow = nr, ncol = nc)
image(1:nr, 1:nc, d, xlab = "", ylab = "", axes = F,col=gray(0:255/255))
mtext(side = 3, "Control", line = 0.25, cex = 2)
dnu=floor((d/255)^nu*255)
image(1:nr, 1:nc, dnu, xlab = "", ylab = "", axes = F,col=gray(0:255/255))
mtext(side = 3, paste("nu =",nu), line = 0.25, cex = 2)
}
