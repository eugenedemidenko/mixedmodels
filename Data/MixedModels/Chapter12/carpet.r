carpet <-
function()
{
dump("carpet", "c:\\MixedModels\\Chapter12\\carpet.r")
d <- scan("c:\\MixedModels\\Chapter12\\carpetc.pgm",what="")
d <- as.numeric(d[2:length(d)])
K <- d[1];J <- d[2]
y=d[4:length(d)]
carp.dat <- matrix(y, nrow = K, ncol = J)
x1=rep(1:K,times=J);x2=rep(1:J,each=K)
print("Linear model:")
o <- lm(y ~x1 + x2)
print(summary(o))
print("Nonlinear model:")
o <- nls(y ~a1 * sqrt(((a2 - x1)^2 + (x2 - a3)^2)),start = c(a1 = sqrt(0.006), a2 = 2000, a3 = -500))
print(summary(o))
a <- coef(o)
par(mfrow = c(1, 1), err = -1, mar = c(3, 3, 1, 2))
image(1:J, 1:K, t(carp.dat), axes = T, xlab = "", ylab = "",xlim = c(0, 1800), ylim = c(-500, J),col=gray(0:255/255))
points(a[2], a[3], pch = 16, cex = 1.5)
N <- 30;h <- 3000/N
theta <- seq(from = 0, to = 8 * atan(1), by = 0.01)
for(i in 2:(N - 1)) {
	dc <- h * i
     x <- dc * cos(theta)
     y <- dc * sin(theta)
     lines(a[2] + x, a[3] + y, lty = 2)
}
arrows(a[2], a[3], J, 0)
}
