lena <-
function()
{
dump("lena", "c:\\MixedModels\\Chapter12\\lena.r")
d <- scan("c:\\kluwer\\image\\lena.pgm",what="")
N <- 256
d=dr <- matrix(as.numeric(d[12:length(d)]), ncol = N, nrow = N)
for(i in 1:N) d[,i]=dr[,N-i+1]
dv <- as.vector(d)
ix <- 1:N
par(mfrow = c(2, 2), mar = c(1, 1, 2, 1))
image(ix, ix, d, xlab = "", ylab = "", axes = F,col=gray((0:255)/255))
mtext(side = 3, "Original image, EPP=7.5", cex = 1.25, line = 0.1)
md <- median(dv)
d1 <- d2 <- d3 <- d
d1[d < md] <- 0
d1[d >= md] <- 1
image(ix, ix, d1, xlab = "", ylab = "", axes = F,col=gray((0:255)/255))
mtext(side = 3, "Binary image, EPP=1", cex = 1.25, line = 0.1)
dr <- dv[order(dv)]
d2[d < dr[N^2/4]] <- 0
d2[d >= dr[N^2/4] & d < dr[N^2/2]] <- 1
d2[d >= dr[N^2/2] & d < dr[(2 * N^2)/3]] <- 3
d2[d >= dr[(2 * N^2)/3]] <- 4
image(ix, ix, d2, xlab = "", ylab = "", axes = F,col=gray((0:255)/255))
mtext(side = 3, "4 gray scale image, EPP=2", cex = 1.25, line = 0.1)
d3[d < dr[N^2/16]] <- 0
d3[d >= dr[(15 * N^2)/16]] <- 16
for(i in 1:14)
	d3[d >= dr[(i * N^2)/16] & d < dr[((i + 1) * N^2)/16]] <- i
image(ix, ix, d3, xlab = "", ylab = "", axes = F,col=gray((0:255)/255))
mtext(side = 3, "16 gray scale image, EPP=4", cex = 1.25, line = 0.1)
}
