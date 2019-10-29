maple <-
function()
{
dump("maple","c:\\MixedModels\\Chapter11\\maple.r")    
par(mfrow = c(3, 3), mar = c(0, 0, 0, 0))
for(i in 1:9) {
    pf <- paste("c:\\MixedModels\\Chapter11\\maple", as.character(i), ".xy", sep = "")
    d <- matrix(scan(pf), ncol = 2)
    d <- d[!is.na(d[, 1]) & !is.na(d[, 2]),  ]
    n <- nrow(d)
    dy <- min(d[, 2])
    im <- (1:n)[d[, 2] == dy]
    plot(d[, 1],  - d[, 2], type = "l", xlab = "", ylab = "", axes = F)
    points(d[im, 1],  - d[im, 2], pch = 1, cex = 1.5,col=2)
    polygon(x = d[, 1], y =  - d[, 2], col = 3)
    lines(d[, 1],  - d[, 2],lwd=4)
}
}
