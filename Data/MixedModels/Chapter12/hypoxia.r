hypoxia <-
function(job = 1)
{
	dump("hypoxia", "c:\\MixedModels\\Chapter12\\hypoxia.r")
	n <- 1024
	x <- 1:n
	if(job == 0) {
		h1 <- h2 <- matrix(0, ncol = 8, nrow = 256)
		cc <- "c:\\MixedModels\\Chapter\\Hypoxia\\Group1\\_"
		for(i in 1:8) {
			d <- scan(paste(cc, i, "c_1a_p.pgm", sep = ""),what="")
			d <- as.numeric(d[12:length(d)])
			for(g in 0:255) h1[g + 1, i] <- length(d[d == g])
		}
		cc <- "c:\\MixedModels\\Chapter12\\Hypoxia\\Group2\\_"
		for(i in 1:8) {
			d <- scan(paste(cc, i, "_1a_p.pgm", sep = ""),what="")
			d <- as.numeric(d[12:length(d)])
			for(g in 0:255) h2[g + 1, i] <- length(d[d == g])
		}
		#hypoxia.h_hypoxia(0)
		return(list(h1, h2))
	}
	if(job == 1) {
		par(mfrow = c(2, 8), mar = c(0, 1, 0, 0), omi = c(0, 0.25,
			0.25, 0))
		cc <- "c:\\MixedModels\\Chapter12\\Hypoxia\\Group1\\_"
		for(i in 1:8) {
			d <- scan(paste(cc, i, "c_1a_p.pgm", sep = ""),what="")
			d <- matrix(as.numeric(d[12:length(d)]), n, n)
			image(x, x, d, xlab = "", ylab = "", axes = F)
			mtext(side = 3, paste("Rat #", i, sep = ""), line =	0.25, cex = 1.25)
			if(i == 1)
				mtext(side = 2, "Control group", line = 1, cex = 1.5)
		}
		cc <- "c:\\MixedModels\\Chapetr12\\Hypoxia\\Group2\\_"
		for(i in 1:8) {
			d <- scan(paste(cc, i, "_1a_p.pgm", sep = ""),what="")
			d <- matrix(as.numeric(d[12:length(d)]), n, n)
			image(x, x, d, xlab = "", ylab = "", axes = F)
			if(i == 1)	mtext(side = 2, "Hypoxia group", line = 1, cex = 1.5)
		}
	}
	if(job == 2) {
		par(mfrow = c(2, 4), mar = c(4, 4, 3, 0), omi = c(0.1, 0.1, 0, 0))
		h1 <- hypoxia.h[[1]]
		h2 <- hypoxia.h[[2]]
		n2 <- n^2
		J <- n2/2
		for(i in 1:8) {
			f1 <- f2 <- rep(0, 256)
			for(g in 2:256) {
				f1[g] <- f1[g - 1] + h1[g, i]
				f2[g] <- f2[g - 1] + h2[g, i]
			}
			f1 <- c(f1/n2, 1)
			f2 <- c(f2/n2, 1)
			matplot(0:256, cbind(f1, f2), ylim = c(0, 1), type = 
				"l", xlab = "", ylab = "", axes = T, cex = 1)
			lines(0:256, f1, lwd = 5)
			lines(0:256, f2)
			mtext(side = 3, paste("Rat #", i, sep = ""), line = 
				0.75, cex = 1.25)
			mf <- max(abs(f1 - f2))
			lambda <- mf * (sqrt(J) + 0.11/sqrt(J) + 0.12)
			j <- 1:10000
			js <- rep(1, 10000)
			js[seq(from = 2, to = 10000, by = 2)] <- -1
			Q <- 2 * sum(js * exp(-2 * j^2 * lambda^2))
			print(c(i, Q))
			if(i == 1)
				legend(10, 0.20000000000000001, c("Control",
					"Hypoxia"), lty = 1, lwd = c(5, 1),
					cex = 0.80000000000000004)
		}
		mtext(side = 1, "Gray level", outer = T, cex = 1.5, line = -1
			)
		mtext(side = 2, "Cumulative probability", outer = T, cex = 1.5,
			line = -0.20000000000000001)
	}
	if(job == 3) {
		par(mfrow = c(1, 1), mar = c(4, 4, 0, 0), omi = c(0, 0, 0,
			0))
		rog <- matrix(nrow = 255, ncol = 2)
		B <- matrix(nrow = 8, ncol = 2)
		kk <- 1
		Nt <- 8 * n^2
		nn <- rep(0, 2)
		for(i in 1:2) {
			h <- hypoxia.h[[i]]
			k0 <- sum(h[1,  ])
			dg <- h[2:256,  ] %*% rep(1, 8)
			rog[, i] <- log(dg/(Nt - k0))
			nn[i] <- Nt - k0
			B[, i] <- (n^2 - h[1,  ])/h[1,  ]
		}
		matplot(kk:255, rog[kk:255,  ], type = "l", col = 1, xlab = "",
			ylab = "", axes = F)
		lines(kk:255, rog[kk:255, 1], lwd = 1)
		lines(kk:255, rog[kk:255, 2], lwd = 5)
		legend(50, -2, c("Control", "Hypoxia"), lty = 1, lwd = c(5,
			1), cex = 1.5)
		axis(side = 1, at = c(1, 50, 100, 150, 200, 255), cex = 1.25)
		axis(side = 2, at = c(-8, -7, -6, -5, -4, -3, -2), cex = 1.25,
			srt = 90)
		mtext(side = 1, "Gray level", cex = 2, outer = T, line = -1.5
			)
		mtext(side = 2, "r", cex = 2, font = 8, outer = T, line = -0.5
			)
		r1 <- exp(rog[2:255, 1]/(n^2 - nn[1]))
		r2 <- exp(rog[2:255, 2]/(n^2 - nn[2]))
		print(cbind(r1, r2))
		V1 <- (diag(r1, 254, 254) - r1 %*% t(r1))/(n^2 - nn[1])
		V2 <- (diag(r2, 254, 254) - r2 %*% t(r2))/(n^2 - nn[2])
		chi <- t(r1 - r2) %*% solve(V1 + V2) %*% (r1 - r2)
	}
}
