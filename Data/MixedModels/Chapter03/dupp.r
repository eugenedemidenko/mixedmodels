dupp <-function(n)
    {
       d1 <- n^2
       d2 <- (n * (n + 1))/2
       D <- matrix(rep(0, d1 * d2), nrow = d1)
       for(i in 1:n)
          for(j in 1:i) {
             u <- matrix(rep(0, d2), ncol = 1)
             u[(j - 1) * n + i - (j * (j - 1))/2, 1] <- 1
             Tm <- matrix(rep(0, d1), ncol = n)
             if(i == j)
                Tm[i, i] <- 1
             else Tm[i, j] <- Tm[j, i] <- 1
             D <- D + matrix(as.vector(Tm), ncol = 1) %*% t(u)
          }
       return(D)
    }
