RobustMedianML <-
function(b, y, si, nu, maxiter = 100, eps = 0.0001)
{
    n <- length(y)
    ssi2 <- sum(si^2)
    FIS <- matrix(0, 2, 2)
    der <- rep(0, 2)
    lambda <- 1
    for(iter in 1:maxiter) {
        Ei <- exp(nu * (y - b))
        Fi1 <- pnorm( - nu * si - (y - b)/si)
        fi1 <- dnorm( - nu * si - (y - b)/si)
        Fi2 <- pnorm( - nu * si + (y - b)/si)
        fi2 <- dnorm( - nu * si + (y - b)/si)
        asu <- Ei * Fi1 + Fi2/Ei
        loglik <- n * log(nu) + 0.5 * nu^2 * ssi2 + sum(log(asu))
        derb <- (fi1/si * Ei - nu * Ei * Fi1 + nu/Ei * Fi2 - fi2/si/Ei)/asu
        dernu <- 1/nu + nu*si^2-(si*fi1*Ei-(y-b)*Ei*Fi1+(y-b)/Ei*Fi2+si/Ei*fi2)/asu
        der[1] <- sum(derb)
        der[2] <- sum(dernu)
        FIS[1, 1] <- sum(derb^2)
        FIS[1, 2] <- FIS[2, 1] <- sum(derb * dernu)
        FIS[2, 2] <- sum(dernu^2)
        iFIS <- solve(FIS)
        delta <- iFIS %*% der
        lambda <- 1
        for(try in 1:10) {
            bnu <- c(b, nu) + lambda * delta
            b1 <- bnu[1]
            nu1 <- bnu[2]
            Ei <- exp(nu1 * (y - b1))
            Fi1 <- pnorm( - nu1 * si - (y - b1)/si)
            fi1 <- dnorm( - nu1 * si - (y - b1)/si)
            Fi2 <- pnorm( - nu1 * si + (y - b1)/si)
            fi2 <- dnorm( - nu1 * si + (y - b1)/si)
            asu <- Ei * Fi1 + Fi2/Ei
            loglik1 <- n * log(nu1) + 0.5 * nu1^2 * ssi2 + sum(log(asu))
            if(loglik1 > loglik)  break
            lambda <- lambda/2
        }
        bnew <- b1; nunew <- nu1
        if(abs(b - bnew) + abs(nu - nunew) < eps)
            return(cbind(c(b, nu), sqrt(diag(iFIS))))
        b <- bnew
        nu <- nunew
    }
}
