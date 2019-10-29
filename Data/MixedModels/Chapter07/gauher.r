gauher <-
function(K, EPS = 3e-014, PIM4 = 0.75112554446494295, MAXIT = 100)
    {
        x <- w <- rep(0, K)
        m = (K + 1)/2
        for(i in 1:m) {
            if(i == 1) {
                z = sqrt(2 * K + 1) - 1.85575 * (2 * K + 1)^(-0.16667)
            }
            else if(i == 2) {
                z = z - (1.14 * K^0.426)/z
            }
            else if(i == 3) {
                z = 1.86 * z - 0.86 * x[1]
            }
            else if(i == 4) {
                z = 1.91 * z - 0.91 * x[2]
            }
            else {
                z = 2. * z - x[i - 2]
            }
            for(its in 1:MAXIT) {
                p1 = PIM4
                p2 = 0.
                for(j in 1:K) {
                    p3 = p2
                    p2 = p1
                    p1 = z * sqrt(2./j) * p2 - sqrt((j - 1)/j) * p3
                }
                pp = sqrt(2 * K) * p2
                z1 = z
                z = z1 - p1/pp
                if(abs(z - z1) <= EPS)
                    break
            }
            x[i] = z
            x[K + 1 - i] =  - z
            w[i] = 2./(pp * pp)
            w[K + 1 - i] = w[i]
        }
        return(cbind(x, w))
    }
