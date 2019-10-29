dupplus <-function(n)
{
D=dupp(n)
Dplus <- solve(t(D) %*% D) %*% t(D)
return(Dplus)
}
