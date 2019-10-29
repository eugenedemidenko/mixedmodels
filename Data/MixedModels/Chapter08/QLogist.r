QLogist <-
function(a1,a2,a3,a4,x)
{
    dene=a2-a3*x-a4*x^2 
    y=a1/(1+exp(dene))
    return(y)
}
