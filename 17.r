beta.id = Vectorize(function(Low, High, Cover = NA){
  
options(warn = -1)
L <- if(is.character(Low)) as.numeric(substr(Low, 1, nchar(Low)-1)) / 100 else Low
U <- if(is.character(High)) as.numeric(substr(High, 1, nchar(High)-1)) / 100 else High
  
if(L <= 0 || U >= 1) stop("NOTE: The smallest LOWER value that you can choose is \".000001\"AND the largest UPPER value is \".999999\".")
if(L >= U) stop("Put the smaller value for Low, and the larger value for High")
  
coverage  <- if(is.character(Cover)) as.numeric(substr(Cover, 1, nchar(Cover)-1)) / 100 else if(is.na(Cover)) .95 else Cover
  
p1 = (1 - coverage) / 2 
p2 = 1 - p1
  
if( p1 <= 0 || p2 >= 1 || Low > High || p1 > p2 || coverage >= 1 ){
stop("Error: \n\tUnable to find such a prior, make sure you have selected the correct values.") 
    } else {

f.beta <- function(alpha, beta, x, lower = 0, upper = 1){
p <- pbeta((x-lower)/(upper-lower), alpha, beta)
      log(p/(1-p))
    }

delta <- function(fit, actual) sum((fit-actual)^2)
 
objective <- function(theta, x, prob, ...) {
      ab <- exp(theta)
      fit <- f.beta(ab[1], ab[2], x, ...)
      return (delta(fit, prob))
    }
    
x.p <- (function(p) log(p/(1-p)))(c(p1, p2))
    
sol <- nlm(objective, log(c(1e1, 1e1)), x = c(L, U), prob = x.p, lower = 0, upper = 1, typsize = c(1, 1), 
           fscale = 1e-12, gradtol = 1e-12)
    
parms <- as.numeric(exp(sol$estimate))
    
quantiles <- qbeta(p = c(p1, p2), parms[1], parms[2])
    
unequal <- function(a, b, sig = 3) round(a, sig) != round(b, sig)
    
if(unequal(L, quantiles[1]) || unequal(U, quantiles[2])){
      
stop("Error: \n\tUnable to find such a prior, make sure you have selected the correct values.")

  }else{
      
return(c(alpha = parms[[1]], beta = parms[[2]]))    
    }
  }
  
}, c("Low", "High", "Cover"))