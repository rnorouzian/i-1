prob_ab <- function(fun, a, b, domain){
  total.area <- integrate(fun, domain[1], domain[2])[[1]]
  integrate(fun, a, b)[[1]] / total.area
}
invert_prob_ab <- function(fun, a, prob, domain){  
  O <- function(b, fun, a, prob){
  (prob_ab(fun, a, b, domain = domain) - prob)^2
}
  b <- optimize(O, c(a, domain[2]), a = a, fun = fun, prob = prob)$minimum
  return(b)
}

HDI <- function(fun, prob = .95, domain = c(0, 1)){
  mode <- optimize(fun, interval = domain, maximum = TRUE, tol = 1e-12)[[1]]
     O <- function(a, fun, prob, domain){
     b <- invert_prob_ab(fun, a, prob, domain)  
     b - a
}
 abest <- optimize(O, c(0, mode), fun = fun, prob = prob, domain = domain)$minimum
     b <- invert_prob_ab(fun, abest, prob, domain) 
     return(c(abest,b))
}