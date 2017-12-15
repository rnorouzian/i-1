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
  
#===============================================================================================
 
prop.priors <- function(a, b, lo = 0, hi = 1, dist.name, yes = 55, n = 1e2, scale = .1, top = 1.5, show.prior = FALSE){
  
d = dist.name
is.eq = function(...) length(unique(lengths(list(...)))) == 1
ab.df = function(x, y) c(x, rep(rev(x)[1], abs(length(y) - length(x))))

if(is.eq(a, b) & !is.eq(a, b, d)) d = ab.df(d, a) else if(is.eq(a, d) & !is.eq(a, b, d)) b = ab.df(b, a) else if(is.eq(b, d) & !is.eq(a, b, d)) a = ab.df(a, d)

deci = function(x, k = 3) format(round(x, k), nsmall = k)                                                                                                                           
  
  Bi = yes
  pr = show.prior
  loop = length(d)
  CI = matrix(NA, loop, 2)
  mode = numeric(loop)
  h = list()
  
  for(i in 1:loop){
    p = function(x) get(d[i])(x, a[i], b[i])*as.integer(x >= lo)*as.integer(x <= hi)
    prior = function(x) p(x)/integrate(p, lo, hi)[[1]]
    
    if(!pr){     
      likelihood = function(x) dbinom(Bi, n, x)
      k = integrate(function(x) prior(x)*likelihood(x), lo, hi)[[1]]
      posterior = function(x) prior(x)*likelihood(x) / k
      h[i] = list(curve(posterior, ty = "n", ann = FALSE, yaxt = "n", xaxt = "n", add = i!= 1, bty = "n", n = 1e3))
      mode[i] = optimize(posterior, c(lo, hi), maximum = TRUE)[[1]]
      CI[i,] = HDI(posterior)
    }
  }
  
  if(!pr){  
    
    plot(CI[, 1:2], rep(1:loop, 2), ty = "n", xlim = 0:1, ylim = c(1, top*loop), ylab = NA, yaxt = "n", xaxt = "n", xlab = "Credible Interval 'B'", font.lab = 2, mgp = c(2, .5, 0))
    abline(h = 1:loop, col = 8, lty = 3)
    axis(1, at = axTicks(1), lab = paste0(axTicks(1)*1e2, "%"), mgp = c(2, .3, 0))
    segments(CI[, 1], 1:loop, CI[, 2], 1:loop, lend = 1, lwd = 4, col = 1:loop, xpd = NA)
    axis(2, at = 1:loop, lab = substring(d, 2), font = 2, las = 1, cex.axis = .8, tck = -.006, mgp = c(2, .3, 0))
    legend("topleft", paste0(substring(d, 2), "(", round(a, 2), ", ", round(b, 2), ")"), pch = 22, pt.bg = 1:loop, col = 1:loop, cex = .7, bg = NA, bty = "n", pt.cex = .6, xpd = NA)
    for(i in 1:loop){
      polygon(x = h[[i]]$x, y = scale*h[[i]]$y +i, col = adjustcolor(i, .4), border = NA, xpd = NA)
    }
    m = scale*sapply(h, function(x) max(x[[2]])) + 1:loop
    segments(mode, 1:loop, mode, m, lty = 3, xpd = NA, lend = 1)  
    points(mode, 1:loop, pch = 21, bg = "cyan", cex = 1.3, col = "magenta", xpd = NA)
    I = deci(CI*1e2 , 2); o = deci(mode*1e2, 2)
    text(mode, 1:loop, paste0(I[,1], "%", "    ", o, "%", "    ", I[,2], "%"), cex = .75, pos = 3, font = 2, xpd = NA)
  }else{
    curve(prior, lo, hi, yaxt = "n", xaxt = "n", ylab = NA, xlab = "Proportion", bty = "n", font.lab = 2, lwd = 2, n = 1e3, main = bquote(Proportion*" ~ "*.(if(lo > 0 || hi < 1) "truncated-")*.(substring(d, 2))(.(round(a, 2)), .(round(b, 2)))))
    axis(1, at = axTicks(1), lab = paste0(axTicks(1)*1e2, "%"), mgp = c(2, .4, 0))
  }
  
}

#==========================================================================================================================
prop.hyper <- function(a, b, lo = 0, hi = 1, dist.name, yes = 55, n = 1e2, show.prior = FALSE, pos = 3){
  
  d = dist.name
  Bi = yes
  pr = show.prior
  d = dist.name
  d = dist.name
  is.eq = function(...) length(unique(lengths(list(...)))) == 1
  ab.df = function(x, y) c(x, rep(rev(x)[1], abs(length(y) - length(x))))
  
  if(is.eq(a, b) & !is.eq(a, b, d)) d = ab.df(d, a) else if(is.eq(a, d) & !is.eq(a, b, d)) b = ab.df(b, a) else if(is.eq(b, d) & !is.eq(a, b, d)) a = ab.df(a, d)
  
  
  deci = function(x, k = 3) format(round(x, k), nsmall = k)
  
  loop = length(a)
  CI = matrix(NA, loop, 2)
  mode = numeric(loop)
  
  for(i in 1:loop){
    p = function(x) get(d[i])(x, a[i], b[i])*as.integer(x >= lo)*as.integer(x <= hi)
    prior = function(x) p(x)/integrate(p, lo, hi)[[1]]
    if(!pr){
      like = function(x) dbinom(Bi, n, x)
      k = integrate(function(x) prior(x)*like(x), lo, hi)[[1]]
      posterior = function(x) prior(x)*like(x) / k
      mode[i] = optimize(posterior, c(lo, hi), maximum = TRUE)[[1]]
      CI[i,] = HDI(posterior)
    }
  }
  
  if(!pr){
    graphics.off()
    original.par = par(no.readonly = TRUE)
    on.exit(par(original.par))
    
    par(mgp = c(2.2, .3, 0), mar = c(5.1, 4.1, 4.1, 3))   
    plot(CI[, 1:2], rep(1:loop, 2), ty = "n", xlim = c(0, 1), ylim = c(1, 1.01*loop), ylab = NA, yaxt = "n", xaxt = "n", xlab = "Credible Interval 'Proportion'", font.lab = 2)
    abline(h = 1:loop, col = 8, lty = 3)
    axis(1, at = axTicks(1), lab = paste0(axTicks(1)*1e2, "%"))
    segments(CI[, 1], 1:loop, CI[, 2], 1:loop, lend = 1, col = "red4", xpd = NA)
    points(mode, 1:loop, pch = 21, bg = "red4", cex = .8, col = "red4", xpd = NA)
    axis(2, at = 1:length(a), lab = deci(a), font = 2, las = 1, cex.axis = .8, tick = FALSE, mgp = c(2, .2, 0))
    axis(4, at = 1:length(b), lab = deci(b), font = 2, las = 1, cex.axis = .8, tick = FALSE, mgp = c(2, .2, 0))
    text(par('usr')[1:2], par('usr')[4], c("A", "B"), pos = 3, cex = 1.5, xpd = NA, font = 2)
    I = deci(CI*1e2 , 2); o = deci(mode*1e2, 2)
    text(mode, 1:loop, paste0("[", I[,1], "%", ",  ", o, "%", ",  ", I[,2], "%", "]"), cex = .75, pos = pos, xpd = NA)
  }else{
    curve(prior, lo, hi, yaxt = "n", ylab = NA, xaxt = "n", xlab = "Proportion", bty = "n", font.lab = 2, lwd = 2, n = 1e3, main = bquote(Proportion*" ~ "*.(if(lo > 0 || hi < 1) "truncated-")*.(substring(d, 2))(.(round(a, 2)), .(round(b, 2)))))
    axis(1, at = axTicks(1), lab = paste0(axTicks(1)*1e2, "%"), mgp = c(2, .4, 0))
  }
}

#===================================================================================================================

ab.prop.hyper <- function(a, b, lo = 0, hi = 1, dist.name, add = FALSE, 
                          yes = 55, n = 1e2, col = 1, show.prior = FALSE){
  d = dist.name
  Bi = yes
  pr = show.prior    
  is.v = function(x) length(x) > 1
  d = dist.name
  is.eq = function(...) length(unique(lengths(list(...)))) == 1
  ab.df = function(x, y) c(x, rep(rev(x)[1], abs(length(y) - length(x))))
  
  if(is.eq(a, b) & !is.eq(a, b, d)) d = ab.df(d, a) else if(is.eq(a, d) & !is.eq(a, b, d)) b = ab.df(b, a) else if(is.eq(b, d) & !is.eq(a, b, d)) a = ab.df(a, d)
  
  if(is.v(a) & pr || is.v(b) & pr) message("\tNote: You can see only '1 prior' at a time.")
  if(add & pr) message("\tNote: 'add' only works for overlying 'Credible Intervals' to compare them.")
  
  deci = function(x, k = 3) format(round(x, k), nsmall = k)
  
  loop = length(d)
  CI = matrix(NA, loop, 2)
  mode = numeric(loop)
  
  for(i in 1:loop){
    p = function(x) get(d[i])(x, a[i], b[i])*as.integer(x >= lo)*as.integer(x <= hi)
    prior = function(x) p(x)/integrate(p, lo, hi)[[1]]
    
    if(!pr){    
      likelihood = function(x) dbinom(Bi, n, x)
      k = integrate(function(x) prior(x)*likelihood(x), lo, hi)[[1]]
      posterior = function(x) prior(x)*likelihood(x) / k
      mode[i] = optimize(posterior, c(lo, hi), maximum = TRUE, tol = 1e-20)[[1]]
      CI[i,] = HDI(posterior)
    }
  }
  
  if(!add & !pr){
    plot(rep(1:loop, 2), CI[, 1:2], type = "n", ylim = 0:1, xlim = c(1, loop), xlab = ifelse(d == "dbeta", expression(alpha), "Prior Parameter 'A'"), xaxt = "n", yaxt = "n", ylab = "Credible Interval (Proportion)", font.lab = 2, mgp = c(2.3, .3, 0), cex.lab = 1.2)
    abline(v = 1:loop, col = 8, lty = 3)
    axis(2, at = axTicks(2), lab = paste0(axTicks(2)*1e2, "%"), mgp = c(2, .4, 0), las = 1)
    axis(3, at = 1:length(b), lab = deci(b), font = 2, las = 1, cex.axis = .8, mgp = c(2, .2, 0))
    text(mean(par('usr')[1:2]), 1.06*par('usr')[4], ifelse(d == "dbeta", expression(beta), "Prior Parameter 'B'"), pos = 3, cex = 1.2, xpd = NA, font = 2)
    axis(1, at = 1:length(a), lab = round(a, 3), font = 2, las = 1, cex.axis = .8, mgp = c(2, .3, 0))
  }
  
  if(!pr){
    segments(1:loop, CI[, 1], 1:loop, CI[, 2], lend = 1, col = col)  
    lines(1:loop, mode, col = col, lty = 3)
    points(1:loop, mode, pch = 21, bg = col, cex = .8, col = col, xpd = NA)
  }
  
  if(!add & pr){curve(prior, lo, hi, yaxt = "n", ylab = NA, xaxt = "n", xlab = "Proportion", bty = "n", font.lab = 2, lwd = 2, n = 1e3, main = bquote(Proportion*" ~ "*.(if(lo > 0 || hi < 1) "truncated-")*.(substring(d, 2))(.(round(a, 2)), .(round(b, 2)))))
    axis(1, at = axTicks(1), lab = paste0(axTicks(1)*1e2, "%"), mgp = c(2, .4, 0))
  }
}

#====================================================================================================================

