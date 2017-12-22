HDI <- function(FUN, lower = 0, upper = 1, level = .95, eps = 1e-3){
  
  if(!is.function(FUN)) stop("Error: 'FUN' must be a function.")
  if(length(formals(FUN)) > 1) stop("Error: 'FUN' must be a 'single-argument' function.")
  if(1 <= level || level <= 0) stop("Error: 'level' must be between '0' and '1'.")
  x <- formals(FUN)
  fun <- function(x) FUN(x)

  lower = min(lower, upper) ; upper = max(lower, upper)
  posterior = function(x) fun(x)/integrate(fun, lower, upper)[[1]]
  mode = optimize(posterior, c(lower, upper), maximum = TRUE, tol = 1e-12)[[1]]
  inverse.posterior <- function(x, side = "left") {
    target <- function(y) posterior(y) - x
    ur <- switch(side,
                 left = try(uniroot(target, interval = c(lower, mode))),
                right = try(uniroot(target, interval = c(mode, upper))))
    if(inherits(ur, "try-error")) stop("inverse.posterior failed: You may change prior hyperparameters or extend limit?")
    return(ur$root)
  }
  areafun <- function(h) {
    i1 <- inverse.posterior(h, "left")
    i2 <- inverse.posterior(h, "right")
    return(integrate(posterior, i1, i2)[[1]])
  }
  post.area <- 1
  find.lims <- function(a) {
    ur <- uniroot(function(h) areafun(h) / post.area - a,
                  c(eps, posterior(mode) - eps))
    return(ur$root)
  }
  f <- find.lims(level)
  return(c(inverse.posterior(f, "left"),
           inverse.posterior(f, "right")))
}

#==================================================================================================================

hdi <- function(x, y, level = .95){
if(1 <= level || level <= 0) stop("Error: 'level' must be between '0' and '1'.")
  dx <- diff(x)
  areas <- dx * .5 * (head(y, -1) + tail(y, -1))
  peak <- which.max(areas)
  range <- c(peak, peak)
  found <- areas[peak]
  while(found < level) {
    if(areas[range[1]-1] > areas[range[2]+1]) {
      range[1] <- range[1]-1
      found <- found + areas[range[1]-1]
    } else {
      range[2] <- range[2]+1
      found <- found + areas[range[2]+1]
    }
  }
  val<-x[range]
  attr(val, "indexes") <- range
  attr(val, "area") <- found
  return(val)
}

#==================================================================================================================

beta.id <- Vectorize(function(Low, High, Cover = NA){
  
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
    
parm <- as.numeric(exp(sol$estimate))
    
q <- qbeta(p = c(p1, p2), parm[[1]], parm[[2]])
    
is.df <- function(a, b, sig = 3) round(a, sig) != round(b, sig)
    
if(is.df(L, q[1]) || is.df(U, q[2])){
      
stop("Error: \n\tUnable to find such a prior, make sure you have selected the correct values.")

  }else{
      
return(c(alpha = parm[[1]], beta = parm[[2]]))    
    }
  } 
}, c("Low", "High", "Cover"))
  
#===============================================================================================

cauchy.id <- Vectorize(function(Low, High, Cover = NA){

options(warn = -1)

coverage  <- if(is.character(Cover)) as.numeric(substr(Cover, 1, nchar(Cover)-1)) / 100 else if(is.na(Cover)) .95 else Cover
  
p1 = (1 - coverage) / 2
p2 = 1 - p1

if(p1 <= 0 || p2 >= 1 || Low > High || p1 > p2 || coverage >= 1) {

stop("\n\tUnable to find such a prior, make sure you have selected the correct values.")
  
} else {
  
f <- function(x) {   
    y <- c(Low, High) - qcauchy(c(p1, p2), location = x[1],  scale = x[2])
  }
  
parm <- optim(c(1, 1), function(x) sum(f(x)^2), control = list(reltol = (.Machine$double.eps)))[[1]]
}

q <- qcauchy(c(p1, p2), parm[[1]], parm[[2]])

is.df = function(a, b, sig = 4) round(a, sig) != round(b, sig)

if(is.df(Low, q[1]) || is.df(High, q[2])) {
  
stop("\n\tUnable to find such a prior, make sure you have selected the correct values")
  
} else { 
  
return(c(mode = parm[[1]], scale = parm[[2]])) 
  }
}, c("Low", "High", "Cover"))    
  
#===============================================================================================

normal.id <- Vectorize(function(Low, High, Cover = NA){

options(warn = -1)
  
coverage <- if(is.character(Cover)) as.numeric(substr(Cover, 1, nchar(Cover)-1)) / 100 else if(is.na(Cover)) .95 else Cover
  
p1 <- (1 - coverage) / 2 
p2 <- 1 - p1
  
q <- c(Low, High)  
alpha <- c(p1, p2)
  
is.df <- function(a, b, sig = 4) (round(a, sig) != round(b, sig))
  
if( p1 <= 0 || p2 >= 1 || q[1] >= q[2] || p1 >= p2 ) {

stop("\n\tUnable to find such a prior, make sure you have selected the correct values.")
  
} else {
    
beta <- qnorm(alpha)
    
parm <- solve(cbind(1, beta), q)
    
q <- qnorm(c(p1, p2), parm[[1]], parm[[2]])
}

if(is.df(Low, q[[1]]) || is.df(High, q[[2]])) {
  
  stop("\n\tUnable to find such a prior, make sure you have selected the correct values.")
} else {
  
  return(c(mean = parm[[1]], sd = parm[[2]]))
  
  }
}, c("Low", "High", "Cover"))
  
#===============================================================================================
  
prop.priors <- function(a, b, lo = 0, hi = 1, dist.name, yes = 55, n = 1e2, scale = .1, top = 1.5, show.prior = FALSE, bottom = 1){
  
d = dist.name
eq <- function(...){ lapply(list(...), function(x) c(x, rep(rev(x)[1], max(lengths(list(...))) - length(x)))) }
I = eq(a, b, d, lo, hi)
a = I[[1]] ; b = I[[2]] ; d = I[[3]] ; lo = I[[4]] ; hi = I[[5]]
                           
deci <- function(x, k = 3) format(round(x, k), nsmall = k)                                                                                                                           
  
  Bi = yes
  pr = show.prior
  loop = length(d)
  CI = matrix(NA, loop, 2)
  mode = numeric(loop)
  h = list()
  
  for(i in 1:loop){
    p = function(x) get(d[i])(x, a[i], b[i])*as.integer(x >= lo[i])*as.integer(x <= hi[i])
    prior = function(x) p(x)/integrate(p, lo[i], hi[i])[[1]]
    
    if(!pr){     
      likelihood = function(x) dbinom(Bi, n, x)
      k = integrate(function(x) prior(x)*likelihood(x), lo[i], hi[i])[[1]]
      posterior = function(x) prior(x)*likelihood(x) / k
      h[i] = list(curve(posterior, type = "n", ann = FALSE, yaxt = "n", xaxt = "n", add = i!= 1, bty = "n", n = 1e3))
      mode[i] = optimize(posterior, c(lo[i], hi[i]), maximum = TRUE)[[1]]
      CI[i,] = HDI(posterior)
    }
  }
  
  if(!pr){  
    
    plot(CI[, 1:2], rep(1:loop, 2), type = "n", xlim = 0:1, ylim = c(bottom*1, top*loop), ylab = NA, yaxt = "n", xaxt = "n", xlab = "Credible Interval (Proportion)", font.lab = 2, mgp = c(2, .5, 0))
    abline(h = 1:loop, col = 8, lty = 3)
    axis(1, at = axTicks(1), lab = paste0(axTicks(1)*1e2, "%"), mgp = c(2, .3, 0))
    segments(CI[, 1], 1:loop, CI[, 2], 1:loop, lend = 1, lwd = 4, col = 1:loop, xpd = NA)
    axis(2, at = 1:loop, lab = substring(d, 2), font = 2, las = 1, cex.axis = .8, tck = -.006, mgp = c(2, .3, 0))
    legend("topleft", rev(paste0(substring(d, 2), "(", round(a, 2), ", ", round(b, 2), ")")), pch = 22, pt.bg = loop:1, col = loop:1, cex = .7, pt.cex = .6, bg = 0, box.col = 0, xpd = NA, x.intersp = .5)
    box()
    for(i in 1:loop){
      polygon(x = h[[i]]$x, y = scale*h[[i]]$y +i, col = adjustcolor(i, .6), border = NA, xpd = NA)
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

prop.hyper <- function(a, b, lo = 0, hi = 1, dist.name, yes = 55, n = 1e2, show.prior = FALSE, pos = 3, top = 1.01){
  
d = dist.name
Bi = yes
pr = show.prior
  
eq <- function(...) { lapply(list(...), function(x) c(x, rep(rev(x)[1], max(lengths(list(...))) - length(x)))) }
I = eq(a, b, d, lo, hi)
a = I[[1]] ; b = I[[2]] ; d = I[[3]] ; lo = I[[4]] ; hi = I[[5]]
  
  deci <- function(x, k = 3) format(round(x, k), nsmall = k)
  
  loop = length(a)
  CI = matrix(NA, loop, 2)
  mode = numeric(loop)
  
  for(i in 1:loop){
    p = function(x) get(d[i])(x, a[i], b[i])*as.integer(x >= lo[i])*as.integer(x <= hi[i])
    prior = function(x) p(x)/integrate(p, lo[i], hi[i])[[1]]
    if(!pr){
      likelihood = function(x) dbinom(Bi, n, x)
      k = integrate(function(x) prior(x)*likelihood(x), lo[i], hi[i])[[1]]
      posterior = function(x) prior(x)*likelihood(x) / k
      mode[i] = optimize(posterior, c(lo[i], hi[i]), maximum = TRUE)[[1]]
      CI[i,] = HDI(posterior)
    }
  }
  
  if(!pr){
    graphics.off()
    original.par = par(no.readonly = TRUE)
    on.exit(par(original.par))
    
    par(mgp = c(2.2, .3, 0), mar = c(5.1, 4.1, 4.1, 3))   
    plot(CI[, 1:2], rep(1:loop, 2), type = "n", xlim = c(0, 1), ylim = c(1, top*loop), ylab = NA, yaxt = "n", xaxt = "n", xlab = "Credible Interval (Proportion)", font.lab = 2)
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
  
eq <- function(...){ lapply(list(...), function(x) c(x, rep(rev(x)[1], max(lengths(list(...))) - length(x)))) }
I = eq(a, b, d, lo, hi)
a = I[[1]] ; b = I[[2]] ; d = I[[3]] ; lo = I[[4]] ; hi = I[[5]]
  
  if(is.v(a) & pr || is.v(b) & pr) message("\tNote: You can see only '1 prior' at a time.")
  if(add & pr) message("\tNote: 'add' only works for overlying 'Credible Intervals' to compare them.")
  
  deci <- function(x, k = 3) format(round(x, k), nsmall = k)
  
  loop = length(d)
  CI = matrix(NA, loop, 2)
  mode = numeric(loop)
  
  for(i in 1:loop){
    p = function(x) get(d[i])(x, a[i], b[i])*as.integer(x >= lo[i])*as.integer(x <= hi[i])
    prior = function(x) p(x)/integrate(p, lo[i], hi[i])[[1]]
    
    if(!pr){    
      likelihood = function(x) dbinom(Bi, n, x)
      k = integrate(function(x) prior(x)*likelihood(x), lo[i], hi[i])[[1]]
      posterior = function(x) prior(x)*likelihood(x) / k
      mode[i] = optimize(posterior, c(lo[i], hi[i]), maximum = TRUE)[[1]]
      CI[i,] = HDI(posterior)
    }
  }
  
  if(!add & !pr){
    plot(rep(1:loop, 2), CI[, 1:2], type = "n", ylim = 0:1, xlim = c(1, loop), xlab = "Prior Parameter 'A'", xaxt = "n", yaxt = "n", ylab = "Credible Interval (Proportion)", font.lab = 2, mgp = c(2.3, .3, 0), cex.lab = 1.2)
    abline(v = 1:loop, col = 8, lty = 3)
    axis(2, at = axTicks(2), lab = paste0(axTicks(2)*1e2, "%"), mgp = c(2, .4, 0), las = 1)
    axis(3, at = 1:length(b), lab = deci(b), font = 2, las = 1, cex.axis = .8, mgp = c(2, .2, 0))
    text(mean(par('usr')[1:2]), 1.06*par('usr')[4], "Prior Parameter 'B'", pos = 3, cex = 1.2, xpd = NA, font = 2)
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

d.priors <- function(t, n1, n2 = NA, m, s, lo = -Inf, hi = Inf, dist.name, scale = 1, margin = 7, top = .8, show.prior = FALSE, LL = -5, UL = 5, bottom = 1){
  
  d = dist.name 
  pr = show.prior
  eq <- function(...){ lapply(list(...), function(x) c(x, rep(rev(x)[1], max(lengths(list(...))) - length(x)))) }
  I = eq(m, s, d, lo, hi)
  m = I[[1]] 
  s = I[[2]] 
  d = I[[3]] 
  lo = I[[4]] 
  hi = I[[5]]
  
  deci <- function(x, k = 3) format(round(x, k), nsmall = k)                           
  loop = length(d) 
  CI = matrix(NA, loop, 2)
  mode = numeric(loop)
  mean = numeric(loop)
  sd = numeric(loop)
  from = numeric(loop)
  to = numeric(loop) 
  h = list()
  
  N = ifelse(is.na(n2), n1, n1 * n2 / (n1 + n2))
  df = ifelse(is.na(n2), n1 - 1, n1 + n2 - 2)   
  
  options(warn = -1)
  
  for(i in 1:loop){
    p = function(x) get(d[i])(x, m[i], s[i])*as.integer(x >= lo[i])*as.integer(x <= hi[i])
    prior = function(x) p(x)/integrate(p, lo[i], hi[i])[[1]]
    
    if(!pr){     
      likelihood = function(x) dt(t, df, x*sqrt(N))
      k = integrate(function(x) prior(x)*likelihood(x), lo[i], hi[i])[[1]]
      posterior = function(x) prior(x)*likelihood(x) / k
      mean[i] = integrate(function(x) x*posterior(x), lo[i], hi[i])[[1]]
      sd[i] = sqrt(integrate(function(x) x^2*posterior(x), lo[i], hi[i])[[1]] - mean^2)
      from[i] = mean - margin * sd
      to[i] = mean + margin * sd
      mode[i] = optimize(posterior, c(from, to), maximum = TRUE, tol = 1e-10)[[1]]
      CI[i,] = HDI(posterior, LL, UL)
      h[i] = list(curve(posterior, from, to, type = "n", ann = FALSE, yaxt = "n", xaxt = "n", add = i!= 1, bty = "n", n = 5e2))
    }
  }
  
  if(!pr){
    f = sapply(h, function(x) max(x[[2]])) + 1:loop
    plot(CI[, 1:2], rep(1:loop, 2), type = "n", xlim = c(min(from), max(to)), ylim = c(bottom*1, top*max(f)), ylab = NA, yaxt = "n", xlab = bquote(bold("Credible Interval "(delta))), font.lab = 2, mgp = c(2, .5, 0))
    abline(h = 1:loop, col = 8, lty = 3)
    legend("topleft", rev(paste0(substring(d, 2), "(", round(m, 2), ", ", round(s, 2), ")")), pch = 22, pt.bg = loop:1, col = loop:1, cex = .7, pt.cex = .6, bg = 0, box.col = 0, xpd = NA, x.intersp = .5)
    box()
    segments(CI[, 1], 1:loop, CI[, 2], 1:loop, lend = 1, lwd = 4, col = 1:loop)
    axis(2, at = 1:loop, lab = d, font = 2, las = 1, cex.axis = .8, tick = FALSE, mgp = c(2, .5, 0)) ; axis(3, mgp = c(2, .5, 0))
    
    for(i in 1:loop){
      polygon(x = h[[i]]$x, y = scale*h[[i]]$y +i, col = adjustcolor(i, .6), border = NA, xpd = NA)
    }
    a = scale*(f-1:loop)+1:loop
    segments(mode, 1:loop, mode, a, lty = 3, xpd = NA, lend = 1)
    points(mode, 1:loop, pch = 21, bg = "cyan", cex = 1.1, col = 4, xpd = NA)
    I = deci(CI) ; o = deci(mode)
    text(c(CI[,1], o, CI[,2]), 1:loop, c(I[,1], o, I[,2]), pos = 3, font = 2, cex = .8, xpd = NA)
  }else{
    curve(prior, -6, 6, yaxt = "n", ylab = NA, xlab = bquote(bold("Effect Size "(delta))), bty = "n", font.lab = 2, lwd = 2, n = 1e3, main = bquote(delta*" ~ "*.(if(lo > -Inf || hi < Inf) "truncated-")*.(substring(d, 2))(.(round(m, 2)), .(round(s, 2)))), mgp = c(2, .5, 0))
  }
}

#========================================================================================================================

d.hyper <- function(t, n1, n2 = NA, m, s, lo = -Inf, hi = Inf, dist.name, LL = -4, UL = 4, pos = 3, show.prior = FALSE, top = 1.01){

  d = dist.name 
 pr = show.prior
  eq <- function(...){ lapply(list(...), function(x) c(x, rep(rev(x)[1], max(lengths(list(...))) - length(x)))) }
  I = eq(m, s, d, lo, hi)
  m = I[[1]] ; s = I[[2]] ; d = I[[3]] ; lo = I[[4]] ; hi = I[[5]] 
  
  deci <- function(x, k = 3) format(round(x, k), nsmall = k)
  
  N = ifelse(is.na(n2), n1, n1 * n2 / (n1 + n2))
  df = ifelse(is.na(n2), n1 - 1, n1 + n2 - 2)   
  
  options(warn = -1)
  original.par = par(no.readonly = TRUE)
  on.exit(par(original.par))
  
  loop = length(m)
  CI = matrix(NA, loop, 2)
  mode = numeric(loop)
  mean = numeric(loop)
  sd = numeric(loop)
  from = numeric(loop)
  to = numeric(loop)
  
  for(i in 1:loop){
    p = function(x) get(d[i])(x, m[i], s[i])*as.integer(x >= lo[i])*as.integer(x <= hi[i])
    prior = function(x) p(x)/integrate(p, lo[i], hi[i])[[1]]
    
if(!pr){   
    likelihood = function(x) dt(t, df, x*sqrt(N))
    k = integrate(function(x) prior(x)*likelihood(x), lo[i], hi[i])[[1]]
    posterior = function(x) prior(x)*likelihood(x) / k
    mode[i] = optimize(posterior, c(LL, UL), maximum = TRUE, tol = 1e-10)[[1]]
    mean[i] = integrate(function(x) x*posterior(x), lo[i], hi[i])[[1]]
    sd[i] = sqrt(integrate(function(x) x^2*posterior(x), lo[i], hi[i])[[1]] - mean^2)
    CI[i,] = HDI(posterior, LL, UL)
    from[i] = mean - 7 * sd
    to[i] = mean + 7 * sd 
     }
  }
  
if(!pr){   
  par(mgp = c(2, .5, 0), mar = c(5.1, 4.1, 4.1, 3))   
  plot(CI[, 1:2], rep(1:loop, 2), type = "n", xlim = c(min(from), max(to)), ylim = c(1, top*loop), ylab = NA, yaxt = "n", xlab = bquote(bold("Credible Interval "(delta))), font.lab = 2)
  abline(h = 1:loop, col = 8, lty = 3)
  segments(CI[, 1], 1:loop, CI[, 2], 1:loop, lend = 1, col = "red4")  
  points(mode, 1:loop, pch = 21, bg = "red4", cex = .8, col = "red4", xpd = NA)
  axis(2, at = 1:length(m), lab = deci(m), font = 2, las = 1, cex.axis = .8, tick = FALSE, mgp = c(2, .2, 0))
  axis(4, at = 1:length(s), lab = deci(s), font = 2, las = 1, cex.axis = .8, tick = FALSE, mgp = c(2, .2, 0))
  text(par('usr')[1:2], par('usr')[4], c("M", "S"), pos = 3, cex = 1.5, xpd = NA, font = 2)
  I = deci(CI) ; o = deci(mode)
  text(mode, 1:loop, paste0("[", I[,1], ",  ", o, ",  ", I[,2], "]"), pos = pos, cex = .8, xpd = NA)
}else{
  
curve(prior, -6, 6, yaxt = "n", ylab = NA, xlab = bquote(bold("Effect Size "(delta))), bty = "n", font.lab = 2, lwd = 2, n = 1e3, main = bquote(delta*" ~ "*.(if(lo > -Inf || hi < Inf) "truncated-")*.(substring(d, 2))(.(round(m, 2)), .(round(s, 2)))), mgp = c(2, .5, 0))
  }  
}

#===================================================================================================================

ms.d.hyper <- function(t, n1, n2 = NA, m, s, lo = -Inf, hi = Inf, dist.name, add = FALSE, 
                      col = 1, top = 6, margin = 1.01, LL = -5, UL = 5, show.prior = FALSE){
  
  d = dist.name 
  pr = show.prior
  eq <- function(...){ lapply(list(...), function(x) c(x, rep(rev(x)[1], max(lengths(list(...))) - length(x)))) }
  I = eq(m, s, d, lo, hi)
  m = I[[1]] ; s = I[[2]] ; d = I[[3]] ; lo = I[[4]] ; hi = I[[5]]
    
  deci <- function(x, k = 3) format(round(x, k), nsmall = k)

  if(add & pr) message("\tNote: 'add' only works for overlying 'Credible Intervals' to compare them.")

  N = ifelse(is.na(n2), n1, n1 * n2 / (n1 + n2))
  df = ifelse(is.na(n2), n1 - 1, n1 + n2 - 2)   
  
  options(warn = -1)
  loop = length(m) 
  CI = matrix(NA, loop, 2)
  mode = numeric(loop)
  mean = numeric(loop)
  sd = numeric(loop)
  from = numeric(loop)
  to = numeric(loop)
  
  for(i in 1:loop){
    p = function(x) get(d[i])(x, m[i], s[i])*as.integer(x >= lo[i])*as.integer(x <= hi[i])
    prior = function(x) p(x)/integrate(p, lo[i], hi[i])[[1]]
    
if(!pr){    
    likelihood = function(x) dt(t, df, x*sqrt(N))
    k = integrate(function(x) prior(x)*likelihood(x), lo[i], hi[i])[[1]]
    posterior = function(x) prior(x)*likelihood(x) / k
    mode[i] = optimize(posterior, c(LL, UL), maximum = TRUE, tol = 1e-10)[[1]]
    mean[i] = integrate(function(x) x*posterior(x), lo[i], hi[i])[[1]]
    sd[i] = sqrt(integrate(function(x) x^2*posterior(x), lo[i], hi[i])[[1]] - mean^2)
    CI[i,] = HDI(posterior, LL, UL)
    from[i] = mean - top * sd
    to[i] = mean + top * sd
}
  }
  
  if(!add & !pr){
    plot(rep(1:loop, 2), CI[, 1:2], type = "n", ylim = c(min(from), max(to)), xlim = c(1, margin*loop), xlab = "Prior Parameter 'M'", xaxt = "n", ylab = bquote(bold("Credible Interval "(delta))), font.lab = 2)
    abline(v = 1:loop, col = 8, lty = 3, mgp = c(2, .5, 0))
    axis(3, at = 1:length(s), lab = round(s, 3), font = 2, las = 1, cex.axis = .8, mgp = c(2, .2, 0))
    text(mean(par('usr')[1:2]), 1.06*par('usr')[4], "Prior Parameter 'S'", pos = 3, cex = 1, xpd = NA, font = 2)
    axis(1, at = 1:length(m), lab = round(m, 3), font = 2, las = 1, cex.axis = .8, mgp = c(2, .5, 0))
  }
  
  if(!pr){
  segments(1:loop, CI[, 1], 1:loop, CI[, 2], lend = 1, col = col)  
  lines(1:loop, mode, col = col, lty = 3)
  points(1:loop, mode, pch = 21, bg = col, cex = .8, col = col, xpd = NA)
  }
  
  if(!add & pr){
    curve(prior, -6, 6, yaxt = "n", ylab = NA, xlab = bquote(bold("Effect Size "(delta))), bty = "n", font.lab = 2, lwd = 2, n = 1e3, main = bquote(delta*" ~ "*.(if(lo > -Inf || hi < Inf) "truncated-")*.(substring(d, 2))(.(round(m, 2)), .(round(s, 2)))), mgp = c(2, .5, 0))
  }
}

#=================================================================================================================

prop.update <- function(n = 100, yes = 55, top = 5, scale = .1, a = 1.5, b = 1.5, dist.name = "dbeta", prior.scale = 1, level = .95, show.prior = FALSE, lo = 0, hi = 1, tol = 1e5){

pri <- show.prior
s <- yes  
d <- dist.name
if(tol < 1e4) stop("'tol' must be '10,000' or larger.")
is.v = function(x) length(x) > 1
if(is.v(d) || is.v(a) || is.v(b)) stop("Error: Choose only 'one' prior knowledge base at a time.")

eq <- function(...){ lapply(list(...), function(x) c(x, rep(rev(x)[1], max(lengths(list(...))) - length(x)))) }
deci <- function(x, k = 3) format(round(x, k), nsmall = k) 
I <- eq(n, s) ; n <- I[[1]] ; s <- I[[2]]
loop <- length(n) 

  props <- seq(lo, hi, 1/tol)
  pr <- get(d)(props, a, b)
  pr <- tol * pr / sum(pr)
  
  original.par = par(no.readonly = TRUE)
  on.exit(par(original.par))
   
  par(mar = c(5, 6.8, 4, 2))
  plot(pr~props, ylim = c(0, top*loop), type = "n", yaxs = "i", ylab = NA, xlab = "Proportion", font.lab = 2, axes = FALSE, mgp = c(2, .4, 0),main = if(pri) bquote(Proportion*" ~ "*.(if(lo > 0 || hi < 1) "truncated-")*.(substring(d, 2))(.(round(a, 2)), .(round(b, 2)))) else NA)
  axis(1, at = axTicks(1), lab = paste0(axTicks(1)*1e2, "%"), mgp = c(2, .3, 0))

if(!pri){  
  abline(h = 1:loop+1, col = 8, lty = 3)
  axis(2, at = 0:loop+1, lab = c("Base knowledge", paste0("Study ", 1:loop)), las = 1, font = 2, cex.axis = .9, mgp = c(2, .2, 0), tick = FALSE, xpd = NA)
}  
  polygon(x = c(lo, props, hi), y = prior.scale*c(0, pr, 0), col = adjustcolor(8, .8))
  
  I = hdi(x = props, y = pr, level = level)
  
  m = props[which.max(pr)]
  y = prior.scale*(pr[which.max(pr)])
  segments(I[1], 0, I[2], 0, lend = 1, lwd = 4, xpd = NA)
  points(m, 0, pch = 19, xpd = NA, cex = 1.4)
  segments(m, 0, m, y, lty = 3)
  text(c(.85*I[1], m, I[2]), 0, paste0(round(c(I[1], m, I[2])*1e2, 4), "%"), pos = 3, cex = .8, font = 2, xpd = NA)

if(!pri){
  for(i in 1:loop) {
    ps <- dbinom(s[i], n[i], props) * pr
    ps <- tol * ps / sum(ps)
polygon(y = scale*ps+i+1, x = props, col = adjustcolor(i+1, .5), border = NA, xpd = NA)
I = hdi(x = props, y = ps, level = level)
m = props[which.max(ps)]
q = deci(I*1e2 , 2); o = deci(m*1e2, 2)
y = ps[which.max(ps)]*scale + (i+1)
segments(I[1], i+1, I[2], i+1, lend = 1, lwd = 3, col = i +1)
segments(m, i+1, m, y, lty = 3, xpd = NA)
text(m, i+1, paste0(q[1], "%", "     ", o, "%", "     ", q[2], "%"), pos = 3, cex = .7, font = 2)
points(m, i+1, pch = 21, bg = "cyan", col = "magenta")

    pr <- ps
    }
  }
}

#=======================================================================================================================

d.update <- function(t, n1, n2 = NA, top = 5, scale = .1, m = 0, s = 1, dist.name, prior.scale = 1, level = .95, show.prior = FALSE, lo = -2, hi = 2, tol = 1e4){
  
  pri <- show.prior
  d <- dist.name
  if(tol < 1e4) stop("'tol' must be '10,000' or larger.")
  is.v = function(x) length(x) > 1
  if(is.v(d)|| is.v(m) || is.v(s)) stop("Error: Choose only 'one' prior knowledge base at a time.")
  
  eq <- function(...){ lapply(list(...), function(x) c(x, rep(rev(x)[1], max(lengths(list(...))) - length(x)))) }
  deci <- function(x, k = 3) format(round(x, k), nsmall = k) 
  I <- eq(t, n1, n2) 
  t <- I[[1]]  
  n1 <- I[[2]]  
  n2 <- I[[3]] 
  loop <- length(t) 
  
  ds <- seq(lo, hi, 1/tol)
  pr <- get(d)(ds, m, s)
  pr <- tol * pr / sum(pr)
  
  original.par = par(no.readonly = TRUE)
  on.exit(par(original.par))
  
  par(mar = c(5, 6.8, 4, 2))
  plot(pr~ds, ylim = c(0, top*loop), xlim = c(lo, hi), type = "n", xaxs = "i", yaxs = "i", ylab = NA, xlab = bquote(bold("Effect Size"~ (delta))), font.lab = 2, mgp = c(2, .4, 0),main = if(pri) bquote("Effect Size"*" ~ "*.(if(lo > -Inf || hi < Inf) "truncated-")*.(substring(d, 2))(.(round(m, 2)), .(round(s, 2)))) else NA, yaxt = "n", bty = "n")
  
  if(!pri){  
    abline(h = 1:loop+1, col = 8, lty = 3)
    axis(2, at = 0:loop+1, lab = c("Base knowledge", paste0("Study ", 1:loop)), las = 1, font = 2, cex.axis = .9, mgp = c(2, .2, 0), tick = FALSE, xpd = NA)
  }  

  polygon(x = c(lo, ds, hi), y = prior.scale*c(0, pr, 0), col = adjustcolor(8, .8))
  
  I = hdi(x = ds, y = pr, level = level)
  
  mode = ds[which.max(pr)]
  y = prior.scale*(pr[which.max(pr)])
  segments(I[1], 0, I[2], 0, lend = 1, lwd = 4, xpd = NA)
  points(mode, 0, pch = 19, xpd = NA, cex = 1.4)
  segments(mode, 0, mode, y, lty = 3)
  text(c(.85*I[1], mode, I[2]), 0, paste0(round(c(I[1], mode, I[2]), 3)), pos = 3, cex = .8, font = 2, xpd = NA)
  
  N = ifelse(is.na(n2), n1, n1 * n2 / (n1 + n2))
  df = ifelse(is.na(n2), n1 - 1, n1 + n2 - 2) 
  
  options(warn = -1)
  if(!pri){
    for(i in 1:loop) {
     
      ps <- dt(t[i], df[i], ds*sqrt(N[i])) * pr
      ps <- tol * ps / sum(ps)
      polygon(y = scale*ps+i+1, x = ds, col = adjustcolor(i+1, .5), border = NA, xpd = NA)
      I = hdi(x = ds, y = ps, level = level)
      mode = ds[which.max(ps)]
      q = deci(I, 3); o = deci(mode, 3)
      y = ps[which.max(ps)]*scale + (i+1)
      segments(I[1], i+1, I[2], i+1, lend = 1, lwd = 3, col = i +1)
      segments(mode, i+1, mode, y, lty = 3, xpd = NA)
      text(mode, i+1, paste0(q[1], "     ", o, "     ", q[2]), pos = 3, cex = .7, font = 2)
      points(mode, i+1, pch = 21, bg = "cyan", col = "magenta")
      
      pr <- ps
    }
  }
}

