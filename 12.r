source("https://raw.githubusercontent.com/izeh/i/master/i.r")

ms.d.hyper <- function(t, n1, n2 = NA, m, s, lo = -Inf, hi = Inf, dist.name, 
                       margin = 1.01, add = FALSE, col = 1, top = 6, LL = -8, UL = 8){
  
is.v = function(x) length(x) > 1
d = dist.name   
  
if(is.v(d)) stop("Error: Only 'm' or 's' can be a vector with length > 1.")  
if(is.v(m) & is.v(s)) stop("Error: Explore 'm' and 's' one at a time.")
f <- if(is.v(m)) function(x) get(d)(x, m[i], s) else function(x) get(d)(x, m, s[i])
  
 N = ifelse(is.na(n2), n1, n1 * n2 / (n1 + n2))
df = ifelse(is.na(n2), n1 - 1, n1 + n2 - 2)   
  
options(warn = -1)
  
loop = if(is.v(m)) length(m) else length(s) 
  CI = matrix(NA, loop, 2)
mode = numeric(loop)
mean = numeric(loop)
  sd = numeric(loop)
from = numeric(loop)
  to = numeric(loop)
  
for(i in 1:loop){
        p = function(x) f(x)*as.integer(x >= lo)*as.integer(x <= hi)
    prior = function(x) p(x)/integrate(p, lo, hi)[[1]]
     like = function(x) dt(t, df, x*sqrt(N))
        k = integrate(function(x) prior(x)*like(x), lo, hi)[[1]]
posterior = function(x) prior(x)*like(x) / k
  mode[i] = optimize(posterior, c(LL, UL), maximum = TRUE, tol = 1e-20)[[1]]
  mean[i] = integrate(function(x) x*posterior(x), lo, hi)[[1]]
    sd[i] = sqrt(integrate(function(x) x^2*posterior(x), lo, hi)[[1]] - mean^2)
   CI[i,] = HDI(posterior, LL, UL)
  from[i] = mean - top * sd
    to[i] = mean + top * sd 
}
  
if(!add){
  plot(rep(1:loop, 2), CI[, 1:2], type = "n", ylim = c(min(from), max(to)), xlim = c(1, margin*loop), xlab = if(is.v(m)) "Prior Mean" else if(is.v(s)) "Prior SD" else NA, xaxt = "n", ylab = bquote(bold("Credible Interval "~ (delta))), font.lab = 2, mgp = c(2, .5, 0))
  abline(v = 1:loop, col = 8, lty = 3)
  legend("topleft", if(is.v(m)) "The effect of\nPrior Mean" else if(is.v(s)) "The effect of\nPrior Width (SD)" else " ", bty = "n", text.font = 2, cex = .8)  
  at = if(is.v(m))length(m) else if(is.v(s))length(s)
  if(is.v(m) || is.v(s)) axis(1, at = 1:at, lab = if(is.v(m)) m else if(is.v(s)) s, font = 2, las = 1, cex.axis = .8, mgp = c(2, .5, 0))
}
  segments(1:loop, CI[, 1], 1:loop, CI[, 2], lend = 1, col = col)  
  lines(1:loop, mode, col = col, lty = 3)
  points(1:loop, mode, pch = 21, bg = col, cex = .8, col = col, xpd = NA)
}
