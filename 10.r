source("https://raw.githubusercontent.com/izeh/i/master/i.r")

d.priors <- function(n1, n2 = NA, t, m, s, lo = -Inf, hi = Inf, dist.name = c("dnorm", "dcauchy", "dlogis"),
                     scale = .35, margin = 7, top = .71, LL = -9, UL = 9){
  
  is.v = function(x) length(x) > 1
  if(is.v(m) || is.v(s)) stop("Error: Only 'dist.name' can be a vector with length > 1.")
  
   d = dist.name
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
original.par = par(no.readonly = TRUE)
on.exit(par(original.par))
   
for(i in 1:loop){
         p = function(x) get(d[i])(x, m, s)*as.integer(x >= lo)*as.integer(x <= hi)
     prior = function(x) p(x)/integrate(p, lo, hi)[[1]]
likelihood = function(x) dt(t, df, x*sqrt(N))
         k = integrate(function(x) prior(x)*likelihood(x), lo, hi)[[1]]
 posterior = function(x) prior(x)*likelihood(x) / k
   mode[i] = optimize(posterior, c(LL, UL), maximum = TRUE, tol = 1e-20)[[1]]
   mean[i] = integrate(function(x) x*posterior(x), lo, hi)[[1]]
     sd[i] = sqrt(integrate(function(x) x^2*posterior(x), lo, hi)[[1]] - mean^2)
    CI[i,] = HDI(posterior, LL, UL)
   from[i] = mean - margin * sd
     to[i] = mean + margin * sd  
      h[i] = list(curve(posterior, from, to, type = "n", ann = FALSE, yaxt = "n", xaxt = "n", add = i!= 1, bty = "n", n = 5e2))
}

f = sapply(h, function(x) max(x[[2]])) + 1:loop
par(mgp = c(2, .5, 0))                                             
plot(CI[, 1:2], rep(1:loop, 2), type = "n", xlim = c(min(from), max(to)), ylim = c(1, top*max(f)), ylab = NA, yaxt = "n", xlab = "Credible Interval 'd'", font.lab = 2)
abline(h = 1:loop, col = 8, lty = 3)
segments(CI[, 1], 1:loop, CI[, 2], 1:loop, lend = 1, lwd = 4, col = 1:loop)
axis(2, at = 1:loop, lab = d, font = 2, las = 1, cex.axis = .8, tick = FALSE) ; axis(3)

for(i in 1:loop){
  polygon(x = h[[i]]$x, y = scale*h[[i]]$y +i, col = adjustcolor(i, .4), border = NA, xpd = NA)
}
  a = scale*(f-1:loop)+1:loop
  segments(mode, 1:loop, mode, a, lty = 3, xpd = NA, lend = 1)
  points(mode, 1:loop, pch = 21, bg = "cyan", cex = 1.1, col = 4, xpd = NA)
  deci = function(x, k = 3) format(round(x, k), nsmall = k)
  I = deci(CI) ; o = deci(mode)
  text(c(CI[,1], o, CI[,2]), 1:loop, c(I[,1], o, I[,2]), pos = 3, font = 2, cex = .8, xpd = NA)

structure(list(post.mode = deci(mode, 6), post.mean = deci(mean, 6), post.sd = deci(sd, 6)), class = "power.htest")
}
