source("https://raw.githubusercontent.com/izeh/i/master/i.r")

prop.priors <- function(a, b, lo = 0, hi = 1, dist.name, yes = 55, n = 1e2, scale = .1, top = 1.5, show.prior = FALSE){
  
  is.v = function(x) length(x) > 1
  if(is.v(a) || is.v(b)) stop("Error: Only 'dist.name' can be a vector with length > 1.")
  deci = function(x, k = 3) format(round(x, k), nsmall = k)
   d = dist.name ; Bi = yes ; pr = show.prior
loop = length(d)
  CI = matrix(NA, loop, 2)
mode = numeric(loop)
   h = list()

for(i in 1:loop){
         p = function(x) get(d[i])(x, a, b)*as.integer(x >= lo)*as.integer(x <= hi)
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
legend("topleft", paste0(substring(d, 2), "(", round(a, 2), ", ", round(b, 2), ")"), pch = 22, pt.bg = 1:loop, col = 1:loop, cex = .7, bg = NA, bty = "n", pt.cex = .6)
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
