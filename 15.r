source("https://raw.githubusercontent.com/izeh/i/master/i.r")

prop.hyper <- function(a, b, lo = 0, hi = 1, dist.name, yes = 55, n = 1e2, show.prior = FALSE, pos = 3){
  
  d = dist.name
  Bi = yes
  pr = show.prior
  is.v = function(x) length(x) > 1
  if(is.v(d)) stop("Error: Only 'a' or 'b' can be a vector with length > 1.")
  if(is.v(a) & is.v(b)) stop("Error: Explore 'a' and 'b' one at a time.")
  
  deci = function(x, k = 3) format(round(x, k), nsmall = k)
  
  f <- if(is.v(a)) function(x) get(d)(x, a[i], b) else function(x) get(d)(x, a, b[i])
  
  loop = if(is.v(a)) length(a) else length(b) 
  CI = matrix(NA, loop, 2)
  mode = numeric(loop)
  
  for(i in 1:loop){
    p = function(x) f(x)*as.integer(x >= lo)*as.integer(x <= hi)
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
    
    par(mgp = c(2, .3, 0))   
    plot(CI[, 1:2], rep(1:loop, 2), ty = "n", xlim = c(0, 1), ylim = c(1, 1.01*loop), ylab = if(is.v(a)) "Prior 'A' parameter" else if(is.v(b)) "Prior 'B' parameter" else NA, yaxt = "n", xaxt = "n", xlab = "Credible Interval 'Proportion'", font.lab = 2)
    abline(h = 1:loop, col = 8, lty = 3)
    axis(1, at =  axTicks(1), lab = paste0(axTicks(1)*1e2, "%"))
    segments(CI[, 1], 1:loop, CI[, 2], 1:loop, lend = 1, col = "red4", xpd = NA)
    points(mode, 1:loop, pch = 21, bg = "red4", cex = .8, col = "red4", xpd = NA)
    lab = if(is.v(a)) round(a, 3) else if(is.v(b)) round(b, 3)
    if(is.v(a) || is.v(b)) axis(2, at = 1:if(is.v(a))length(a) else if(is.v(b))length(b), lab = lab, font = 2, las = 1, cex.axis = .8, tck = -.006)
    I = deci(CI*1e2 , 2); o = deci(mode*1e2, 2)
    text(mode, 1:loop, paste0("[", I[,1], "%", ", ", o, "%", ", ", I[,2], "%", "]"), cex = .7, pos = pos, xpd = NA)
  }else{
    curve(prior, lo, hi, yaxt = "n", ylab = NA, xaxt = "n", xlab = "Proportion", bty = "n", font.lab = 2, lwd = 2, n = 5e2, main = bquote(Proportion*" ~ "*.(substring(d, 2))(.(round(a, 2)), .(round(b, 2)))))
    axis(1, at = axTicks(1), lab = paste0(axTicks(1)*1e2, "%"), mgp = c(2, .4, 0))
  }
}