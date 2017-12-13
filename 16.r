source("https://raw.githubusercontent.com/izeh/i/master/i.r")

ab.prop.hyper <- function(a, b, lo = 0, hi = 1, dist.name, add = FALSE, yes = 55, n = 1e2, col = 1, show.prior = FALSE){
  
  d = dist.name ; Bi = yes ; pr = show.prior    
  is.v = function(x) length(x) > 1
  if(is.v(d)) stop("Error: Only 'a' or 'b' can be a vector with length > 1.")
  if(is.v(a) & is.v(b)) stop("Error: Explore 'a' and 'b' one at a time.")
  if(is.v(a) & pr || is.v(b) & pr) message("\tNote: You can see only '1 prior' at a time.")
  if(add & pr) message("\tNote: 'add' only works for overlying 'Credible Intervals' to compare them.")
  
  deci = function(x, k = 3) format(round(x, k), nsmall = k)
  
  f <- if(is.v(a)) function(x) get(d)(x, a[i], b) else function(x) get(d)(x, a, b[i])
  
  loop = if(is.v(a)) length(a) else length(b) 
  CI = matrix(NA, loop, 2)
  mode = numeric(loop)
  
  for(i in 1:loop){
    p = function(x) f(x)*as.integer(x >= lo)*as.integer(x <= hi)
    prior = function(x) p(x)/integrate(p, lo, hi)[[1]]
    
    if(!pr){    
      likelihood = function(x) dbinom(Bi, n, x)
      k = integrate(function(x) prior(x)*likelihood(x), lo, hi)[[1]]
      posterior = function(x) prior(x)*likelihood(x) / k
      mode[i] = optimize(posterior, c(lo, hi), maximum = TRUE, tol = 1e-20)[[1]]
      CI[i,] = HDI(posterior)
    }
  }
  
  if(!add & pr){curve(prior, lo, hi, yaxt = "n", xaxt = "n", ylab = NA, xlab = "Proportion", bty = "n", font.lab = 2, lwd = 2, n = 1e3, main = bquote(Proportion*" ~ "*.(substring(d, 2))(.(round(a, 2)), .(round(b, 2)))))
    axis(1, at = axTicks(1), lab = paste0(axTicks(1)*1e2, "%"), mgp = c(2, .4, 0))
  }
  
  if(!add & !pr){
    plot(rep(1:loop, 2), CI[, 1:2], type = "n", ylim = 0:1, xlim = c(1, loop), xlab = if(is.v(a)) "Prior 'A' parameter" else if(is.v(b)) "Prior 'B' parameter" else NA, xaxt = "n", yaxt = "n", ylab = "Credible Interval (Proportion)", font.lab = 2, mgp = c(2.3, .3, 0))
    abline(v = 1:loop, col = 8, lty = 3)
    axis(2, at = axTicks(2), lab = paste0(axTicks(2)*1e2, "%"), mgp = c(2, .4, 0), las = 1)
    legend("topleft", if(is.v(a)) "The effect of\nPrior 'A' parameter" else if(is.v(b)) "The effect of\nPrior 'B' parameter" else " ", bty = "n", text.font = 2, cex = .8)  
    at = if(is.v(a))length(a) else if(is.v(b))length(b)
    if(is.v(a) || is.v(b)) axis(1, at = 1:at, lab = if(is.v(a)) round(a, 3) else if(is.v(b)) round(b, 3), font = 2, las = 1, cex.axis = .8, mgp = c(2, .3, 0))
  }
  
  if(!pr){
    segments(1:loop, CI[, 1], 1:loop, CI[, 2], lend = 1, col = col)  
    lines(1:loop, mode, col = col, lty = 3)
    points(1:loop, mode, pch = 21, bg = col, cex = .8, col = col, xpd = NA)
  }
}
