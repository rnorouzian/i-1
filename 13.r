d.hyper <- function(t, n1, n2 = NA, m, s, lo = -Inf, hi = Inf, dist.name = "dcauchy", LL = -3, UL = 3){
  
  is.v = function(x) length(x) > 1
  if(is.v(m) & is.v(s)) stop("Error: Explore 'm' and 's' one at a time.")  
  if(is.v(dist.name)) stop("Error: Explore only one 'dis.name' at a time.")
  
  deci = function(x, k = 3) format(round(x, k), nsmall = k)
  
  d = dist.name
  f <- if(is.v(m)) function(x) get(d)(x, m[i], s) else function(x) get(d)(x, m, s[i])
  
  N = ifelse(is.na(n2), n1, n1 * n2 / (n1 + n2))
  df = ifelse(is.na(n2), n1 - 1, n1 + n2 - 2)   
  
  options(warn = -1)
  original.par = par(no.readonly = TRUE)
  on.exit(par(original.par))
  
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
    likelihood = function(x) dt(t, df, x*sqrt(N))
    k = integrate(function(x) prior(x)*likelihood(x), lo, hi)[[1]]
    posterior = function(x) prior(x)*likelihood(x) / k
    mode[i] = optimize(posterior, c(LL, UL), maximum = TRUE, tol = 1e-20)[[1]]
    mean[i] = integrate(function(x) x*posterior(x), lo, hi)[[1]]
    sd[i] = sqrt(integrate(function(x) x^2*posterior(x), lo, hi)[[1]] - mean^2)
    CI[i,] = HDI(posterior, LL, UL)
    from[i] = mean - 7 * sd
    to[i] = mean + 10 * sd 
  }
  par(mgp = c(2, .5, 0))   
  plot(CI[, 1:2], rep(1:loop, 2), type = "n", xlim = c(min(from), max(to)), ylim = c(1, 1.01*loop), ylab = if(is.v(m)) "Prior Mean" else if(is.v(s)) "Prior SD" else NA, yaxt = "n", xlab = "Credible Interval 'd'", font.lab = 2)
  abline(h = 1:loop, col = 8, lty = 3)
  segments(CI[, 1], 1:loop, CI[, 2], 1:loop, lend = 1, col = "red4")  
  legend("topleft", if(is.v(m)) "The effect of\nPrior Mean" else if(is.v(s)) "The effect of\nPrior Width (SD)" else " ", bty = "n", text.font = 2, cex = .7)  
  points(mode, 1:loop, pch = 21, bg = "red4", cex = .8, col = "red4", xpd = NA)
  lab = if(is.v(m)) m else if(is.v(s)) s
  if(is.v(m) || is.v(s)) axis(2, at = 1:if(is.v(m))length(m) else if(is.v(s))length(s) else NA, lab = lab, font = 2, las = 1, cex.axis = .8, tick = FALSE, mgp = c(1, .1, 0))
  
  I = deci(CI) ; o = deci(mode)
  
  text(.7*par('usr')[2], 1:loop, paste0("[", I[,1], ",  ", o, ",  ", I[,2], "]"), cex = .8, xpd = NA)
}