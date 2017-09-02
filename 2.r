CI.bi = function(n, p, n.sim, ylabel = FALSE){
  
  fun = function(){
    x = rbinom(1, n, p)
   pe = x/n
   CI = binom.test(x, n, p)[[4]]
   c(CI[1], CI[2], pe)
  }
  
    sim = t(replicate(n.sim, fun()))
capture = sim[ ,1] <= p & p <= sim[ ,2]
  
  par(mgp = c(2, .2, 0), tck = -.015)
  plot(sim[, 1:2], rep(1:n.sim, 2), ty = "n", ylab = NA, yaxt = "n", xaxt = "n", xlab = "Proportion of (B)", font.lab = 2, bty = "n")
  
  abline(h = 1:n.sim, col = 8, lty = 3)
  abline(v = p, lty = 2, col = 2)
  segments(sim[ ,1], 1:n.sim, sim[ ,2], 1:n.sim, lend = 1, col = ifelse(capture, 1, 2))
  axis(1, at = axTicks(1), labels = paste0(axTicks(1)*1e2, "%"))
  axis(1, at = p, col.axis = 2, col = 2, font = 2)
  if(ylabel) axis(2, at = 1:n.sim, labels = paste0("Repeat ", rev(1:n.sim)), font = 2, las = 1, cex.axis = .8, tck = -.006)
  points(sim[, 3], 1:n.sim, pch = 19, col = ifelse(capture, 1, 2), cex = ifelse(n.sim > 50, .6, .65))
  
  cat("Coverage =", mean(capture)*1e2, "%", "\nLong-run proportion =", mean(sim[, 3])*1e2, "%")
}
# Example of use:
CI.bi(n = 100, p = .75, n.sim = 20, ylabel = TRUE)