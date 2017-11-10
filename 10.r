p.CI <- function(s, n, conf.level = .95) { ## s = no. of agreements
                                           ## n = total no. of parents
  
  par(mfrow = c(2, 1), mar = c(3, 4.1, 3, .5), las = 1, lend = 1)
  alpha = (1 - conf.level)/2
  
  L <- if(s == 0) 0 else qbeta(alpha, s, n - s + 1)
  U <- if(s == n) 1 else qbeta(1 - alpha, s + 1, n - s)
  
  curve(dbeta(x, s, n - s + 1), n = 1e4, lwd = 2, col = 4, ylab = "Density", main = "Beta Distributions", xlab = NA)
  curve(dbeta(x, s + 1, n - s), n = 1e4, lwd = 2, col = 2, add = TRUE)
  
  abline(v = c(L, U), col = c(4, 2), lty = 3)
  
  segments(L, 0, U, 0, lwd = 4) ; points(s/n, 0, pch = 19, cex = 2)
  
  text(c(L, U, s/n), 0, signif(c(L, U, s/n), 2), pos = 3, font = 2, xpd = TRUE)
  
  plot(x <- 0:n, dbinom(x, s = n, p = L), ty = "h", lwd = 2, col = 4, ylab = "Probability", main = "Binomial Distributions")
  lines(x <- 0:n, dbinom(x, s = n, p = U), ty = "h", lwd = 2, col = 2, lty = 2, xpd = TRUE)
  
  return(c(L, U))
}
#Example of use:
p.CI(s = 55, n = 100)
