biling.survey = function(N = 40, pBi = .5, ...){
  
  ResponseSequence = sample(x = 0:1, prob = c(1 - pBi, pBi), size = N, replace = TRUE)  # Generate N responses from parents
                                                                                        # (B = 1, M = 0)
  Prop = cumsum(ResponseSequence) / 1:N   # Compute running proportion of B as each parent responds
   
  par(las = 1, tck = -.02, font.lab = 2, cex.lab = 1, ...)  
  plot.ts(Prop, ty = "o", ylim = c(0, 1), yaxt = "n", pch = 21, bg = 3, xlab = "Number of Parents", 
          ylab = "Proportion of (B)")
  
  axis(2, at = seq(0, 1L, len = 6), lab = paste0(seq(0, 1e2, len = 6), "%"))
  
  ResponseLetters = paste( c("M","B")[ResponseSequence[1L:1e1] + 1L], collapse = "")
  
  Display = paste0("Response Sequence = ", ResponseLetters, ". . .")
  
  text(N, c(1, .95), c(Display, paste0("Proportion of B in ", N, " responses = ", Prop[N]*1e2, "%")), adj = c(1, .5), col = "red4", font = 2, cex = .8)
}
#Example of use:
biling.survey(N = 100, pBi = .75)
