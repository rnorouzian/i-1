biling.survey = function(N = 40, pBi = .5, ...){
  
  Responses = sample(0:1, size = N, prob = c(1 - pBi, pBi), replace = TRUE)  # Generate N responses from parents
                                                                                 # (B = 1, M = 0)
  Prop = cumsum(Responses) / 1:N   # Compute running proportion of B as each parent responds
   
  par(las = 1, tck = -.02, font.lab = 2, ...)  
  plot.ts(Prop, ty = "o", ylim = 0:1, yaxt = "n", pch = 21, bg = 3, xlab = "Number of Parents", 
          ylab = "Proportion of (B)")
  
  axis(2, at = seq(0, 1L, len = 6), lab = paste0(seq(0, 1e2, len = 6), "%"))
  
  ResponseSeq = paste(c("M", "B")[Responses[1L:1e1] + 1L], collapse = "")
  
  Display = paste0("Response Sequence = ", ResponseSeq, ". . .")
  
  text(N, c(1, .95), c(Display, paste0("Proportion of B in ", N, " responses = ", Prop[N]*1e2, "%")), adj = c(1, .5), col = "red4", font = 2, cex = .8)
}
#Example of use:
biling.survey(N = 100, pBi = .75)
