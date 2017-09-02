

ttest = function(n1 = 5, n2 = 5, 
                 min.score = 0, 
                 max.score = 20, 
                 simulation = FALSE, 
                 n.sim = 1, 
                 sim.time, 
                 paired = FALSE,
                 subjects = TRUE,
                 descriptives = TRUE,
                 correlation = .5,
                 researcher.tool = FALSE,
                 effect.size = .5,
                 ...
){
  
  ttest.sim = function(){
    
    options(warn = -1)
    
    if(min.score >= max.score){
      stop("\n\tYour \"min.score\" must be smaller than your \"max.score\".")  }
    
    if(paired & n1 != n2){ 
      warning(message("\n\tIn a paired design, Your \"n1\" must be equal to your \"n2\".\n")) ; n2 = n1 }
    
    if(!researcher.tool){
      
      beta = qnorm(c(1e-10, .9999999999))
         q = c(min.score, max.score)
      
      
  mu.sigma = solve(cbind(1L, beta), q)
      
      mean = mu.sigma[[1]]
        sd = mu.sigma[[2]]
      
     coeff = effect.size*sd
        
        aa = mean + .5*mean
        bb = mean + .3*mean
        cc = aa - coeff
     
   mean.g2 = min(bb, cc)
   mean.g1 = mean.g2 + coeff
      
    TRUE.d = (mean.g1 - mean.g2) / sd
    }
    
    
    if(n1 < 2L) { n1 = 2L } ; if(n2 < 2L) { n2 = 2L }
    
    if(!paired & !researcher.tool){  
      
      y = as.vector(unlist(mapply(FUN = rnorm, n = c(n1, n2), mean = c(mean.g1, mean.g2), sd = c(sd, sd) )))
    }
    
    if(paired & !researcher.tool){
      
      cor.pop = correlation # true correlation between normal pops.
      
      
      mu <- c(0, 0)
      cov.pop <- matrix(c(1, cor.pop, cor.pop, 1), nrow = 2)
      
      if(!require(MASS)) install.packages("MASS")
      library(MASS)
      mvnorm.mat <- mvrnorm(n1, Sigma = cov.pop, mu = mu)
      
      a <- mvnorm.mat[ , 1] * sd + mean.g1
      b <- mvnorm.mat[ , 2] * sd + mean.g2
      
      y = c(a, b)
    }
    
    if(researcher.tool & !paired){ 
      
      message("\nNOTE: When using \"researcher.tool\" no effect size other than \"0\" could be used due to the Binomial data-generating processes in place." )
      
      p1 = .5 ; p2 = .5
      
      y = as.vector(unlist(mapply(FUN = rbinom, n = c(n1, n2), size = c(max.score, max.score), prob = c(p1, p2))))
    
      m1 = max.score*p1
      m2 = max.score*p2
      sd1 = sqrt(m1*(1-p1))
      sd2 = sqrt(m2*(1-p2))
      TRUE.d = (m1-m2) / sqrt((((n1 - 1)*((sd1^2))) + (n2 - 1)*((sd2^2)))/((n1+n2)-2))
    }
    
    
    if(researcher.tool & paired){
      
      message("\nNOTE: When using \"researcher.tool\" no effect size other than \"0\" could be used due to the Binomial data-generating processes in place." )
      
      a = function(rho, p, q) {
        rho * sqrt(p*q*(1-p)*(1-q)) + (1-p)*(1-q)
      }
      #
      # Specify the parameters.
      #
      trials = max.score
           p = .5
           q = .5
         rho = correlation
         
         m1 = max.score*p
         m2 = max.score*q
         sd1 = sqrt(m1*(1-p))
         sd2 = sqrt(m2*(1-q))
         TRUE.d = (m1-m2) / sqrt((((n1 - 1)*((sd1^2))) + (n2 - 1)*((sd2^2)))/((n1+n2)-2))
      
      a.0 <- a(rho, p, q)
      prob <- c(`(0,0)`= a.0, `(1,0)`=1-q-a.0, `(0,1)`= 1-p-a.0, `(1,1)`= a.0+p+q-1)
      if (min(prob) < 0) {
        print(prob)
        stop("Error: a probability is negative.")
      }
      #
      # Illustrate generation of correlated Binomial variables.
      #
      n <- n1
      u <- sample.int(4, n * trials, replace=TRUE, prob=prob)
      y <- floor((u-1)/2)
      x <- 1 - u %% 2
      bi1 <- colSums(matrix(x, nrow = trials)) # Sum in groups of `max.score`
      bi2 <- colSums(matrix(y, nrow = trials)) # Sum in groups of `max.score`
      
      y = c(bi1, bi2)
    }
    
    
    if(!paired) {
      
      groups = factor(rep(1:2, times = c(n1, n2)), labels = c("Treatment", "Control") )
      
    } else { 
      
      groups = factor(rep(1:2, times = c(n1, n2)), labels = c("Post-Test", "Pre-Test") )
    }
    
    
    mean.g1 = if(!paired)  mean(y[groups == "Treatment"])  else   mean(y[groups == "Post-Test"])  
    mean.g2 = if(!paired)  mean(y[groups == "Control"])    else   mean(y[groups == "Pre-Test"])    
    
    sd.g1 = if(!paired) sd(y[groups == "Treatment"]) else sd(y[groups == "Post-Test"])
    sd.g2 = if(!paired) sd(y[groups == "Control"])   else sd(y[groups == "Pre-Test"])
    
    
    groups.for.t = factor( rep(1:2, times = c(n1, n2)) )
    
    test = t.test(y ~ groups.for.t, var.equal = TRUE, paired = ifelse(paired, TRUE, FALSE))
    
    t.value = unname(test$statistic) ; p.value = test$p.value
    
    N = ifelse(paired, n1, (n1*n2)/(n1+n2) )
    
    Cohend = t.value / sqrt(N) 
    
    lab1 = if(n1 < 10 || n2 < 10) paste0("subj #", rev(1L:n1)) else c(paste0("subj #", rev(n1)[1]), paste0(rep(".", n1 - 2)), paste0("subj #", 1L))
    lab2 = if(n1 < 10 || n2 < 10) paste0("subj #", rev(1L:n2)) else c(paste0("subj #", rev(n2)[1]), paste0(rep(".", n2 - 2)), paste0("subj #", 1L))
    
    original.par = par(no.readonly = TRUE)
    on.exit(par(original.par))
    
    if(subjects) {
      
      par(font.lab = 2, xaxt = "n", ...)
      dotchart(y, groups = groups, color = c(4, 2)[groups], 
               font = 2, pch = 19, gcolor = c(4, 2), xlab = "Participants' Scores",
               pt.cex = ifelse(n1 <= 20 || n2 <= 20, 1.5, .8), labels = c(lab1, lab2), main = NA,
               cex.main = 2)
      
    } else {
      
      par(font.lab = 2, xaxt = "n", ...)
      dotchart(y, groups = groups, color = c(4, 2)[groups], 
               font = 2, pch = 19, gcolor = c(4, 2), xlab = "Participants' Scores",
               pt.cex = ifelse(n1 <= 20 || n2 <= 20, 1.5, .8), labels = NA, main = NA)
    }
    
par(xaxt = "s") ; if(researcher.tool) axis(1, at = min(y):max(y), font = 2) else axis(1, font = 2)

    gpos = rev(cumsum(rev(tapply(groups, groups, length)) + 2) - 1)
    
    decimal <- function(x, k){
      
      if(is.character(x)) {
        return(x)
      }else{
        format(round(x, k), nsmall = k, scientific =
                 ifelse(x >= 1e5 || x <= -1e5 || x <= 1e-5 & x >= -1e-5, TRUE, FALSE) )
      }
    }
    
    pars = par("usr")  
    
    segments( c(mean.g2, mean.g1), c(pars[3], pars[4]), c(mean.g2, mean.g1), rep( gpos[[2]], 2), lty = 2,
              col = c(2, 4) )
    
    arrows(mean.g2, gpos[[2]], mean.g1, gpos[[2]], code = 3, length = .08, col = "darkgreen")
    
    
    mean.diff = mean.g1 - mean.g2
    
    text( (mean.g1+mean.g2)/2, gpos[[2]], bquote(bold("Mean diff." == .(decimal((mean.diff), 2)))), font = 2, pos = 3, col = "green4", cex = 1.15 )
    
    legend("topright", legend = bquote(bold("Cohen's"~ bolditalic(d) == .(decimal(Cohend, 2)) )), bty = "n", text.col = "red4", cex = 1.15, bg = NA)
    legend("bottomright", legend = bquote("t.value" == .(as.numeric(t.value))))
    
    if(descriptives) {
      
      legend("topleft", legend = bquote(bold(Mean == .(decimal(mean.g1, 2)))), text.col = 4, bty = "n", bg = NA)
      
      legend("topleft", legend = bquote(bold(sd == .(decimal(sd.g1, 2)))), text.col = 4, bty = "n", bg = NA,
             inset = .03, adj =  c(.2, 0.5) )
      
      
      legend("bottomleft", legend = bquote(bold(Mean == .(decimal(mean.g2, 2)))), text.col = 2, bty = "n", bg = NA, 
             inset = .03, adj = .1)
      legend("bottomleft", legend = bquote(bold(sd == .(decimal(sd.g2, 2)))), text.col = 2, bty = "n", bg = NA,
             adj =  c(-.1, 0.5))
    }
    
    invisible(list(Cohend = Cohend, mean.diff = mean.diff, t.value = t.value, TRUE.d = TRUE.d, p.value = p.value))
  }
  
 
  if(simulation){
    
    if(missing(sim.time)) sim.time = .7
    
    
       Cohend.sim = numeric(n.sim)
    mean.diff.sim = numeric(n.sim)
      t.value.sim = numeric(n.sim)
    
    
    for(i in 1:n.sim){
      
      a = ttest.sim()
      
         Cohend.sim[i] <- a$Cohend
      mean.diff.sim[i] <- a$mean.diff
        t.value.sim[i] <- a$t.value
      
      Sys.sleep(sim.time)  
      
    }
    
    cbind("TRUE.Cohen.d" = a$TRUE.d, "Cohen.d" = decimal(Cohend.sim, 2), "Mean.diff" = decimal(mean.diff.sim, 2), 
          t.value = decimal(t.value.sim, 2), "Cohend.sd" = decimal(sd(Cohend.sim), 2), 
          "Mean.diff.sd" = decimal(sd(mean.diff.sim), 2),
          "Cohend.Range" = decimal(max(Cohend.sim) - min(Cohend.sim), 2), "Mean.diff.Range" = 
            decimal(max(mean.diff.sim) - min(mean.diff.sim), 2))
    
  } else {
    
    b = ttest.sim()

    library(ReporteRs)
    
    doc = addPlot(docx(), fun = ttest.sim, vector.graphic = TRUE, width = 5.5, height = 6,  
               par.properties = parCenter(), editable = TRUE)
    
     writeDoc(doc, file = "t.testsim3.docx" )
    
    
    list(t.value = b$t.value, Cohend = b$Cohend, TRUE.Cohen.d = b$TRUE.d, p.value = b$p.value)
    }
}


ttest(n1 = 30, n2 = 30, sim.time = .7, min.score = 0, max.score = 25, 
      simulation = FALSE, n.sim = 30, paired = FALSE,
      researcher.tool = FALSE, effect.size = .5)