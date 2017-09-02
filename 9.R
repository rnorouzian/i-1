d.CI.sim <- function(d, n1, n2 = NA, conf.level = .95, n.sim = 5, ylabel = FALSE){
  
CI.d <- function(da = d, n1a = n1, n2a = n2, conf.levela = conf.level){
  
alpha = (1 - conf.levela)/2
    N = ifelse(is.na(n2a), n1a, (n1a * n2a)/(n1a + n2a))
   df = ifelse(is.na(n2a), n1a - 1, (n1a + n2a) - 2)
 d.SE = 1/sqrt(N)  ;   t = da/d.SE
    
f <- function (ncp, alpha, q, df) {
    abs(suppressWarnings(pt(q = t, df = df, ncp, lower.tail = FALSE)) - alpha)
  }
    
a = lapply(14:ifelse(da!= 0, da*sqrt(N)+5, 30), function(x) c(-x, x))
    
CI = matrix(NA, length(a), 2)
    
for(i in 1:length(a)){    
CI[i,] = sapply(c(alpha, 1-alpha),
function(x) optimize(f, interval = a[[i]], alpha = x, q = t, df = df, tol = 1e-12)[[1]]*d.SE)
  }  
    
CI[which.max(ave(1:nrow(CI), do.call(paste, round(data.frame(CI), 3)), FUN = seq_along)), ]
  }
    
fun <- function(db = d, n1b = n1, n2b = n2, conf.levelb = conf.level){ 
  
   N = ifelse(is.na(n2b), n1b, (n1b * n2b)/(n1b + n2b))
  df = ifelse(is.na(n2b), n1b - 1, (n1b + n2b) - 2)
d.SE = 1/sqrt(N)  ;   t = db/d.SE
    
  ds = rt(1, df, t)*d.SE
c(CI = CI.d(da = ds, n1a = n1, n2a = n2, conf.levela = conf.level), ds = ds)
  }
  
  sim <- t(replicate(n.sim, fun()))
  capture = sim[ ,1] <= d & d <= sim[ ,2]
  
  original.par = par(no.readonly = TRUE) ; on.exit(par(original.par))
  par(mgp = c(2, .2, 0), tck = -.015)  
                
  plot(sim[, 1:2], rep(1:n.sim, 2), ty = "n", ylab = NA, yaxt = "n", xlab = "Effect Size", font.lab = 2)
  axis(1, at = d, col.axis = 2, col = 2, font = 2)
  abline(h = 1:n.sim, col = 8, lty = 3)
  if(ylabel) axis(2, at = 1:n.sim, labels = paste0("Repeat ", rev(1:n.sim)), font = 2, las = 1, cex.axis = .8, tck = -.006)
  abline(v = d, lty = 2, col = 2) 
  segments(sim[ ,1], 1:n.sim, sim[ ,2], 1:n.sim, lend = 1, col = ifelse(capture, 1, 2))
  points(sim[, 3], 1:n.sim, pch = 19, col = ifelse(capture, 1, 2), cex = ifelse(n.sim > 50, .6, .65))
  
  noquote(paste0("Coverage = ", mean(capture)*1e2, "%")) 
}
# Example of use:
d.CI.sim(d = .5, n1 = 20, n.sim = 20, ylabel = TRUE)
