
pre.c = rnorm(100)

pos.c = rnorm(100)

pre.t = rnorm(100)

pos.t = rnorm(100, 5)


data.frame(subjects = rep(1:200, each = 2), 
                  y = c(pre.c, pos.c, pre.t, pos.t), 
                  time = rep(0:1, 200), 
                  group = rep(c(0, 1), each = 200))
