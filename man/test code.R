dyn.load("normal_appro.so")
  pp=matrix(c(.1, .2, .7, .2, .5, .3, .7, .1, .2, .5, .2, .3), nrow=4, byrow=T)
pp <- pp[,-3]
mm=ncol(pp) # m categories
nn=nrow(pp) # n people
x_vec <- c(10,10,10,10,10,10,10,10,10)
temp=.C("pmd_normal", as.integer(nn), as.integer(mm), as.double(pp),as.integer(x_vec))
  
mm <- 10
nn <- 100
pp <- pmatrix(nn,mm)
pp <- pp[,-10]
mm=ncol(pp) # m categories
nn=nrow(pp) # n people
x_vec <- c(10,10,10,5,10,10,10,15,1)
temp=.C("pmd_normal", as.integer(nn), as.integer(mm), as.double(pp),as.integer(x_vec))



pp <- pmatrix(5,5)
dpmn(NULL,pp)
