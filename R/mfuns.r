################################################################################
dpmn=function(kk=NULL, pp)
{
    dyn.load("pm-fft.so")
    mm=ncol(pp) # m categories
    nn=nrow(pp) # n people

    nn.vec=rep(nn+1, mm-1)
    l.vec=rep(0, mm-1)
    cn.vec=cumprod(nn.vec)
    cn.vec=c(1, cn.vec[-(mm-1)])
    cn.vec=cn.vec[length(cn.vec):1]
    cn.vec=as.integer(cn.vec)

    nnt=prod(nn.vec)
    res0=double(nnt)

  #browser()

  tmp=.C("pmn_mdfft", as.double(res0), as.integer(nnt), as.integer(nn), as.integer(mm), as.double(pp), as.integer(nn.vec), as.integer(l.vec), as.integer(cn.vec))
  res0=tmp[[1]]
  #example an_array[k + 27 * (j + 12 * i)]
  #print(round(res0, 9))

  res=array(0, nn.vec)

  res.expr="res[idx[1]"
  if(mm>=3)
  {
    for(i in 2:(mm-1))
    {
       res.expr=paste0(res.expr, ", idx[", i, "]")
    }
  }
  res.expr=paste0(res.expr, "]=res0[i]")

  #browser()

  #print(nnt)

  for(i in 1:nnt)
  {
    idx=l.vec.compute(k=i, cn.vec=cn.vec, m=mm)
    #print(idx)
    eval(parse(text=res.expr))
  }

  res=round(res, 10)
  return(res)
}
################################################################################
l.vec.compute=function(k, cn.vec, m)
{
  k=k-1
  l.vec=rep(0, m-1)
  for(i in 1:(m-1))
  {
    aa=k%%cn.vec[i]
    bb=(k-aa)/cn.vec[i]
    l.vec[i]=bb
    k=aa
  }
  l.vec=l.vec+1
  return(l.vec)
}
################################################################################
pmatrix <- function(n,m){
  p <- matrix(0,nrow = n,ncol = m,byrow = T)
  for (i in 1:n) {
    r <- runif(m)
    r <- r/sum(r) #generate row
    p[i,] <- r
  }
  return(p)
}
##################################################################################
dpm_sim = function(kk=NULL,pp,t){
  dyn.load("simulation.so")
  mm=ncol(pp) # m categories
  nn=nrow(pp) # n people
  nn.vec=rep(nn+1, mm-1)
  l.vec=rep(0, mm-1)
  cn.vec=cumprod(nn.vec)
  cn.vec=c(1, cn.vec[-(mm-1)])
  cn.vec=cn.vec[length(cn.vec):1]
  cn.vec=as.integer(cn.vec)
  nnt=prod(nn.vec)
  res0=double(nnt)
  temp=.C("pmd_simulation", as.double(res0),
          as.integer(nnt) , as.integer(nn), as.integer(mm), as.double(pp),
          as.integer(nn.vec), as.integer(l.vec), as.integer(cn.vec), as.double(t))
  res=round(temp[[1]],10)
  return(res)
}
#####################################################################################
dpm_noraml = function(kk=NULL,pp,x_vec){
  dyn.load("normal_appro.so")
  mm=ncol(pp) # m categories
  mm = mm - 1
  nn=nrow(pp) # n people
  x_vec = x_vec[1:(length(x_vec)-1)]

  temp=.C("pmd_normal", res0 = as.double(0) ,as.integer(nn), as.integer(mm), as.double(pp),as.integer(x_vec))
  return(temp$res0)
}
