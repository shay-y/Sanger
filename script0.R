rm(list=ls())

fdr.threshold1 = function (pv, pv.swapped, fdr = 0.05) 
{

  ts = (pmin(pv, pv.swapped)) #CAN MAKE IT MORE EFFICIENT SINCE TAKE ONLY MAX AND UNIQUE
#   eps = min(min(ts/2),1e-16)
#   ts = ts+eps
  ratio = sapply(ts, function(t)((1+sum(pv.swapped<=t & pv.swapped<=pv))/max(1,sum(pv<=t & pv<pv.swapped))))
  ok = which(ratio<=fdr)
  
  ifelse(length(ok) > 0, max(ts[ok]), -Inf)
}

#with binom pi0=0.6, n=20, swapping more powerful 0.8 versus 0.6

sim1 =function(pi0,p1,m=1000,B = 1000)
{

R = NULL
S = NULL
fdp = NULL
R.BH = NULL
S.BH= NULL
fdp.BH = NULL
is.null = rep(0,m)
is.null[1:(pi0*m)]= rep(1,pi0*m)
#lambda = 0.3
alpha=0.05
#mu = rep(0,m)
#mu[is.null==0] = 2
n=20 
m0 =m*pi0
for ( b in 1:B){
x= c(rbinom(m0,n,prob=0.5),rbinom(m-m0,n,prob=p1))
pv = 1-pbinom(x-1,n,1/2 )
pv.swapped = pbinom(x,n,1/2 )
#head(cbind(x, pv, pv.swapped, pv+pv.swapped-1, dbinom(x, n, 0.5)))
out = fdr.threshold1(pv,pv.swapped, fdr = alpha)
print(c(b,out))
R[b] = sum(pv<=out & pv<pv.swapped)
S[b]= sum(pv[is.null==0]<=out & pv[is.null==0]<pv.swapped[is.null==0])
print(c(b,out, R[b], S[b]))

fdp[b] = (R[b]-S[b])/max(1,R[b])

out = p.adjust(pv, method = "BH")
R.BH[b] = sum(out<=alpha)
S.BH[b]= sum(out[is.null==0]<=alpha)
fdp.BH[b] = (R.BH[b]-S.BH[b])/max(1,R.BH[b])
}
return(list(R=R, S=S, fdp = fdp,R.BH=R.BH, S.BH=S.BH, fdp.BH = fdp.BH, pi0=pi0, p1=p1, m=m ))
}

m=1000
pi0 =  0.8
p1= 0.7
out =sim1(pi0,p1, m ,B = 1000)


##############################################################################

fdr.threshold2 = function (X, X.swapped, fdr = 0.05) 
{
  ts = (pmax(X, X.swapped)) #CAN MAKE IT MORE EFFICIENT SINCE TAKE ONLY MAX AND UNIQUE
  #  ratio = sapply(ts, function(t)((1+sum(pv.swapped<=t ))/max(1,sum(pv<=t ))))
  ratio = sapply(ts, function(t)((1+sum(X.swapped>=t & X.swapped>=X))/max(1,sum(X>=t & X>X.swapped))))
  ok = which(ratio<=fdr)
  
  ifelse(length(ok) > 0, min(ts[ok]), Inf)
}

swap_x <-  function(x,n,q = 0.05)
{
  x_swap <- n - x 
  S <- x_swap >= x
  t_vec <- unique(sort(pmax(x,x_swap)))
  t_min <- -1
  for (ti in 1:length(t_vec))
  {
    if ( (1 + sum(x_swap[S]>=t_vec[ti]))/max(sum(x[!S]>=t_vec[ti]),1) <= q )
    {
      t_min <- t_vec[ti]
      break
    }
  }
  R <- rep(F,length(p))
  R[!S] <- x[!S]>=t_min
  return(R)
}

#with binom pi0=0.6, n=20, swapping more powerful 0.8 versus 0.6

simDebug =function(pi0,p1,m=1000,B = 1000)
{
  
  R = NULL
  S = NULL
  fdp = NULL
  R.BH = NULL
  S.BH= NULL
  fdp.BH = NULL
  is.null = rep(0,m)
  is.null[1:(pi0*m)]= rep(1,pi0*m)
  #lambda = 0.3
  alpha=0.05
  #mu = rep(0,m)
  #mu[is.null==0] = 2
  n=20 
  m0 =m*pi0
  for ( b in 1:B){
    x= c(rbinom(m0,n,prob=0.5),rbinom(m-m0,n,prob=p1))
   x.swapped = n-x
    # pv = 1-pbinom(x-1,n,1/2 )
  #  pv.swapped = pbinom(x,n,1/2 )
    #head(cbind(x, pv, pv.swapped, pv+pv.swapped-1, dbinom(x, n, 0.5)))
    out = fdr.threshold2(x,x.swapped, fdr = alpha)
    print(c(b,out))
    R[b] = sum(x>=out & x>x.swapped)
    S[b]= sum(x[is.null==0]>=out & x[is.null==0]>x.swapped[is.null==0])
    print(c(b,out, R[b], S[b]))
    
    fdp[b] = (R[b]-S[b])/max(1,R[b])
  }
  return(list(R=R, S=S, fdp = fdp, pi0=pi0, p1=p1, m=m ))
}

m=1000
pi0 =  0.8
p1= 0.7
outDebug =simDebug(pi0,p1, m ,B = 1000)
summary(outDebug$fdp)
head(outDebug$fdp)
summary(out$fdp.BH)

#  ------------------------------------------------------------------------
  pi0 =  0.8
  p1= 0.7

  is.null = rep(0,m)
  is.null[1:(pi0*m)]= rep(1,pi0*m)
  alpha=0.05
  n=20 
  m0 =m*pi0
  x= c(rbinom(m0,n,prob=0.5),rbinom(m-m0,n,prob=p1))
  pv = 1-pbinom(x-1,n,1/2 )
  pv.swapped = pbinom(x,n,1/2 )
  
  x.swapped = n-x
  (out1 = fdr.threshold1(pv,pv.swapped, fdr = alpha))
  (out2 = fdr.threshold2(x,x.swapped, fdr = alpha))
  
  (sum(pv<=out1 & pv<pv.swapped))
  (sum(pv[is.null==0]<=out1 & pv[is.null==0]<pv.swapped[is.null==0]))
  
  (sum(x>=out2 & x>x.swapped))
  (sum(x[is.null==0]>=out2 & x[is.null==0]>x.swapped[is.null==0]))
  
fdr = 0.05
ts1 = (pmin(pv, pv.swapped))
ratio1 = sapply(ts1, function(t)((1+sum(pv.swapped<=ts1[293]) & pv.swapped<=pv))/max(1,sum(pv<=ts1[293] & pv<pv.swapped))))
which(ratio1<=fdr)

r1 = sapply(ts1, function(t)(sum(pv<=t & pv<pv.swapped)))
r2 = sapply(ts2, function(t)(sum(x>=t & x>x.swapped)))
which(r1!=r2)

sum(pv<=ts1[39]) & pv<pv.swapped)
sum(x>=ts2[39]) & x>x.swapped)
pmin(pv[39], pv.swapped[39])

print(cbind(x.swapped,ts2[293],pv.swapped,pv[293])[x.swapped>=ts2[293],][1,],digits=20)
ts2[293]
pv.swapped = pbinom(x,n,1/2)

pv
print(ts1[293],digits=20)

  