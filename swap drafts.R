

p <- c(0.01,0.02,0.1,0.5,0.9,1)
f_n <- c(0.005,0.017,0.02,0.02,0.1,0.1)
q <- 0.05

swap <-  function(p,f_n,q = 0.05)
{
  p_swap <- 1 - p + f_n 
  S <- p_swap <= p
  t_vec <- p
  t_vec[S] <- p_swap[S]
  t_vec <- unique(sort(t_vec))
  t_max <- -1
  for (ti in length(t_vec):1)
  {
    if ( (1 + sum(p_swap[S]<=t_vec[ti]))/max(sum(p[!S]<=t_vec[ti]),1) <= q )
    {
      t_max <- t_vec[ti]
      break
    }
  }
  R <- rep(F,length(p))
  R[!S] <- p[!S]<=t_max
  return(R)
}


r <- swap(p = tbl_p$p,f_n = tbl_p$f_n)
sum(r)

(1 + sum(ps[S]<=t))/max(sum(p[!S]<=t),1)<=q
sum(p[!S]<=t) 

rank(t_vec[S])
o <- order(t_vec)
t_vec[o]                           # 0.01 0.02 0.10 0.10 0.20 0.50
p_rank      <- cumsum(S[o])    # 0 0 1 2 2 2 ; 0 0 0 1 2 2 ; 0 0 1 1 2 2 
p_swap_rank <- cumsum((!S)[o]) # 1 2 2 2 3 4 ; 1 2 3 3 3 4 ; 1 2 3 3 3 4 

t_vec <- (1 + p_rank) / pmax(1,p_swap_rank)        # 1.00 0.50 1.00 1.50 1.00 0.75
                                          # 1.0000000 0.5000000 0.3333333 0.6666667 1.0000000 0.7500000
R <- rep(F,length(p))
if (any(t_vec <= q))
{
  TT <- t_vec[max(which(t_vec <= q))]
  R[!S] <- p[!S]<=TT
}
  
return(R)


p_p_swap[o]

