i=4
a
mm <- n_u_minus_n[,i]+n[,i]
nn <- N_minus_n_t[,i] + n_t[,i] - mm
kk  <- n_t[,i]
BB  <- pmin(mm,kk)
stopifnot(n[,i]<=BB)
xx  <- pmin(n[,i],BB)

# get pv (right tail) and alpha star:
p[,i]  <- dhyper(x = xx, m = mm, n = nn, k = kk) + phyper(q = xx, m = mm, n = nn, k = kk, lower.tail = FALSE)
alpha_star[,i] <- dhyper(x=BB, m = mm, n = nn, k = kk)

# get realized randomized pv
f_n <- dhyper(x = xx, m = mm, n = nn, k = kk)
kappa <- apply(cbind(BB,mm,nn,kk,f_n),1,
               function(x) sum(dhyper(x = 0:x[1], m = x[2], n = x[3], k = x[4])==x[5]))
subset  <- rep(T,nrow(p))
m_tag <- sum(subset)
lpl <- (lambda < p[subset,i]) & (p[subset,i] <= lambda + kappa[subset]*f_n[subset])

EU_ecdf_lambda <- sum(p[subset,i]<=lambda)/m_tag +
                  sum(1- (p[lpl,i]-lambda) / (kappa[lpl]*f_n[lpl]) )/m_tag
       
EU_Pi_0_hat <- (1-EU_ecdf_lambda)/(1-lambda)

B <- 10000
p_rand_sim <- matrix(NA,nrow = nrow(n), ncol = B)
for (j in 1:ncol(p_rand_sim))
{
  u <- runif(n=length(f_n))
  p_rand_sim[,j] <- p[,i]-u*kappa*f_n
}
(1-mean(p_rand_sim<=lambda))/(1-lambda)



