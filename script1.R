# functions: --------------------------------------------------------------

# print a table:
get_222_ftable <- function(row)
{
  z <- array(c(row$n111,row$n211,row$n121,row$n221,row$n112,row$n212,row$n122,row$n222), dim=c(2, 2, 2))
  #addmargins(z) 
  dimnames(z) <- list("Deformed"=c("Y","N"),"Group"=c("KO","WT"),"Gender"=c("Female","Male"))
  return(ftable(z,row.vars = 1,col.vars = c(3,2)))
}

# get the pdf and support , input: table cells - see example
get_cond_test_pdf_supp <- function (df) 
{
  ## deformed in KO (F and M) ; white balls in hand
  xx1 <- df$n111
  xx2 <- df$n112
  ## deformed in KO and WT(F and M) : total white balls
  mm1 <- df$n111 + df$n121 # = n1.1
  mm2 <- df$n112 + df$n122 # = n1.2
  ## not deformed in KO and WT (F and  M)
  nn1 <- df$n211 + df$n221 # = n2.1
  nn2 <- df$n212 + df$n222 # = n2.2
  ## in KO deformed and not deformed(F and M) : total in hand
  kk1 <- df$n.11
  kk2 <- df$n.12
  
  BB1  <- pmin(mm1,kk1)
  BB2  <- pmin(mm2,kk2)
  AA1  <- pmax(0,kk1-nn1)
  AA2  <- pmax(0,kk2-nn2)
  
  pdf_list <- supp_list <- list()
  for (i in 1:nrow(df))
  {
    x <- dhyper(x = AA1[i]:BB1[i], m = mm1[i], n = nn1[i], k = kk1[i])
    y <- dhyper(x = AA2[i]:BB2[i], m = mm2[i], n = nn2[i], k = kk2[i])
    pdf_list[[i]] <- convolve(x, rev(y), type = "open")
    supp_list[[i]] <- (AA1[i]+AA2[i]):(BB1[i]+BB2[i])
  }
  
  alpha_star <- sapply(pdf_list,function(x) rev(x)[1], simplify = TRUE)
  
  return(data.frame(pdf_list=I(pdf_list),supp_list=I(supp_list),alpha_star=alpha_star))
}

# get the p value, the last pdf value i the p-value and if the value repeat itself more than once (kappa). 
get_cond_test_pvs <- function(pdf_list,supp_list,xx1,xx2)
{
  xx <- xx1+xx2
  ind_pv <- mapply(function(supp,xx) which(rev(supp)==(xx))[1],supp_list,xx)
  p   <-   mapply(function(pdf,ind_pv) sum(rev(pdf)[1:ind_pv]), pdf_list, ind_pv)
  f_n <-   mapply(function(pdf,ind_pv) rev(pdf)[ind_pv], pdf_list, ind_pv)
  kappa <- mapply(function(pdf,f_n) sum(pdf==f_n), pdf_list, f_n)   
  return(data.frame(p = p,f_n = f_n,kappa = kappa))
}

# example -----------------------------------------------------------------
xx <- tbl_sim[which(tbl_sim$Outcome==u[b] & tbl_sim$p<=0.1)[2],]
cbind(xx$supp_list[[1]],xx$pdf_list[[1]])
plot(xx$supp_list[[1]],xx$pdf_list[[1]])

# the table with the originally observed values 
get_222_ftable(xx)

# pdf and supp
pdf_supp <- get_cond_test_pdf_supp(xx)

# p value and f_n (the most inner pdf value in the p-value)
get_cond_test_pvs(pdf_list = pdf_supp$pdf_list,supp_list = pdf_supp$supp_list,
                  xx1 = xx$n111_sim,xx2 = xx$n112_sim)

# finally, in the swap() function - which is equvivalent to your fdr.threshold -
# there is the p_swap computation as follows:

p_swap <- 1 - p + f_n 



