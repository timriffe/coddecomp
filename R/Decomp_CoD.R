# This function tries to get the direct discrete life expectancy
# sensitivity to mx.
# in continous math it's -l(x)e(x), we just need to find the best
# approx with a discrete lifetable

arriaga_cof <- function(mx1,
                        mx2,
                        mx1_causes,
                        mx2_causes,
                        age = 0:(length(mx1) - 1),
                        sex1 = 't',
                        sex2 = sex1,
                        closeout = TRUE){
  ax1 <- mx_to_ax(mx = mx1,
                  age = age,
                  sex = sex1,
                  closeout = closeout)
  ax2 <- mx_to_ax(mx = mx2,
                  age = age,
                  sex = sex2,
                  closeout = closeout)
  qx1 <- mx_to_qx(mx = mx1,
                  ax = ax1,
                  age = age,
                  closeout = closeout)
  qx2 <- mx_to_qx(mx = mx2,
                  ax = ax2,
                  age = age,
                  closeout = closeout)
  lx1 <- qx_to_lx(qx1)
  lx2 <- qx_to_lx(qx2)
  dx1 <- lx_to_dx(lx1)
  dx2 <- lx_to_dx(lx2)
  Lx1 <- ald_to_Lx(ax = ax1,
                   lx = lx1,
                   dx = dx1)
  Lx2 <- ald_to_Lx(ax = ax2,
                   lx = lx2,
                   dx = dx2)
  Tx1 <- rcumsum(Lx1)
  Tx2 <- rcumsum(Lx2)
  
  direct = lx1 * (Lx2 / lx2 - Lx1 / lx1)
  indirect = lead(Tx2, default = 0) * (lx1 / lx2 - lead(lx1, default = 0) / lead(lx2, default = 0))
  N <- length(mx1)
  indirect[N] = lx1[N] * (Tx2[N] / lx2[N] - Tx1[N] / lx1[N])
  
  cc = direct + indirect
  
  cc2 <- (cc*(mx1_causes - mx2_causes))/(mx1 - mx2)
  
  return(cc2)
}

arriaga_instantaneous_cof <- function(mx1,
                                      mx2,
                                      mx1_causes,
                                      mx2_causes,
                                      age = 0:(length(mx1) - 1),
                                      sex1 = 't',
                                      sex2 = sex1,
                                      closeout = TRUE){
  mx_mean  <- (mx1 + mx2) / 2
  s_mean <- sen_arriaga_instantaneous(mx_mean, age = age)
  
  mx_causes <- mx2_causes - mx1_causes
  
  cc <- s_mean*mx_causes
  
  return(cc)
}

