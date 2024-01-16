# This function tries to get the direct discrete life expectancy
# sensitivity to mx.
# in continous math it's -l(x)e(x), we just need to find the best
# approx with a discrete lifetable
sen_e0_mx <- function(mx){
  # ax <-
  lx <- mx_to_lx(mx)
  dx <- lx_to_dx(lx)
  qx <- dx / lx
  ax <- mq_to_ax(mx,qx)
  Lx <- ald_to_Lx(ax,lx,dx)
  ex <- lL_to_ex(lx,Lx)
  # ex <- mx_to_ex(mx)
  N <- length(mx)
  x <- (1:N) - 1 
  #ex2 <- mean2(ex)
  ex2 <- approx(x,ex,xout = x + ax)$y
  Lx2 <- approx(x,Lx,xout = x + ax)$y
  dx2 <- approx(x,dx,xout = x + ax)$y
  ax2 <- approx(x,ax,xout = x + ax)$y
  ex2[N] <- ex[N]
  -Lx * ex2
  #-Lx2 * ex2 - dx2*(1-ax2) *ex2 #+ direct
}

# direct = (lx - dx * ax)
# plot(direct[1:100] - Lx[1:100])
# plot(sen_e0_mx(mx))
# plot(sen_e0_mx(mx) - numDeriv::grad(mx_to_e0, mx))
# 
# plot(sen_e0_mx(mx) - numDeriv::grad(mx_to_e0, mx))
# lines(numDeriv::grad(mx_to_e0, mx))
# 
# 
# plot(dx * approx(x=0:100,ex,xout = 0:100 + ax)$y)



# e0_qx_simplified <- function(qx){
#   N  = length(qx)
#   px = 1 - qx
#   lx = c(1, cumprod(px),0)
#   # Lx = (lx + lead(lx, default = 0)) / 2
#   # Lx[N] = sum(Lx[N:length(Lx)])
#   # Lx <- Lx[1:N]
#   sum(lx)
# }
# 
# sen_arriaga_qx_simplified <- function(qx1,qx2){
#   N  = length(qx1)
#   px1 = 1 - qx1
#   lx1 = c(1, cumprod(px1),0)
#   Lx1 = (lx1 + lead(lx1, default = 0)) / 2
#   Lx1[N] = sum(Lx1[N:length(Lx1)])
#   Lx1 = Lx1[1:N]
#   Tx1 = rev(cumsum(rev(Lx1)))
#   lx1 = lx1[1:N]
#   ex1 = Tx1 / lx1
#   dx1 = qx1 * lx1
#   
#   px2 = 1 - qx2
#   lx2 = c(1, cumprod(px2),0)
#   Lx2 = (lx2 + lead(lx2, default = 0)) / 2
#   Lx2[N] = sum(Lx2[N:length(Lx2)])
#   Lx2 = Lx2[1:N]
#   Tx2 = rev(cumsum(rev(Lx2)))
#   lx2 = lx2[1:N]
#   ex2 = Tx2 / lx2
#   dx2 = qx2 * lx2
#   
#   #direct = lx1 * (Lx2 / lx2 - Lx1 / lx1)
#   direct = lx1 * ((lx2 - dx2 * .5) / lx2 - (lx1 - dx1 * .5)/ lx1)
#   N = length(mx1)
#   indirect = 
#     lead(Tx2, default = 0) * 
#     (lx1 / lx2 - lead(lx1, default = 0) / 
#        lead(lx2, default = 0))
#   indirect[N] = lx1[N] * (ex2[N] - ex1[N])
#   #indirect = ifelse(is.na(indirect),0,indirect),
#   cc = direct + indirect
#   delta = qx2 - qx1
#   # watch out for 0s in delta denominator
#   sen = cc / delta 
#   sen
# }
# 
# # well, pseudo-instantaneous anyway, and symmetrical 
# # This is designed to work with ltre()
# sen_arriaga_instantaneous_qx_simplified <- function(qx, perturb = 1e-6){
#   qx1 <- qx * (1 / (1 - perturb))
#   qx2 <- qx * (1 - perturb)
#   (sen_arriaga_qx_simplified(qx1 = qx1, qx2 = qx2) + 
#       sen_arriaga_qx_simplified(qx1 = qx2, qx2 = qx1 )) / 2
# }
# logit <- function(x){
#   log(x/(1-x))
# }
# expit <- function(x){
#   exp(x) / (1+exp(x))
# }
# 
# sen_arriaga_instantaneous_qx_simplified_logit <- function(qx, perturb = 1e-6){
#   qx1 <- expit(logit(qx) + perturb)
#   qx2 <- expit(logit(qx) - perturb)
#   (sen_arriaga_qx_simplified(qx1 = qx1, qx2 = qx2) + 
#       sen_arriaga_qx_simplified(qx1 = qx2, qx2 = qx1 )) / 2
# }
# 
#  qx <- readHMDweb("USA","mltper_1x1",password = Sys.getenv("pw"), username = Sys.getenv("us")) |> 
#    filter(Year == 2000) |> 
#    pull(qx) %>%
#    '['(1:101)
#  
# sai  <- sen_arriaga_instantaneous_qx_simplified(qx)
# # sail <- sen_arriaga_instantaneous_qx_simplified_logit(qx)
# sn  <-  numDeriv::grad(fun=e0_qx_simplified, qx)
# 
# sen_e0_qx_simplified <- function(qx){
#   N  = length(qx)
#   N1 = N+1
#   px = 1 - qx
#   lx = c(1, cumprod(px),0)
#   Lx = (lx + lead(lx, default = 0)) / 2
#   Lx[N1] = sum(Lx[N1:length(Lx)])
#   Lx = Lx[1:N1]
#   ex = rev(cumsum(rev(Lx))) / lx[1:N1]
#   ex2 = (ex + lead(ex, default = 0)) / 2
#   ex2 = ex2[1:N]
#   lx = lx[1:N]
#   Lx[N] = sum(Lx[N:N1])
#   Lx = Lx[1:N]
#   dx = qx * lx
# 
#   -(lx[1:N]+dx/2) * ex2[1:N] 
#   #-Lx[1:N]  * ex[1:N]
#   #-lx[1:N] * ex[1:N]
# }
# sqs = sen_e0_qx_simplified(qx)
# resid = sqs - sn
# plot(resid)
# lines(
#        (lx[1:N] - px * lx[1:N])
#      )
# plot(sai - sn[1:100])

