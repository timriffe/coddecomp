# This function tries to get the direct discrete life expectancy
# sensitivity to mx.
# in continous math it's -l(x)e(x), we just need to find the best
# approx with a discrete lifetable

#' @title A direct approximation of the sensitivity of life expectancy at birth to changes in mortality.
#' @description This direct lifetable-based calculation requires a few approximations to get a usable value whenever we're working with discrete data. In continuous notation, we know that the sensitivity \eqn{s(x)}
#' \deqn{s(x) = -l(x)e(x)}
#' but it is not obvious what to use from a discrete lifetable. In this implementation, we use \eqn{L(x)} and an \eqn{a(x)}-weighted average of successive \eqn{e(x)} values, specifically, we calculate:
#' \deqn{s(x) = -L(x) \cdot \left e(x) \cdot (1 - a(x)) + e(x+1) * a(x)\right)}
#' This seems to be a very good approximation for ages >0, but we still have a small, but unaccounted-for discrepancy in age 0, at least when comparing with also-imperfect numerical derivatives.
#' @importFrom data.table shift
#' @export
#' @examples
#' x <- 0:100
#' mx <- 0.001 * exp(x * 0.07)
#' sl <-  sen_e0_mx_lt(mx,age=x,sex='t',closeout=TRUE)
#' sn <- numDeriv::grad(mx_to_e0, mx, age=x, sex = 't', closeout=TRUE)
#' \dontrun{
#' plot(x,sl)
#' lines(x,sn)
#' }
#' # examine residuals:
#' sl - sn
#' # Note discrepancies in ages >0 are due to numerical precision only
#' \dontrun{
#' plot(x, sl - sn, main = "still uncertain what accounts for the age 0 discrepancy")
#' }
sen_e0_mx_lt <- function(mx, 
                      age = 0:(length(mx)-1), 
                      sex = 't',
                      closeout = TRUE){
  ax <- mx_to_ax(mx = mx, 
                 age = age, 
                 sex = sex, 
                 closeout = closeout)
  qx <- mx_to_qx(mx, ax)
  lx <- qx_to_lx(qx)
  dx <- lx_to_dx(lx)
  Lx <- ald_to_Lx(ax = ax,
                  lx = lx,
                  dx = dx)
  ex <- lL_to_ex(lx = lx, Lx = Lx)
  # ex <- mx_to_ex(mx)
  N <- length(mx)

  # This is the current-best approximation,
  # but still requires an unknown adjustment for age 0
  exs <- shift(ex, n = -1, fill = ex[N]) 
  ex2 <- ex * (1 - ax) + exs * ax

  ex2[N] <- ex[N]
 
  sen <- -Lx * ex2 

  sen
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

