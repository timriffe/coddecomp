###' This script contains functions that implement and adapt the 
###' Arriaga cause-of-death life expectancy decomposition approach
###'
###' @title Arriaga cause-of-death life expectancy decomposition approach
###' @description Following the notation given in Preston et al (2000), Arriaga's ##cause-of-death life expectancy decomposition method can written as:
###' \deqn{_{n}\Delta_{x}^{i} = _{n}\Delta_{x} \frac{ _{n}m^{i}_{x2} - ##_{n}m^{i}_{x1} }{ _{n}m_{x2} - _{n}m_{x1}  } }
###' where \eqn{_{n}\Delta_{x}} is the contribution of rate differences in age ##\eqn{x} to the difference in life expectancy implied by `mx1` and `mx2`. 
###' Additionally, `mx1_causes` and `mx2_causes` corresponds to a matrix with the ##mortality rates of every considered causes-of-death that the user wants to ##decompose. 
###' @details A little-known property of this decomposition method is that it is ##directional, in the sense that we are comparing a movement of `mx1` to `mx2`, and ##this is not exactly symmetrical with a comparison of `mx2` with `mx1`. 
###' Note also, if decomposing in reverse from the usual, you may need a slight ##adjustment to the closeout value in order to match sums properly (see examples ##for a demonstration).
###' 
###' @param mx1 numeric vector of the mortality rates (central death rates) for ##population 1
###' @param mx2 numeric vector of the mortality rates (central death rates) for ##population 2
###' @param age integer vector of the lower bound of each age group (currently only ##single ages supported)
###' @param mx1_causes numeric matrix of the cause-of-death mortality rates ##(central death rates) for population 1  
###' @param mx2_causes numeric matrix of the cause-of-death mortality rates ##(central death rates) for population 2
###' @param sex1 character either the sex for population 1: Male (`"m"`), Female ##(`"f"`), or Total (`"t"`)
###' @param sex2 character either the sex for population 2: Male (`"m"`), Female ##(`"f"`), or Total (`"t"`) assumed same as `sex1` unless otherwise specified.
###' @param closeout logical. Default `TRUE`. Shall we use the HMD Method Protocol ##to close out the `ax` and `qx` values? See details.
###' @details setting `closeout` to `TRUE` will result in value of `1/mx` for the ##final age group, of `ax` and a value of 1 for the closeout of `qx`.
###' @return `cc` numeric vector with one element per age group, and which sums to ##the total difference in life expectancy between population 1 and 2.
###' @importFrom data.table shift
###' @export
###' @references
###' \insertRef{arriaga1989changing}{coddecomp}
###' \insertRef{preston2000demography}{coddecomp}
###' @examples
###' a <- .001
###' b <- .07
###' x <- 0:100
###' mx1 <- a * exp(x * b)
###' mx2 <- a/2 * exp(x * b)
###' mx1_causes <- matrix(c(mx1/4, mx1/4, mx1/2), length(x), 3) 
###' mx2_causes <- matrix(c(mx2/3, mx2/3, mx2/3), length(x), 3)
###' 
###' cod1 <- arriaga_cof(mx1 = mx1, mx2 = mx2, age = x
###'             mx1_causes = mx1_causes, 
###'             mx2_causes = mx2_Causes)
###' e01 <- mx_to_e0(mx1, age = x)
###' e02 <- mx_to_e0(mx2, age = x)
###' 
###' (delta <- e02 - e01)
###' sum(cod1)
###'
###' # asymmetrical with a decomposition in the opposite direction
###' cod2 <- -arriaga_cof(mx1 = mx2, mx2 = mx1, age = x
###'                      mx1_causes = mx2_causes, 
###'                      mx2_causes = mx1_Causes)
###' sum(cod2);sum(ccod1);sum(delta)
###' 
###' #We shoudl include a function to plot the decomposition methods through the ##causes.
###' 
##arriaga_cof <- function(mx1,
##                        mx2,
##                        mx1_causes,
##                        mx2_causes,
##                        age = 0:(length(mx1) - 1),
##                        sex1 = 't',
##                        sex2 = sex1,
##                        closeout = TRUE){
##  ax1 <- mx_to_ax(mx = mx1,
##                  age = age,
##                  sex = sex1,
##                  closeout = closeout)
##  ax2 <- mx_to_ax(mx = mx2,
##                  age = age,
##                  sex = sex2,
##                  closeout = closeout)
##  qx1 <- mx_to_qx(mx = mx1,
##                  ax = ax1,
##                  age = age,
##                  closeout = closeout)
##  qx2 <- mx_to_qx(mx = mx2,
##                  ax = ax2,
##                  age = age,
##                  closeout = closeout)
##  lx1 <- qx_to_lx(qx1)
##  lx2 <- qx_to_lx(qx2)
##  dx1 <- lx_to_dx(lx1)
##  dx2 <- lx_to_dx(lx2)
##  Lx1 <- ald_to_Lx(ax = ax1,
##                   lx = lx1,
##                   dx = dx1)
##  Lx2 <- ald_to_Lx(ax = ax2,
##                   lx = lx2,
##                   dx = dx2)
##  Tx1 <- rcumsum(Lx1)
##  Tx2 <- rcumsum(Lx2)
##  
##  direct = lx1 * (Lx2 / lx2 - Lx1 / lx1)
##  indirect = lead(Tx2, default = 0) * (lx1 / lx2 - lead(lx1, default = 0) / lead##(lx2, default = 0))
##  N <- length(mx1)
##  indirect[N] = lx1[N] * (Tx2[N] / lx2[N] - Tx1[N] / lx1[N])
##  
##  cc = direct + indirect
##  
##  cc2 <- (cc*(mx1_causes - mx2_causes))/(mx1 - mx2)
##  
##  return(cc2)
##}
###' This script contains the Instantaneous Arriaga cause-of-death
###' life expectancy decomposition 
###'
###' @title Instantaneous Arriaga cause-of-death life expectancy decomposition ##approach
###' @description We present the next method to decompose the contribution of every ##cause of death through the ages to the differences in the life expectancy:
###' \deqn{_{n}\Delta_{x}^{i} = _{n}\Delta_{x} \frac{ _{n}m^{i}_{x2} - ##_{n}m^{i}_{x1} }{ _{n}m_{x2} - _{n}m_{x1}  } }
###' where \eqn{_{n}\Delta_{x}} is the contribution of rate differences in age ##\eqn{x} to the difference in life expectancy implied by `mx1` and `mx2` and ##obtained from the Arriaga-Instantaneous. 
###' Additionally, `mx1_causes` and `mx2_causes` corresponds to a matrix with the ##mortality rates of every considered causes-of-death that the user wants to ##decompose. 
###' @details A little-known property of this decomposition method is that it is ##directional, in the sense that we are comparing a movement of `mx1` to `mx2`, and ##this is not exactly symmetrical with a comparison of `mx2` with `mx1`. 
###' Note also, if decomposing in reverse from the usual, you may need a slight ##adjustment to the closeout value in order to match sums properly (see examples ##for a demonstration).
###' 
###' @param mx1 numeric vector of the mortality rates (central death rates) for ##population 1
###' @param mx2 numeric vector of the mortality rates (central death rates) for ##population 2
###' @param age integer vector of the lower bound of each age group (currently only ##single ages supported)
###' @param mx1_causes numeric matrix of the cause-of-death mortality rates ##(central death rates) for population 1  
###' @param mx2_causes numeric matrix of the cause-of-death mortality rates ##(central death rates) for population 2
###' @param sex1 character either the sex for population 1: Male (`"m"`), Female ##(`"f"`), or Total (`"t"`)
###' @param sex2 character either the sex for population 2: Male (`"m"`), Female ##(`"f"`), or Total (`"t"`) assumed same as `sex1` unless otherwise specified.
###' @param closeout logical. Default `TRUE`. Shall we use the HMD Method Protocol ##to close out the `ax` and `qx` values? See details.
###' @details setting `closeout` to `TRUE` will result in value of `1/mx` for the ##final age group, of `ax` and a value of 1 for the closeout of `qx`.
###' @return `cc` numeric vector with one element per age group, and which sums to ##the total difference in life expectancy between population 1 and 2.
###' @importFrom data.table shift
###' @export
###' @references
###' \insertRef{arriaga1989changing}{coddecomp}
###' \insertRef{preston2000demography}{coddecomp}
###' @examples
###' a <- .001
###' b <- .07
###' x <- 0:100
###' mx1 <- a * exp(x * b)
###' mx2 <- a/2 * exp(x * b)
###' mx1_causes <- matrix(c(mx1/4, mx1/4, mx1/2), length(x), 3) 
###' mx2_causes <- matrix(c(mx2/3, mx2/3, mx2/3), length(x), 3)
###' 
###' cod1 <- arriaga_cof(mx1 = mx1, mx2 = mx2, age = x
###'             mx1_causes = mx1_causes, 
###'             mx2_causes = mx2_Causes)
###' cod1_inst <- arriaga_cof(mx1 = mx1, mx2 = mx2, age = x
###'             mx1_causes = mx1_causes, 
###'             mx2_causes = mx2_Causes)
###' e01 <- mx_to_e0(mx1, age = x)
###' e02 <- mx_to_e0(mx2, age = x)
###' 
###' (delta <- e02 - e01)
###' sum(cod1); sum(cod1_inst)
###' 
###' #There are some differences between the Arriaga-Cause-of-death LE ##decomposition
###' #and our Arriaga-instantaneous proposal but they appear in the 3 decimal. 
###' 
###' #We shoudl include a function to plot the decomposition methods through the ##causes.
###' 
###Esta seria la descomposicion de Arriaga por causas de muerte
##arriaga_instantaneous_cof <- function(mx1,
##                                      mx2,
##                                      mx1_causes,
##                                      mx2_causes,
##                                      age = 0:(length(mx1) - 1),
##                                      sex1 = 't',
##                                      sex2 = sex1,
##                                      closeout = TRUE){
##  mx_mean  <- (mx1 + mx2) / 2
##  s_mean <- sen_arriaga_instantaneous(mx_mean, age = age)
##  
##  mx_causes <- mx2_causes - mx1_causes
##  
##  cc <- s_mean*mx_causes
##  
##  return(cc)
##}