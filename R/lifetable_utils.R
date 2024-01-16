# helper for Tx
rcumsum <- function(x){
  rev(cumsum(rev(x)))
}

# assumes step function for mx
# mx_to_lx <- function(mx){
#   mx[is.na(mx)] <- 0
#   lx <- exp(-cumsum(mx))
#   lx <- c(1,lx)
#   lx[1:length(mx)]
# }
#' @importFrom DemoTools lt_rule_ak_m0_a0

mx_to_ax <- function(mx, age = 0:(length(mx1)-1), sex = "t", closeout = TRUE){
  stopifnot(all(diff(age) == 1))
  ax <- rep(.5, length(mx))
  sex <- sex |> tolower() |> substr(1,1)
  sex <- if_else(sex == "b","t",sex)
  stopifnot(sex %in% c("m","f","t"))
  if (min(age) == 0){
    if (sex != "t"){
      a0 <- lt_rule_ak_m0_a0(M0 = mx[1], Sex = sex)
    } else {
      a0m <- lt_rule_ak_m0_a0(M0 = mx[1], Sex = "m")
      a0f <- lt_rule_ak_m0_a0(M0 = mx[1], Sex = "f")
      a0  <- (a0m + a0f) / 2
    }
    ax[1] <- a0
  }
  if (closeout){
    ax[length(ax)] <- 1 / mx[length(mx)]
  }
  ax
}

mx_to_qx <- function(mx, ax, age = 0:(length(mx1)-1), closeout = TRUE){
  lt_id_ma_q(nMx = mx, 
             nax = ax, 
             AgeInt = rep(1, length(mx)), 
             closeout = closeout)
}

qx_to_lx <- function(qx){
  lt_id_q_l(nqx = qx, radix = 1)
}

lx_to_dx <- function(lx){
  lt_id_l_d(lx)
}


lx_to_dx <- function(lx){
  -diff(c(lx,0))
}


ald_to_Lx <- function(ax,lx,dx){
  lt_id_lda_L(nax = ax,
              lx = lx,
              ndx = dx,
              AgeInt = rep(1,length(ax)))
}

lL_to_ex <- function(lx, Lx){
  Tx <- rcumsum(Lx) 
  ex <- Tx / lx
  ex
}

mx_to_ex <- function(mx, age, sex = 't'){
  ax <- mx_to_ax(mx = mx, 
                 age = age, 
                 sex = sex)
  qx <- mx_to_qx(mx, ax)
  lx <- qx_to_lx(qx)
  dx <- lx_to_dx(lx)
  Lx <- ald_to_Lx(ax = ax,
                  lx = lx,
                  dx = dx)
  ex <- lL_to_ex(lx = lx, Lx = Lx)
  ex
}

mx_to_e0 <- function(mx, age, sex = 't'){
  ex <- mx_to_ex(mx = mx, age = age, sex = sex) 
  ex[1]
}
