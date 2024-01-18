# testing Arriaga with composition.
# A little script I wrote for Michael Lachanski. And the world

# simiulate 2 groups at two time points,
# Gompertz mort, with parameters a and b

# we will use functions from a package in development coddecomp
# install like so:
# remotes::install_github("timriffe/coddecomp")
library(coddecomp)

gomp_mx <- function(a,b,x){
  a * exp(b * x)
}
a_g1t1 <- .001
a_g2t1 <- .0015

a_g1t2 <- .0009
a_g2t2 <- .0017

b_g1t1 <- .069
b_g2t1 <- .071

b_g1t2 <- .068
b_g2t2 <- .072

x <- 0:100

# calculate rates
mx_g1t1 <- gomp_mx(a_g1t1,b_g1t1,x)
mx_g2t1 <- gomp_mx(a_g2t1,b_g2t1,x)
mx_g1t2 <- gomp_mx(a_g1t2,b_g1t2,x)
mx_g2t2 <- gomp_mx(a_g2t2,b_g1t2,x)

# calculate e0 for groups
e0g1t1 <- mx_to_e0(mx_g1t1, age = x)
e0g2t1 <- mx_to_e0(mx_g2t1, age = x)
e0g1t2 <- mx_to_e0(mx_g1t2, age = x)
e0g2t2 <- mx_to_e0(mx_g2t2, age = x)

# changing initial composition
pi1 <- c(.4,.6)
pi2 <- c(.7,.3)

# weight for aggregate e0
(e01 <- e0g1t1 * pi1[1] + e0g2t1 * pi1[2])
(e02 <- e0g1t2 * pi2[1] + e0g2t2 * pi2[2])

# check values
e01;e0g1t1;e0g2t1
e02;e0g1t2;e0g2t2

# Note, Arriaga is directional, in the sense that it matters whether you decompose t1 w respect to t2 or vice versa. Not only the sign changes, but the magnitudes change. In this case, differences are not extreme, but still we take the average. For temporal change there is a case to be made for prefering "forward" decomposition, but for group comparisons we might prefer the average (or symmetrical version). Here I average them.
g1_arriaga_1 <- arriaga(mx1 = mx_g1t1,mx2 = mx_g1t2, age = x, closeout = FALSE)
g1_arriaga_2 <- -arriaga(mx1 = mx_g1t2,mx2 = mx_g1t1, age = x, closeout = FALSE)
g1_arriaga <- (g1_arriaga_1 + g1_arriaga_2) / 2

g2_arriaga_1 <- arriaga(mx1 = mx_g2t1, mx2 = mx_g2t2, age = x, closeout = FALSE)
g2_arriaga_2 <- -arriaga(mx1 = mx_g2t2, mx2 = mx_g2t1, age = x, closeout = FALSE)
g2_arriaga <- (g2_arriaga_1 + g2_arriaga_2) / 2

# take a look at within-group Arriaga decompositions.
plot(x, g1_arriaga, type = 'l', ylim = c(-.02,.05))
lines(x, g2_arriaga, col = "red")
abline(h=0)

# within-group changes, check sums
(g1_d <- sum(g1_arriaga) )
(e0g1_d  <- e0g1t2 - e0g1t1)

(g2_d <- sum(g2_arriaga))
(e0g2_d  <- e0g2t2 - e0g2t1)

# total change
(e0_chg  <- e02 - e01)

# the weighting scheme is a Kitagawa situation:
sum(pi1 * c(e0g1t1,e0g2t1))

# ergo:
avge0 <- (c(e0g1t1,e0g2t1) + c(e0g1t2,e0g2t2)) / 2
chge0 <- (c(e0g1t2,e0g2t2) - c(e0g1t1,e0g2t1))

# compositional effect
# we sum because the individual values are meaningless.
(comp_part <- sum((pi2 - pi1) * (avge0))) 
# for this component the values are meaningful
(within_part  <- (chge0) * ((pi1+pi2)/2 ))

# check sums:
comp_part + sum(within_part)
e0_chg

# so, the effect of compositional change is:
comp_part

# the age-specific effects of rate changes in group 1 are:
plot(x, g1_arriaga)

# the age-specific effects of rate changes in group 2 are:
plot(x, g2_arriaga)

# note that:
e0_chg 

# does not equal
sum(g1_arriaga) + sum(g2_arriaga) + comp_part

# For interpreting overall leverage, do we simply rescale?
# rescaling is tricky because a vector can have both negatives and
# positives. But I'll try?

g1_arriaga_for_total <- g1_arriaga * (within_part[1] / sum(g1_arriaga))
g2_arriaga_for_total <- g2_arriaga * (within_part[2] / sum(g2_arriaga))

sum(g1_arriaga_for_total) + sum(g2_arriaga_for_total) + comp_part
e0_chg 

plot(x, g1_arriaga_for_total, ylim = c(-.01,.025), type='l',col="blue",main = "proposed Arriaga decomposition with composition component\nThis is the non-compositional part:")
lines(x,g2_arriaga_for_total, col = "red")
abline(h=0)

# and now a Horiuchi comparison, to check consistency.

e0_composed_vec <- function(vec){
  mx1 <- vec[1:101]
  mx2 <- vec[102:202]
  pii <- vec[203:204]
  
  e01 <- mx_to_e0(mx1, age = 0:100, closeout = FALSE)
  e02 <- mx_to_e0(mx2, age = 0:100, closeout = FALSE)
  
  sum(c(e01,e02) * pii)
}
vec1 <- c(mx_g1t1, mx_g2t1, pi1)
vec2 <- c(mx_g1t2, mx_g2t2, pi2)

e0_composed_vec(vec1)
e01
e0_composed_vec(vec2)
e02

dh <- DemoDecomp::horiuchi(e0_composed_vec, vec1, vec2, N = 20)

# check sum
sum(dh)
e02-e01

g1_horiuchi_for_total <- dh[1:101]
g2_horiuchi_for_total <- dh[102:202]

# compare composition attribution
(comp_part_h <- sum(dh[203:204]))
comp_part

# compare g1 sums:
sum(g1_horiuchi_for_total)
sum(g1_arriaga_for_total)

# compare g2 sums:
sum(g2_horiuchi_for_total)
sum(g2_arriaga_for_total)

# compare age patterns
plot(g2_horiuchi_for_total, type = 'l')
lines(g2_arriaga_for_total, col = "red", lty = 2)

plot(g1_horiuchi_for_total, type = 'l')
lines(g1_arriaga_for_total, col = "red", lty = 2)
