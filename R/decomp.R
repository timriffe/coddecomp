
# methods list:
#' @title decompose life expectancy using a variety of methods
#' @description A variety of exact or sympototically exact life expectancy decomposition methods are implemented. These include the lifetable response experiment, using three potential sensitivity functions (lifetable-based `ltre_lt`, an instantaneous Arriaga-derived sensitivity `ltre_arriaga_instantaneous`, or numerical sensitivity `ltre_numerical`), the directional Arriaga method `arriaga`, a symmetrical Arriaga method `arriaga_symmetrical`, or the linear intregal (a.k.a pseudo continuous) method of Horiuchi et al `horiuchi`. These all give similar results, except for directional `arriaga`, which is the most different.
#' @param mx1 numeric. age-structured mortality rates for population 1
#' @param mx2 numeric. age-structured mortality rates for population 2, must be same length as `mx1`
#' @param age integer. lower bound of each age group, must be same length as `mx1`
#' @param sex1 character. `"m"`,`"f"`, or `"t"`, affects a0 treatment.
#' @param sex2 character. `"m"`,`"f"`, or `"t"`, affects a0 treatment.
#' @param method character. `"ltre_lt"`,`"ltre_arriaga_instantaneous"`,`"ltre_numerical"`,`"arriaga"`, `"arriaga_symmetrical"`, or `"horiuchi"`
#' @param N integer. used for LTRE or Horiuchi methods
#' @param closeout logical. Do we handle closeout, or truncate at top age_
#' @param ... optional arguments passed to `numDeriv::grad()`
#' @export
# needs examples
decomp_e0 <- function(mx1, 
                       mx2, 
                       age, 
                       sex1, 
                       sex2, 
                       method = c("ltre_lt","ltre_arriaga_instantaneous","ltre_numerical","arriaga", "arriaga_symmetrical", "horiuchi"), 
                       N = 20,
                       closeout = TRUE,
                       ...){
   method <- tolower(method)
   method <- match.arg(method = method,
                       choices = c("ltre_lt","ltre_arriaga_instantaneous", "ltre_numerical","arriaga","arriaga", "arriaga_symmetrical", "horiuchi"))
   
   # if sex1 != sex2, we re-call the function twice and take the average...
   if ((sex1 != sex2) & !method %in% c("arriaga","arriaga_symmetrical")){
     .sex1 <- sex1
     .sex2 <- sex2
     cc1 <- decomp_e0(mx1 = mx1, 
                      mx2 = mx2, 
                      age = age, 
                      sex1 = .sex1, 
                      sex2 = .sex1, 
                      method = method,
                      closeout = closeout,
                      N=N,
                      ...)
     cc2 <- decomp_e0(mx1 = mx1, 
                      mx2 = mx2, 
                      age = age, 
                      sex1 = .sex2, 
                      sex2 = .sex2, 
                      method = method,
                      closeout = closeout,
                      N=N,
                      ...)
     cc <- (cc1 + cc2) / 2
   }
   
   # method based on lifetable sensitivity
   if (method == "ltre_lt"){

       cc <- ltre(func = mx_to_e0,
                   pars1 = mx1,
                   pars2 = mx2,
                   dfunc = sen_e0_mx_lt,
                   sex = sex1,
                   age = age,
                   closeout = closeout,
                   N = N)
   }
   
   # method based on lifetable sensitivity
   if (method == "ltre_arriaga_instantaneous"){
     cc <- ltre(func = mx_to_e0,
                pars1 = mx1,
                pars2 = mx2,
                dfunc = sen_arriaga_instantaneous,
                sex = sex1,
                age = age,
                closeout = closeout,
                N = N)
   }
   
   if (method == "ltre_numerical"){
     cc <- ltre(func = mx_to_e0,
                pars1 = mx1,
                pars2 = mx2,
                dfunc = sen_num,
                sex = sex1,
                age = age,
                closeout = closeout,
                N = N)
   }
   
   if (method == "arriaga"){
     cc <- arriaga(mx1 = mx1,
                   mx2 = mx2,
                   sex1 = sex1,
                   sex2 = sex2,
                   age = age,
                   closeout = closeout)
   }
   
   if (method == "arriaga_symmetrical"){
     cc1 <- arriaga(mx1 = mx1,
                   mx2 = mx2,
                   sex1 = sex1,
                   sex2 = sex2,
                   age = age,
                   closeout = closeout)
     cc2 <- arriaga(mx1 = mx2,
                    mx2 = mx1,
                    sex1 = sex2,
                    sex2 = sex1,
                    age = age,
                    closeout = closeout)
     cc <- (cc1 - cc2) / 2
   }
   
   if (method == "horiuchi"){
     cc <- horiuchi(func = mx_to_e0,
                    pars1 = mx1,
                    pars2 = mx2,
                    N = N,
                    sex = sex1,
                    age = age,
                    closeout = closeout)
   }
   
   return(cc)
 }















