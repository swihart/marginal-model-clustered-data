
## code from the OP Question: edit `data` to `d` 
require(geepack)
d = read.csv(url("http://folk.uio.no/mariujon/data.csv"))
fit = geeglm(moden ~  power, id = defacto, data=d, corstr = "exchangeable", family=binomial)
summary(fit)
# plot(moden ~ power, data=d)
# x = 0:2500
# y = predict(fit, newdata=data.frame(power = x), type="response" )
# lines(x,y)


#######################################
## gee - arguments / tricks.  Nothing's working...
#######################################
- glm(moden ~  power, family = binomial, data=d)$deviance / 2

strt <- coef(glm(moden ~  power, family = binomial, data=d))
strt
## error
geeglm(moden ~  power, id = defacto, data=d, corstr = "exchangeable", family=binomial,
       start=strt)
## same as geeglm() in OP
geese(moden ~ power, id = defacto, data=d, corstr = "exchangeable", family=binomial,
      b=strt)

yvec <- as.numeric(c(d$moden))
Xmat <- matrix(NA, nrow=length(d$power), ncol=2)
Xmat[,1] <- 1
Xmat[,2] <- as.numeric(c(d$power))
idvec <- as.numeric(c(d$defacto))
## error:
geese.fit(y=yvec, 
          x=Xmat,
          id = idvec,
          family=binomial(),
          b=strt)


#######################################
## marginalized random intercept models
#######################################

## I'm so sorry but these methods use attach()
attach(d)

## Heagerty (1999) 
# marginally specifies a logit link and has a nonlinear conditional model
# the following code will not run if lnMLE is not successfully installed.  
# See https://faculty.washington.edu/heagerty/Software/LDA/MLV/
library(lnMLE)
L_N <- logit.normal.mle(meanmodel = moden ~ power,
                        logSigma= ~1,
                        id=defacto,
                        model="marginal",
                        data=d,
                        beta=strt,
                        r=10 ## see help, r:  Number of Gauss-Hermite quadrature points. 
                             ## The user may choose r=3, 5, 10, 20, or 50. The default value is r=20. 
                             ## above 10 really slowed for your dataset
                        ) 
print.logit.normal.mle(L_N)
## variance of normal distribution
##exp(L_N$alphas[1])^2
## intraclass correlation
##exp(L_N$alphas[1])^2 / (pi^2/3 + exp(L_N$alphas[1])^2)

library("stabledist")
library("gnlm")
library("repeated")
## "gnlmix4MMM.R" is a modified gnlmix() function.  The following are the edits to
## Lindsey's original gnlmix() in order to get a function
## that worked for MMMs where the conditional model is a random intercept model:
## To make gnlmix4MMM.R, edit gnlmix.R of repeated package in two places:
##  np <- npl + nps + 1         becomes   np <- npl
##   p <- c(pmu, pshape, pmix)  becomes    p <- pmu
## These are the only edits required to estimate a MMM, which puts the variance
## component into the non-linear predictor of the mean.  Warnings/errors still might
## be generated about Mixing Dispersion Parameters and Correlations because of these
## modifications.  These are superficial:  The Mixing Dispersion parameter is now
## listed in the location parameters and the correlations can be obtained by
## directly accessing fitted.obj$corr; the error relates to assigning dimnames.  See examples below.
source("gnlmix4MMM.R")

## need outcomes to be in two column format
y <- cbind(d$moden,(1-d$moden))

## Wang and Louis 2003
## Wang and Louis (2003) logit-logit-bridge model...in R!  Their paper provided SAS code.
## the "mu" specification is especially busy, but a work around was not easily found.
## A careful look shows structure:  mu = 1/(1+exp(-Delta - random_intercept)
## 1/(1+exp(-eta)) is the conditional inverse logit-link (logistic CDF)
## eta = Delta + random_intercept
## where Delta can be simplified to a scaling of marginal parameters a0 and a1
## because marginalizing the conditional logistic model with Bridge-distributed random intercept
## yields a logistic distribution with inflated variance; since logistic distribution is closed under
## scaling, it can be scaled back to the standard logistic to yield the marginal link.  
## which allows conditional specification:              a0_c + a1_c*black
## to be represented in terms of marginal coefficients:(a0   + a1  *black)*sqrt(1+3/pi/pi*exp(pmix))
## and a0, a1 are the marginal parameters.  pmix is the variance on the log scale.
## where random_intercept is the probabilty integral transform of rand ~ N(0, exp(pmix))
## to give an entity with the Bridge distribution.  Using pnorm(rand/sqrt(exp(pmix))) gives
## a Uniform(0,1) that is non-linear tranformed by the inverse CDF of the Bridge distibution.

start.time <- Sys.time()
LLB  <- gnlmix4MMM(   y = y,
                      distribution = "binomial",
                      mixture = "normal",
                      random = "rand",
                      nest = defacto,
                      mu = ~ 1/(1+exp(-(a0 + a1*power)*sqrt(1+3/pi/pi*exp(pmix)) - sqrt(1+3/pi/pi*exp(pmix))*log(sin(pi*pnorm(rand/sqrt(exp(pmix)))/sqrt(1+3/pi/pi*exp(pmix)))/sin(pi*(1-pnorm(rand/sqrt(exp(pmix))))/sqrt(1+3/pi/pi*exp(pmix)))))),
                      pmu = c(strt, log(1)),
                      pmix = log(1),
                      gradtol = 1e-5,
                      steptol = 1e-5,
                      iterlim = 1e5,
                      eps = 1e-4,
                      points = 6)
end.time <- Sys.time()
time.taken <- end.time - start.time
time.taken

print("call")
LLB$call
print("likelihood")
LLB$maxlike
print("iterations")
LLB$iterations
print("code: 1 -best 2-ok 3,4,5 - problem")
LLB$code
print("coefficients")
LLB$coeff
print("se")
LLB$se
print("corr")
LLB$corr
## print everything, even including some slots that are empty due to
## edits performed on gnlmix() to render gnlmix4MMM()
print(LLB)
## variance of bridge distribution
##exp(LLB$coeff[3])
## intraclass correlation
##exp(LLB$coeff[3]) / (pi^2/3 + exp(LLB$coeff[3]))





## Caffo and Griswold (2006) 
## marginal link-probit-normal model, with a marginal logit-link.
## Variance pmix is specified on log scale; random intercept rand is an entity to be
## specified in the mean mu.
## mu is complicated, but working from the outside in, left to right:
## pnorm() is the conditional inverse link (aka conditional CDF)
## The resulting integral of marginalizing a conditional probit model with normally distributed
## random intercept is a normal distribution with additive means and additive variances.
## Therefore the Delta function can be written as the inverse normal with variance 1+tau^2 applied
## to the marginal inverse link (aka marginal distribution, logistic: 1/(1+exp())).
## Due to Normal distributions being closed under scaling, we can equivalently write
## Delta as a standard Normal inverse (qnorm()) and scale the operators, which we do with tau^2 = exp(pmix).
start.time <- Sys.time()
LPN  <- gnlmix4MMM(   y = y,
                      distribution = "binomial",
                      mixture = "normal",
                      random = "rand",
                      nest = defacto,
                      mu = ~pnorm(qnorm(1/(1+exp(-a0 - a1*power)))*sqrt(1+exp(pmix)) + rand),
                      pmu = c(strt, log(1)),
                      pmix = log(1),
                      gradtol = 1e-5,
                      steptol = 1e-5,
                      iterlim = 1e5,
                      eps = 1e-4,
                      points = 6)

print("call")
LPN$call
print("likelihood")
LPN$maxlike
print("iterations")
LPN$iterations
print("code: 1 -best 2-ok 3,4,5 - problem")
LPN$code
print("coefficients")
LPN$coeff
print("se")
LPN$se
print("corr")
LPN$corr
## print everything, even including some slots that are empty due to
## edits performed on gnlmix() to render gnlmix4MMM()
print(LPN)
## variance of normal distribution
exp(LPN$coeff[3])
## intraclass correlation
exp(LPN$coeff[3]) / (1 + exp(LPN$coeff[3]))

end.time <- Sys.time()
time.taken <- end.time - start.time
time.taken




## The SSS method, currently not behaving well :-(
start.time <- Sys.time()
alp <- 1.89
gam <- 1.2
SSS  <- gnlmix4MMM(   y = y,
                      distribution = "binomial",
                      mixture = "normal",
                      random = "rand",
                      nest = defacto,
                      mu = ~pstable( (a0 + a1*power)*(gam/(gam^alp+ (exp(pmix))^alp)^(1/alp)  )^(-1) + qstable( min(max(pnorm(rand/sqrt(exp(pmix))), 1e-200), 1-1e-16), alp,0,exp(pmix),0,0), alp,0,gam,0,0),
                      pmu = c(strt, log(1)),
                      pmix = log(1),
                      gradtol = 1e-5,
                      steptol = 1e-5,
                      iterlim = 1e5,
                      eps = 1e-4,
                      points = 6)
end.time <- Sys.time()
time.taken <- end.time - start.time
time.taken


print("call")
SSS$call
print("likelihood")
SSS$maxlike
print("iterations")
SSS$iterations
print("code: 1 -best 2-ok 3,4,5 - problem")
SSS$code
print("coefficients")
SSS$coeff
print("se")
SSS$se
print("corr")
SSS$corr
## print everything, even including some slots that are empty due to
## edits performed on gnlmix() to render gnlmix4MMM()
print(SSS)

end.time <- Sys.time()
time.taken <- end.time - start.time
time.taken

rbind("L_N"=L_N$beta, "LLB" = LLB$coefficients[1:2], "LPN"=LPN$coefficients[1:2])

rbind("L_N"=L_N$logL, "LLB" = -LLB$maxlike, "LPN"=-LPN$maxlike)
