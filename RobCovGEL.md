Covariate adjustment using GEL and robust estimators
================
Mara Delesa-Velina
2024-01-29

-   [Introduction](#introduction)
-   [Data example](#data-example)
-   [Simulation experiment](#simulation-experiment)
-   [References](#references)

## Introduction

This is the accompanying R code for my paper “Covariate adjustment for
robust estimators” (preprint). In this paper I present a methodology to
carry out covariate adjustment based on generalised empirical likelihood
and robust estimators, in particular, the trimmed means and smoothed
Huber estimator.

The implementation of generalized empirical likelihood method is based
on R package `gmm` (Chausse, 2010). Smoothed Huber estimator is computed
using the package `smoothmest`.

``` r
# Load the libraries
library(WRS2)
library(ggpubr)
library(tidyverse)
library(smoothmest)
library(gmm)
library(lmtest)
library(mvtnorm)
```

For the trimmed means methods, a custom code is created and available in
this repository:

``` r
source("trimmed_GEL_v1.R")
```

## Data example

We demonstrate our methods on the data sample `WRS2::electric`. Data
contains reading test scores for 192 school classes of four different
grades. A half of classes of each grade were randomized into the
treatment group, that was watching an educational TV show “The electric
company”, and control group that did not. At the beginning and at the
end of the school year, students in all the classes were given a reading
test, yielding the pre-test and post-test scores. The average test
scores per class were recorded. We are interested in possible
differences between treatment and control group post-test scores with
the pre-test score being a covariate.

The data is presented in Figure 1. We can observe that the classical
ANCOVA assumptions do not hold. First, data is not likely to be normally
distributed, it is left skewed and containing outliers. We can observe
within group heterogeneity since the dispersion of regression errors
decrease with the increase of the pre-test score. In addition,
non-parallel regression lines indicate an interaction effect between the
group and the covariate.

![Visualization of pupils’ test scores for electric dataset. On the
left, continuous and dashed lines represent LOESS and linear smoother,
respectively, for each of the control and treatment
groups.](RobCovGEL_files/figure-gfm/unnamed-chunk-4-1.png)

### Results based on linear models

We first compare post-test scores based on four linear models: (1)
ANOVA, (2) ANCOVA with pre-test score as covariate, (3) ANCOVA full
model with pre-test score as a covariate and covariate-group
interaction, (4) ANCOVA full model as above but robust confidence
intervals are calculated based on heteroscedasticity-consistent
covariance matrix sandwich estimate with `vcovHC()`.

The results are presented below and correspond to manuscript section 4
Data analysis.

``` r
fmt <- function(v) {
  #Function for printing confidence intervals in a nice format
  #v - vector containing CI lower and upper values 
  paste0("(", paste(round(v, 3), collapse = ", "), ")")
}
###Normal ANOVA
m0 <- lm(Posttest ~ Group, electric)
cf <- coeftest(m0)
ci <- confint(m0)
rez1 <- data.frame(Method = "ANOVA", 
                   delta = cf["Grouptreatment", 1] %>%  round(3), 
                   delta_ci = fmt(ci["Grouptreatment", ]), 
                   mu = cf["(Intercept)", 1] %>%  round(3), 
                   mu_ci = fmt(ci["(Intercept)", ]),
                   b.x = NA, b.x_ci = NA, b.xy = NA, b.xy_ci = NA)

###ANCOVA
m1 <- lm(Posttest ~ Pretest + Group, electric)
cf <- coeftest(m1)
ci <- confint(m1)
rez2 <- data.frame(Method = "ANCOVA", 
                   delta = cf["Grouptreatment", 1] %>%  round(3), 
                   delta_ci = fmt(ci["Grouptreatment", ]), 
                   mu = cf["(Intercept)", 1] %>%  round(3), 
                   mu_ci = fmt(ci["(Intercept)", ]),
                   b.x = cf["Pretest", 1],
                   b.x_ci = fmt(ci["Pretest", ]),
                   b.xy = NA, b.xy_ci = NA)

###ANCOVA w. interactions
m2 <- (lm(Posttest ~ Pretest * Group, electric))
cf <- coeftest(m2)
ci <- confint(m2)
rez3 <- data.frame(Method = "ANCOVA Inter. ", 
                   delta = cf["Grouptreatment", 1] %>%  round(3), 
                   delta_ci = fmt(ci["Grouptreatment", ]), 
                   mu = cf["(Intercept)", 1] %>%  round(3), 
                   mu_ci = fmt(ci["(Intercept)", ]),
                   b.x = cf["Pretest", 1] %>%  round(3),
                   b.x_ci = fmt(ci["Pretest", ]),
                   b.xy = cf["Pretest:Grouptreatment", 1] %>%  round(3), 
                   b.xy_ci = fmt(ci["Pretest:Grouptreatment", ]))

###ANCOVA w. interactions w.robust SE
cf <- coeftest(m2, vcov = vcovHC(m2))
ci <- confint(cf)
rez4 <- data.frame(Method = "ANCOVA Inter. Rob.", 
                   delta = cf["Grouptreatment", 1] %>%  round(3), 
                   delta_ci = fmt(ci["Grouptreatment", ]), 
                   mu = cf["(Intercept)", 1] %>%  round(3), 
                   mu_ci = fmt(ci["(Intercept)", ]),
                   b.x = cf["Pretest", 1] %>%  round(3),
                   b.x_ci = fmt(ci["Pretest", ]),
                   b.xy = cf["Pretest:Grouptreatment", 1] %>%  round(3), 
                   b.xy_ci = fmt(ci["Pretest:Grouptreatment", ]))
### All ANOVAs
rez <- rbind(rez1, rez2, rez3, rez4)
```

|          |                  |                  |                  |                    |
|:---------|:-----------------|:-----------------|:-----------------|:-------------------|
| Method   | ANOVA            | ANCOVA           | ANCOVA Inter.    | ANCOVA Inter. Rob. |
| delta    | 5.657            | 4.734            | 10.151           | 10.151             |
| delta_ci | (0.653, 10.662)  | (2.445, 7.022)   | (4.812, 15.49)   | (1.024, 19.278)    |
| mu       | 94.321           | 61.558           | 58.890           | 58.890             |
| mu_ci    | (90.782, 97.859) | (58.657, 64.459) | (55.16, 62.62)   | (53.291, 64.489)   |
| b.x      |                  | 0.4600189        | 0.4970000        | 0.4970000          |
| b.x_ci   |                  | (0.426, 0.494)   | (0.45, 0.545)    | (0.437, 0.558)     |
| b.xy     |                  |                  | -0.075           | -0.075             |
| b.xy_ci  |                  |                  | (-0.142, -0.008) | (-0.173, 0.023)    |

Test results for eletric data based on linear models

### GEL for the means

For GEL based estimation, the main R command is `gel` that carries out
the GEL optimization. In addition, we need to set up the moment
estimating conditions in a separate R function and provide it to `gel`
routine. This function always takes two arguments: a parameter vector
`tet` and a data `W` in matrix format. Please see the article (Chausse,
2016) for more examples and details.

First consider estimating mean treatment effect $\Delta$ using GEL and
adjusting for one covariate. We have a data matrix $W=(Y, X, Z)$ and we
have in total three parameters $\theta=(\Delta, \mu_1, \mu_X)$, where
$\mu_1$ is the mean outcome for the treatment group, and $\mu_X$ is the
common mean of the covariate in both groups. Our problem is determined
by four moment conditions

$$ E(g(W, \mu_1, \Delta, \mu_{X})) = E\begin{pmatrix} 
              (Y - \mu_1 - \Delta Z) \\
              Z  (Y - \mu_1 - \Delta Z)\\
              (X_1 -\mu_{X})\\
              Z  (X_1 -\mu_{X})
           \end{pmatrix}
           =0. $$

Note that this is the setup discuss in (Chausse, 2016). We define the
moment conditions for mean treatment effect in function `g.means`.

``` r
#GEL for mean treatment effect
g.means <- function(tet, W) {
  # tet = c(delta, mu, mu.x) 
  # W = (Y, X, Z), where Y-outcome, X-covariate, Z-group indicator
  y <- W[, 1]
  z <- W[, 2]
  x <- W[, 3]
  #Moment conditions:
  m1 <- x - tet[3]
  m2 <- z * (x - tet[3])
  m3 <- (y - tet[2] - tet[1] * z)
  m4 <- z * (y - tet[2] - tet[1] * z)
  f <- cbind(m1, m2, m3, m4)
  return(f)
}
```

Next, we feed our moment estimating function `g.means` to the R command
`gel` that carries out the optimization of the GEL function. In
addition, we need to provide starting values of the $\theta$ parameter,
which we obtain by the corresponding sample means. We also specify the
particular GEL family member we want to use under the `type` argument.
Here and forward we provide an example for the empirical likelihood
estimate, `type="EL"`. Note that other methods used in our paper were
exponential tilting “ET” and continuously updating estimator “CUE”. For
convenience we define vectors $y$, $x$ and $z$ that represent the
post-test score, pre-test score and treatment indicator (treatment group
coded as 1), respectively.

``` r
### Define our data vectors
y <- electric$Posttest
x <- electric$Pretest
z <- as.numeric(electric$Group == "treatment")
W <- cbind(y, z, x)
### Starting values for theta (sample means)
tet.init <- c(mean(y[z == 1]) - mean(y[z == 0]), mean(y[z == 0]), 
              mean(x))
###  Carry out GEL optimization
gel.1 <- gel(g.means, x = cbind(y, z, x), 
             type = "EL", 
             tet = tet.init)
```

Finally, we obtain the confidence intervals for the parameter estimates
using the inverse likelihood ratio method and print out the results.

``` r
### Get Confidence intervals
cf1 <- confint(gel.1, type = "invLR")$test
### Print parameter estimates and confidence intervals
cbind(gel.1$coefficients, cf1)
```

    ##                        0.025     0.975
    ## Theta[1]  4.734823  2.450506  7.091047
    ## Theta[2] 94.811091 91.908506 97.444102
    ## Theta[3] 72.204921 67.266274 76.834103

### GEL for trimmed means

The idea is to trim the data sample with respect to the outcome values
in each of the treatment groups. We then proceed to estimate the trimmed
mean treatment effect, for which we use the same moment conditions
function `gel.means` as above, but on the trimmed sample. Thus we have

$$E(g(W, \mu_1, \Delta, \mu_{X}^{TM})) = E\begin{pmatrix}
              (Y - \mu_1 - \Delta Z) \\
              Z  (Y - \mu_1 - \Delta Z)\\
              (T -\mu^{TM}_{X})\\
              Z  (T -\mu^{TM}_{X})
           \end{pmatrix}
           =0, $$ where $T$ is the random variable having the truncated
distribution with respect to the covariate $X$, and $\mu^{TM}_{X}$
represent the trimmed mean of the covariate $X$. However, the trimming
modifies the limiting distribution of the GEL statistic and has to be
corrected for by a scaling factor.

We first choose the trimming level (here we use 10% trimmed means) and
trim the data. Then we calculate a correction for the confidence level
of trimmed means inverse likelihood ratio test that corresponds to the
scaling factor.

``` r
### Sort the dataset
# Set the trimming level
a <- 0.1
#Sort the dataset wrt to outcome in each treatment group
Ws <- W %>%
  as.data.frame() %>% 
  arrange(y) %>% 
  group_by(z)  %>%
  mutate(keep = trim.index(y, a, a)$keep) 

### Trim the sample
Wt <- Ws %>% 
  filter(keep == T)
yt <- Wt$y
xt <- Wt$x
zt <- Wt$z

### Calculate confidence level for nonscaled chi-square
conflev.TM <-  tmeans.conf.level(Ws$y[Ws$z == 0], Ws$y[Ws$z == 1], 
                                mu.t = 0, alpha = a, 
                                beta = a, conf.level = 0.95)
```

Next, we set the starting values of $\theta$ and optimize the GEL
parameters. Finally, we estimate the confidence interval under the
corrected confidence level. Note that we provide argument `fact=4.5` for
`confint.gel` that controls the span of search for the inversion of the
test in terms of standard deviations. The default values is 3, but we
observed that it might not be enough when data is not normally
distributed.

We do not print the results of the optimization, we report them in Table
2 below for all methods together.

``` r
### Theta starting values as sample trimmed means
tet.init <- c(mean(yt[zt == 1], trim = a) - mean(yt[zt == 0], trim = a),
            mean(yt[zt == 0], trim = a), 
            mean(xt, trim = a))
### Optimize GEL parameters
gel.2 <- gel(g.means, 
             x = cbind(yt, zt, xt),
             tet = tet.init)
### Calculate confidence interval
cf2 <- confint(gel.2, 
               level = conflev.TM,
               type = "invLR", fact = 4.5)$test
#Print results
#cbind(gel.2$coefficients, cf2)
```

### GEL for Huber estimator

We now want to adjust covariates not by their mean values, but by their
M-estimates instead. That is, we expect the smoothed M-estimate of the
covariate in the control and treatment groups being equal. We have again
four moment conditions, but the last two now involve the smoothed
M-estimator’s psi-function $\tilde\psi$, and $\mu^M_{X}$ indicated the
real value of the covariate’s M-estimator:

$$ E(g(W, \theta))=E(g(W, \mu_1, \Delta, \mu_{X}^M)) = E\begin{pmatrix}
              (Y - \mu_1 - \Delta Z) \\
              Z  (Y - \mu_1 - \Delta Z)\\
               \tilde \psi \left( \frac{X - \mu^M_{X}}{\sigma_X} \right)\\
            Z \tilde \psi \left( \frac{X - \mu^M_{X}}{\sigma_X} \right) \\
           \end{pmatrix}
           =0.
$$

In this example, we choose smoothed Huber estimator as an appropriate
M-estimator.

For the starting values of $\theta$ parameter we use the corresponding
sample estimates - sample means for the outcome and smoothed Huber
estimate for the covariate. Note that the smoothed Huber estimate
involves a scale parameter $\sigma_X$, which has to be estimated
beforehand. We estimate it using the mean absolute deviation estimator
(MAD).

``` r
g.hub <- function(tet, W) {
  # tet = c(delta, mu, mu.x) 
  # W = (Y, X, Z), where Y-outcome, X-covariate, Z-group indicator
  y <- W[, 1]
  z <- W[, 2]
  x <- W[, 3]
  #moment conditions. Note that scale estimate s is a global parameter
  m1 <- y - tet[2] - tet[1] * z
  m2 <- z * (y - tet[2] - tet[1] * z)
  m3 <- smpsi((x - tet[3]) / s) 
  m4 <- z * smpsi((x - tet[3]) / s)
  f <- cbind(m1, m2, m3, m4)
  return(f)
}

### Theta starting values
tet.init <- c(mean(y[z == 1]) -  mean(y[z == 0]),
              mean(y[z == 0]), 
              smhuber(x)$mu)
### Scale parameter estimate for the covariate
s <- mad(x)
```

Finally we optimize the GEL parameters, providing the `g.hub` moment
conditions function to the `gel` procedure.

``` r
### Optimize GEL parameters
gel.3 <- gel(g.hub, x = cbind(y, z, x),
               tet = tet.init, type = "EL")
### Calculate CI
cf3 <- confint(gel.3, type = "invLR")$test
###Print results
#cbind(gel.3$coefficients, cf3)
```

### Results based on GEL methods

Finally we print all the results below, reporting the parameter
estimates and their respective confidence intervals for all GEL methods.
They correspond to manuscript section 4 Data analysis.

``` r
#Results of GEL for mean treatment effect
rez.gel1 <- data.frame(Method = "GEL Means", 
                       delta = gel.1$coefficients[1] %>% round(3), 
                       delta_ci = fmt(cf1[1, ]), 
                       mu= gel.1$coefficients[2] %>% round(3), 
                       mu_ci = fmt(cf1[2, ]),
                       mu.x = gel.1$coefficients[3] %>% round(3), 
                       mu.x_ci = fmt(cf1[3, ]))
#Results of GEL for trimmed mean of treatment effect
rez.gel2 <- data.frame(Method = "GEL TM", 
                      delta = gel.2$coefficients[1] %>% round(3), 
                      delta_ci = fmt(cf2[1, ]), 
                      mu= gel.2$coefficients[2] %>% round(3), 
                      mu_ci = fmt(cf2[2, ]),
                      mu.x = gel.2$coefficients[3] %>% round(3), 
                      mu.x_ci = fmt(cf2[3, ]))

#Results for GEL for mean treatment effect adjusted by smoothed Huber estimator
rez.gel3 <- data.frame(Method = "GEL TM", 
                      delta = gel.3$coefficients[1] %>% round(3), 
                      delta_ci = fmt(cf3[1, ]), 
                      mu = gel.3$coefficients[2] %>% round(3), 
                      mu_ci = fmt(cf3[2, ]),
                      mu.x = gel.3$coefficients[3] %>% round(3), 
                      mu.x_ci = fmt(cf3[3, ]))
```

|          |                  |                  |                  |
|:---------|:-----------------|:-----------------|:-----------------|
| Method   | GEL Means        | GEL TM           | GEL TM           |
| delta    | 4.735            | 5.236            | 4.391            |
| delta_ci | (2.451, 7.091)   | (2.691, 8.145)   | (1.824, 6.879)   |
| mu       | 94.811           | 96.785           | 95.035           |
| mu_ci    | (91.909, 97.444) | (93.499, 99.681) | (92.215, 97.594) |
| mu.x     | 72.205           | 74.047           | 80.249           |
| mu.x_ci  | (67.266, 76.834) | (67.404, 79.906) | (75.095, 84.481) |

## Simulation experiment

We show an example of how a simulation experiment can be set up. In
particular, we generate data for Scenario I from manuscript Section 3
Simulation Study. Scenario I involves normally distributed data with
balanced group sample sizes. We demonstrate the use of the three methods
described in the previous section - GEL covariate analysis for (1) mean
treatment effect, (2) trimmed mean treatment effect and (3) for mean
treatment effect when adjusting for smoothed Huber estimator of the
covariate.

First we generate the dataset of sample size $n=200$, involving outcome
and one covariate $(y, x)$ distributed by multivariate normal
distribution. We set a very small Monte Carlo repetition number $B=20$
for demonstration purpose.

``` r
### Set sample and repetition size
set.seed(12345)
B <- 200
n <- 200
### Set covariance structure for y and x
rho <- 0.5
Sigma <- matrix(rep(0.5, 9), nrow = 3)
diag(Sigma) <- 1
### Define treatment and control group indicator z
# d : proportion of observations in group 1 
d <- 0.5
z <- c(rep(0, d * 200), rep(1, (1-d) * 200))
### Generate data
data.norm <- lapply(1:B,  function(i) rmvnorm(n, sigma = Sigma))
```

### GEL for means and trimmed means

Next, we set up the function carrying out the simulation experiment. The
function takes $B$ data replicates as input and estimates for 4
different GEL versions: EL (empirical likelihood), CUE (continupusly
updating estimator), ET (exponential tilting), HD (Hellinger distance)
and also LS (the least squares method) for comparison.

Note that sometimes GEL convergence can fail and no confidence interval
can be produced. To ensure the code does not crash, we define a function
`confint.poss` that returns NA in case of a CI estimation failure. The
number of crashed calculations is tracked in our simulation experiment
and available in output as `num.nas`.

``` r
#object - gel object
#type - GEL method type
#fact - The span of search for the inversion of the test
#level - confidence level
#parm = 1 - we need confidence interval for the 1st parameter (Delta)
confint.poss <- possibly(.f = function(object, type, fact, level=0.95) 
  confint.gel(object, type = type, fact = fact, 
              parm = 1, level=level)$test[1,], 
              otherwise = c(NA, NA))

### Simulation experiment 
sim_1cov <- function(n, z, data, alpha = 0.1, tet0 = c(0, 0, 0), 
conf.nominal = 0.95, fact = 4) {
  ### Input
  # n - sample size
  # z - treatment group indicator (common for all replicates)
  # data - list of B data replicates, where each element
  #    is a matrix with outcome Y in 1st column and covariate X in 2nd column
  # tet0 - theta value under H0 (for bias, RMSE and CI coverage)
  # fact - et factor for CI search span in gel.confint()
  
  ### Output
  # bias, variance, RMSE for each theta element
  # coverage, ci.length - performance of confidence intervals for each of theta elements
  # num.nas - for how many replicates gel.confint() failed to estimate CI
  
  ### Calculations
  # Define arrays to store simulated parameter estimates
  #   and their CIs for 5 different methods (LS, EL, ET, CUE, HD)
  #   Each array has 3 rows for delta, mu and mu_x, respectively
  tet5 <- tet4 <- tet3 <- tet2 <- tet1 <- matrix(0, length(data), 3)
  ci5 <- ci4 <-ci3 <- ci2 <- ci1 <- matrix(0, length(data), 2)
  # Cycle to estimate GEL parameters for each data replicate
  for(i in 1: length(data)) {
    yns <- data[[i]][, 1]
    xns <- data[[i]][, 2]
    i.o <- order(yns)
    ys <- yns[i.o]; zs <- z[i.o]; xs <- xns[i.o]
    ### S2. Trim the sample
    #group indices
    i.z0 <- which(zs==0)
    i.z1 <- which(zs==1)
    #Indices after trimming
    i.tm0 <- trim.index(ys[i.z0], alpha, alpha)$keep.ind
    i.tm1 <- trim.index(ys[i.z1], alpha, alpha)$keep.ind
    #Corrected conf.level for EL 
    conflev.TM <- tmeans.conf.level(ys[i.z0], ys[i.z1], 
                                    mu.t = 0, alpha=alpha, 
                                    beta = alpha, conf.level = conf.nominal)
    #Data  after trimming
    yt <- c(ys[i.z0][i.tm0], ys[i.z1][i.tm1])
    x1t <- c(xs[i.z0][i.tm0], xs[i.z1][i.tm1])
    zt <- c(rep(0, length(i.tm0)), rep(1, length(i.tm1)))
    tet.init <- c(mean(yt[zt == 1]) - mean(yt[zt == 0]),
                  mean(yt[zt == 0]), 
                  mean(x1t))
    # ANCOVA LM
    rez.lm <- lm(yt ~ zt + x1t) # 
    tet1[i, ] <- coefficients(rez.lm)[1:3]
    ci1[i, ] <- confint(rez.lm)[2, ] #
    #GEL, EL
    rez.gel <- gel(g.means, x = cbind(yt, zt, x1t),
                   tet = tet.init)
    tet2[i, ] <- rez.gel$coefficients
    ci2[i, ] <- confint.poss(rez.gel, type = "invLR", 
                             fact = fact, level = conflev.TM)    
    #GEL, CUE
    rez.cue <- gel(g.means, x = cbind(yt, zt, x1t),
                   tet = tet.init,
                   type = "CUE")
    tet3[i, ] <- rez.cue$coefficients
    ci3[i, ] <- confint.poss(rez.cue, type = "invLR", 
                             fact = fact, level = conflev.TM) 
    #GEL, ET
    rez.et <- gel(g.means, x = cbind(yt, zt, x1t),
                  tet = tet.init,
                  type = "ET")
    tet4[i, ] <- rez.et$coefficients
    ci4[i, ] <- confint.poss(rez.et, level = conflev.TM,
                             type = "invLR", fact = fact)
    #GEL, HD
    rez.hd <- gel(g.means, x = cbind(yt, zt, x1t),
                  tet = tet.init,
                  type = "HD")
    tet5[i, ] <- rez.hd$coefficients
    ci5[i, ] <- confint.poss(rez.hd, level = conflev.TM,
                             type = "invLR", fact=fact)
  }
  ### Calculate the simulated performance values- bias, variance, RMSE
  bias <- sapply(list(tet1, tet2, tet3, tet4, tet5), 
                 FUN=function(i) rowMeans(t(i) - tet0, na.rm = FALSE))
  Var <- sapply(list(tet1, tet2, tet3, tet4, tet5), 
                FUN=function(i) diag(var(i, na.rm = FALSE)))
  RMSE <- sapply(list(tet1, tet2, tet3, tet4, tet5), 
                 FUN=function(i) sqrt(rowMeans((t(i) - tet0)^2, na.rm = FALSE)))
  dimnames(Var) <- dimnames(bias) <- dimnames(RMSE) <-
    list(c("delta", "mu1", "mux"), c("LS", "GEL.EL", "GEL.CUE", "GEL.ET", "GEL.HD"))
  ### Calculate the simulated CI performance - coverage and length
  coverage <- array(sapply(list(ci1, ci2, ci3, ci4, ci5), 
                           FUN=function(i) mean(sign((i[, 1] - tet0[1])*i[, 2] - tet0[1])==-1, 
                           na.rm = T)), dim=c(1,5))
  ci.length <- array(sapply(list(ci1, ci2, ci3, ci4, ci5), 
                            FUN=function(i) mean(i[,2] - i[,1],
                            na.rm = T)), dim=c(1,5))
  num.nas = array(sapply(list(tet1, tet2, tet3, tet4, tet5), 
                         FUN=function(i) sum(is.na(i[,1]))), dim=c(1,5))
  dimnames(num.nas) <- dimnames(coverage) <- dimnames(ci.length) <- 
    list(c("delta"),  c("LS", "GEL.EL", "GEL.CUE",  "GEL.ET", "GEL.HD"))
  
  return(list(bias = bias, Variance = Var, RMSE = RMSE, 
              coverage = coverage, ci.length = ci.length, num.nas = num.nas)) 
}
```

Finally, run the function. Note that setting a higher `fact` safeguards
against failing to calculate confidence intervals, but runs longer. Note
that in our example none of the CIs failed for non of the methods, as
`num.nas` is 0.

**GEL for trimmed means**

``` r
simrez1 <- sim_1cov(n, z, data.norm, alpha = 0.1, conf.nominal = 0.95, fact = 4)
simrez1
```

    ## $bias
    ##                LS      GEL.EL     GEL.CUE      GEL.ET      GEL.HD
    ## delta 0.002588005 0.008726018 0.008844041 0.008793116 0.008762487
    ## mu1   0.008837231 0.003264982 0.003256948 0.003256707 0.003259793
    ## mux   0.256157794 0.001238320 0.001218255 0.001213738 0.001222608
    ## 
    ## $Variance
    ##                LS      GEL.EL     GEL.CUE      GEL.ET      GEL.HD
    ## delta 0.007966114 0.017472693 0.017502054 0.017484234 0.017477659
    ## mu1   0.017502382 0.008664294 0.008681239 0.008671141 0.008667398
    ## mux   0.003229040 0.005314090 0.005290636 0.005296562 0.005303804
    ## 
    ## $RMSE
    ##               LS     GEL.EL    GEL.CUE     GEL.ET     GEL.HD
    ## delta 0.08906728 0.13214187 0.13226020 0.13218976 0.13216297
    ## mu1   0.13226098 0.09290658 0.09299699 0.09294295 0.09292302
    ## mux   0.26235417 0.07272588 0.07256491 0.07260545 0.07265521
    ## 
    ## $coverage
    ##         LS GEL.EL GEL.CUE GEL.ET GEL.HD
    ## delta 0.86   0.94   0.955   0.95   0.94
    ## 
    ## $ci.length
    ##              LS    GEL.EL   GEL.CUE    GEL.ET    GEL.HD
    ## delta 0.3895847 0.5275152 0.5411908 0.5276128 0.5258986
    ## 
    ## $num.nas
    ##       LS GEL.EL GEL.CUE GEL.ET GEL.HD
    ## delta  0      0       0      0      0

**GEL for means**

Note, that the same experiment can be used to calculate the GEL for mean
treatment effect, just by setting the trimming constant `alpha` to zero:

``` r
simrez2 <- sim_1cov(n, z, data.norm, alpha = 0, conf.nominal = 0.95, fact = 4)
```

### GEL for Huber estimators

``` r
sim_hub <- function(n, z, data, tet0 = c(0, 0, 0),
                conf.nominal = 0.95, fact = 4)  {
 ### Input
  # n - sample size
  # z - treatment group indicator (common for all replicates)
  # data - list of B data replicates, where each element
  # tet0 - theta value under H0 (for bias, RMSE and CI coverage)
  # fact - et factor for CI search span in gel.confint()
  
  ### Output
  # bias, variance, RMSE for each theta element
  # coverage, ci.length - performance of confidence intrvals for each of theta elements
  # num.nas - for how many replicates gel.confint() failed to estimate CI
  
  ### Calculations
  # Define arrays to store simulated parameter estimates
  #   and their CIs for 5 different methods (LS, EL, ET, CUE, 
  tet5 <- tet4 <- tet3 <- tet2 <- tet1 <- matrix(0, length(data), 3)
  ci5 <- ci4 <-ci3 <- ci2 <- ci1 <- matrix(0, length(data), 2)
  
  for(i in 1: length(data)) {
    XX <- data[[i]] 
    y <- XX[, 1]; x1 <- XX[, 2]
    tet.init <- c(mean(y[z == 1]) - mean(y[z == 0]),
                  mean(y[z == 0]), 
                  smhuber(x1)$mu)
    #Scale estimator for Smooth Huber
    s <- mad(x1)
    #Define moment conditions
    gh <- function(tet, W)
    {
      y <- W[, 1]
      z <- W[, 2]
      x <- W[, 3]
      m1 <- y - tet[2] - tet[1] * z
      m2 <- z * (y - tet[2] - tet[1] * z)
      m3 <- smpsi((x - tet[3]) / s) 
      m4 <- z * smpsi((x - tet[3]) / s)
      f <- cbind(m1, m2, m3, m4)
      return(f)
    }
    ### ANCOVA
    #LM
    rez.lm <- lm(y ~ z + x1) 
    tet1[i, ] <- coefficients(rez.lm)[1:3]
    ci1[i, ] <- confint(rez.lm)[2, ] #ci for z
    #GEL, EL
    rez.gel <- gel(gh,
                   x = cbind(y, z, x1),
                   tet = tet.init)
    tet2[i, ] <- rez.gel$coefficients
    ci2[i, ] <- confint.poss(rez.gel, type = "invLR", fact = fact)    
    #GEL, CUE
    rez.cue <- gel(gh, 
                   x = cbind(y, z, x1),
                   tet = tet.init,
                   type = "CUE")
    tet3[i, ] <- rez.cue$coefficients
    ci3[i, ] <- confint.poss(rez.cue, type = "invLR", fact = fact) 
    
    rez.et <- gel(gh, 
                  x = cbind(y, z, x1),
                  tet = tet.init,
                  type = "ET")
    tet4[i, ] <- rez.et$coefficients
    ci4[i, ] <- confint.poss(rez.et,
                             type = "invLR", fact = fact)
    
    rez.hd <- gel(gh, 
                  x = cbind(y,z, x1),
                  tet = tet.init,
                  type = "HD")
    tet5[i, ] <- rez.hd$coefficients
    ci5[i, ] <- confint.poss(rez.hd, 
                             type = "invLR", fact = fact)
  }
  bias <- sapply(list(tet1, tet2, tet3, tet4, tet5), 
                 FUN = function(i) rowMeans(t(i) - tet0, na.rm = FALSE))
  Var <- sapply(list(tet1, tet2, tet3, tet4, tet5), 
                FUN = function(i) diag(var(i, na.rm = FALSE)))
  RMSE <- sapply(list(tet1, tet2, tet3, tet4, tet5), 
                 FUN = function(i) sqrt(rowMeans((t(i) - tet0)^2, na.rm = FALSE)))
  dimnames(Var) <- dimnames(bias) <- dimnames(RMSE) <-
    list(c("delta", "mu1", "mux"), c("LS", "GEL.EL", "GEL.CUE", "GEL.ET", "GEL.HD"))
  
  coverage <- array(sapply(list(ci1, ci2, ci3, ci4, ci5), 
                           FUN=function(i) mean(sign((i[, 1] - tet0[1])*i[, 2] - tet0[1])==-1, na.rm = T)), dim=c(1,5))
  ci.length <- array(sapply(list(ci1, ci2, ci3, ci4, ci5), 
                            FUN=function(i) mean(i[, 2] - i[, 1], na.rm = T)), dim=c(1, 5))
  num.nas = array(sapply(list(tet1, tet2, tet3, tet4, tet5), 
                         FUN=function(i) sum(is.na(i[, 1]))), dim=c(1, 5))
  dimnames(num.nas) <-    dimnames(coverage) <- dimnames(ci.length) <- 
    list(c("delta"),  c("LS", "GEL.EL", "GEL.CUE",  "GEL.ET", "GEL.HD"))
  
  return(list(bias = bias, Variance = Var, RMSE = RMSE, 
              coverage = coverage, ci.length = ci.length, num.nas = num.nas)) 
}
```

Run the function:

``` r
simrez3 <- sim_hub(n, z, data.norm, 
                conf.nominal = 0.95, fact = 4)
```

### Compare the results

Let’s compare the estimation of $\Delta$ with the three methods above.

``` r
RMSE <- rbind(simrez1$RMSE["delta", ], simrez2$RMSE["delta", ], simrez3$RMSE["delta", ]) %>%  as.data.frame()
RMSE$method <- c("TM", "Means", "Huber")

bias <- rbind(simrez1$bias["delta", ], simrez2$bias["delta", ], simrez3$bias["delta", ])  %>%  as.data.frame()
bias$method <- c("TM", "Means", "Huber")

length <- rbind(simrez1$ci.length, simrez2$ci.length, simrez3$ci.length)  %>%  as.data.frame()
length$method <- c("TM", "Means", "Huber")
```

``` r
RMSE %>% 
  pivot_longer(1:5, names_to = "type", values_to = "RMSE") %>% 
ggplot(aes(x = type, y = RMSE)) +
  geom_jitter(aes(color = method), height = 0, width = 0.1) +
  theme_bw()
```

![Comparison of RMSE of Delta for various covariate adjustment
methods](RobCovGEL_files/figure-gfm/unnamed-chunk-23-1.png)

![Comparison of bias of Delta for various covariate adjustment
methods](RobCovGEL_files/figure-gfm/unnamed-chunk-24-1.png)

![Comparison of the length of CI of Delta for various covariate
adjustment methods](RobCovGEL_files/figure-gfm/unnamed-chunk-25-1.png)

## References

-   Chaussé, P. (2010). Computing Generalized Method of Moments and
    Generalized Empirical Likelihood with R. Journal of Statistical
    Software, 34(11), 1–35. (<https://doi.org/10.18637/jss.v034.i11>)  
-   Chaussé, P., Jin L., and Luta, G. (2016). A Simulation-Based
    Comparison of Covariate Adjustment Methods for the Analysis of
    Randomized Controlled Trials. International Journal of Environmental
    Research and Public Health 13(4), 414.
    (<https://doi.org/10.3390/ijerph13040414>)
