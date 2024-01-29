### # 29.01.2024 Mara Delesa-Velina

trim.index <- function(Y, alpha, beta) {
  ### One-sample trimming and the trimming indices
  ### Input
  # Y: sample to trim
  # a, b: Set lower and upper trimming constants a, b < 0.5
  ### Output: 
  # x - trimmed data; 
  # keep.ind - observation indices of the sorted sample that are contained in the trimmed sample
  # keep - logical vector indicating which observations of original non-sorted sample are kept
  n <- length(Y)
  r <- trunc(n * alpha) + 1
  s <- n - trunc(n * beta)
  m <- s - r + 1
  ys <- sort(Y)
  x <- ys[r : s]
  keep <- order(Y) >= r & order(Y) <= s 
  list(x = x, keep.ind = r : s, keep = keep)
}

a_abhat <-function(Y, alpha, beta) {
  ### Function for calculating the scaling factor a for the limiting distribution of trimmed GEL statistic
  ### Input
  # Y: sample to trim
  # alpha, beta: Set lower and upper trimming constants a, b < 0.5
  ### Output 
  # mu_ab: sample trimmed mean estimate
  # sigma2, tau2: squared sigma and tau constants
  # aconst: scaling factor a
  # m: effective number of observations after trimming
  n <- length(Y)
  r <- trunc(n * alpha) + 1
  s <- n - trunc(n * beta)
  m<- s - r + 1
  Y_tr <- sort(Y)[r : s]
  mu_ab <- mean(Y_tr)
  xi1 <- Y_tr[1]
  xi2 <- tail(Y_tr, 1)
  sigma2 <- sum((Y_tr - mu_ab)^2) / n / (1 - alpha - beta)
  tau2 <- sigma2 / (1 - beta - alpha) + (beta * (1 - beta) * (xi2 - mu_ab)^2 -
          2 * alpha * beta * (xi1 - mu_ab) * (xi2 - mu_ab) + 
          alpha * (1 - alpha) * (xi1 - mu_ab)^2) / 
         (1 - beta - alpha)/(1 - beta - alpha)
  a <- sigma2 / (1 - beta - alpha) / tau2
  return(list("mu_ab" = mu_ab, "sigma2" = sigma2,
              "tau2" = tau2, "aconst" = a, "m" = m))
}


tmeans.conf.level <- function(x, y, mu.t = 0, alpha = 0.1, 
                              beta = 0.1, conf.level = 0.95, Fapprox = F, 
                              q = 2, p = 1){
  ### Function for calculating corrected confidence level for 2-sample comparison corresponding to the scaling factor
  ### Function allows approximation of Chi square distribution by F distribution for small sample sizes
  ## Input
  # x, y - samples 1 and 2 to compare
  # mu.t - trimmed mean value under H0
  # alpha, beta - set lower and upper trimming constants a, b < 0.5
  # conf.level - test confidence level
  # Fapprox - should an F approximation be used
  # q, p - dof for the F test, where 
  # q is number of moment conditions, p is parameter vector dimension
  if (alpha == 0 & beta == 0) {
    conf.level.tmF <- conf.level.tm <- conf.level
  } else {
    n1 <- length(x)
    n2 <- length(y)
    scaling_1 <- a_abhat(x, alpha, beta)
    scaling_2 <- a_abhat(y, alpha, beta)
    m1 <- scaling_1$m
    m2 <- scaling_2$m
    a_scaling <- (n1/m1) * (n2/m2) * (m2 * scaling_1$sigma2 + m1 * scaling_2$sigma2)/
      (n2 * scaling_1$tau2 + n1*scaling_2$tau2)
    #Corrected confidence level for TM GEL test
    conf.level.tm <- pchisq(qchisq(conf.level, q-p) / a_scaling, q - p)
    n0 <- min(m1, m2) - q - p
    conf.level.tmF <- pf(qf(conf.level, p, n0)/a_scaling, 1, n0)
  }
  return(ifelse(Fapprox, conf.level.tmF, conf.level.tm))
}

