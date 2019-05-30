# PC procedure factor model
pcfest <- function(X, r, demean = 0, constant = 1){
                       
  if(demean == 0)
    X <- X

  # demean data
  if(demean == 1)
    X <- scale(X, scale = FALSE)

  # standardized data
  if(demean == 2)
    X <- scale(X, scale = TRUE)
    
  # standardized ahn and horenstein
  if(demean == 3){
    X <- scale(X, scale = FALSE) 
    X <- t(sapply(1:nrow(X), function(x) X[x,] - rowMeans(X)[x]))
    X <- X + mean(X)
  }
  
  # cross-section and times
  n <- ncol(X)
  tt <- nrow(X)
                                                
  # if r
  if(r > 0){
    # estimating the r factors F'F/T=I
    if(tt <= n){
      # singular value decomposition
      XXt <- X%*%t(X)
      SVD <- eigen(XXt)
      index <- order(SVD$values[1:r], decreasing = TRUE)

      # eigen vectros
      L <- SVD$vectors[,index]
      # factors
      Fhat <- L*sqrt(tt)
      # loadings
      lambda <- t(X)%*%Fhat/tt
      # normalization-1
      norm <- round(t(Fhat)%*%Fhat/tt)
    }else{
      # estimating the r factors L'L/N=I
      # singular value decomposition
      XtX <- t(X)%*%X
      SVD <- eigen(XtX)
      index <- order(SVD$values[1:r], decreasing = TRUE)

      # eigen vectors
      L <- SVD$vectors[,index]
      # loadings
      lambda <- L*sqrt(n)
      # factors
      Fhat <- X%*%lambda/n
      # normalization-1
      norm <- round(t(lambda)%*%lambda/n)
    }
    # matrices
    Fhat <- as.matrix(Fhat)
    colnames(Fhat) <- paste("f",1:r, sep = "")
    # common component and idiosyncratic component
    chat <- Fhat%*%t(lambda)
    ehat <- X-chat
  }else{
    # without factors
    norm <- Fhat <- lambda <- chat <- 0
    ehat <- X
  }

  # sum of squared residuals (divided by NT)
  Vfk <- mean(colMeans(ehat^2))
    
  # test of number of factors
  k <- r
  CNT2 <- min(c(n,tt))

  # criteria
  ICP1k <- log(Vfk) + constant*(k*((n+tt)/(n*tt)))*log((n*tt)/(n+tt))
  ICP2k <- log(Vfk) + constant*(k*((n+tt)/(n*tt)))*log(CNT2)
  ICP3k <- log(Vfk) + constant*(k*((log(CNT2)/CNT2)))
  
  # vector
  ICPk <- c(ICP1k, ICP2k, ICP3k)
  names(ICPk) <- c("ICP1k", "ICP2k", "ICP3k")
  
  # result
  result <- list(X, norm, ICPk, Fhat, lambda, ehat, chat, Vfk)
  names(result) <- c("X", "norm", "ICPk", "Fhat", "lambda", "ehat", "chat", 
                    "Vfk")
  class(result) <- "dfm"
                                                 
  return(result)
}

# determinating the number of dynamic factors
dfactest.m <- function(X, r, p = 2, method = 1, delta = 1/10, m = 1,
                demean = 0){

  # raw data
  if(demean == 0)
    X <- X

  # demean data
  if(demean == 1)
    X <- scale(X, scale = FALSE)

  # standardized data
  if(demean == 2)
    X <- scale(X, scale = TRUE)
    
  # standardized ahn and horenstein
  if(demean == 3){
    X <- scale(X, scale = FALSE) 
    X <- t(sapply(1:nrow(X), function(x) X[x,] - rowMeans(X)[x]))
    X <- X + mean(X)
  }
  
  # dimensions
  n <- ncol(X)
  tt <- nrow(X)
  
  # stati factors
  Fhat <- pcfest(X, r = r, norma = method)$Fhat
  
  # VAR(p)
  var1 <- VAR(Fhat, tXpe = "none", p = p)
  e <- resid(var1)
  
  # if covariance matrix
  if(method == 1){
    sigF <- cov(e)
  }
  
  # if correlation matrix
  if(method == 2){
    sigF <- cor(e)
  }
  
  # spectral decomposition
  SVD <- eigen(sigF)
  s <- SVD$values
  
  # squared eigen vectors
  c2 <- s^2
  nf <- sum(c2)
   
  # statistics
  stat <- matrix(0,r-1,2)
  colnames(stat) <- c("stat1", "stat2")
   
  for(k in 1 : (r-1)){   
    stat[k,"stat1"] <- sqrt(c2[k+1]/nf)
    stat[k,"stat2"] <- sqrt(sum(c2[(k+1):r])/nf)
  }
 
  # criterion decisions   
  crit <- m/min(c(n^(1/2-delta),tt^(1/2-delta)))
  Kappa <- stat < crit
  colnames(Kappa) <- c("K3", "K4")
  
  # number of dXnamic factor models
  if(sum(Kappa[,"K3"]) > 0){
    q3 <- min(which(Kappa[,"K3"]))  
  }else{
    q3 <- 0
  }
  
  if(sum(Kappa[,"K4"]) > 0){
    q4 <- min(which(Kappa[,"K4"]))
  }else{
    q4 <- 0
  }
  
  # result
  result <- c(q3, q4)
  names(result) <- c("q3", "q4")
  
  return(result)
}


# PC procedure factor model (non-stationary)
pcbai <- function(X, r, kmax = 8, demean = 0){

  # raw data
  if(demean == 0)
    X <- X

  # demean data
  if(demean == 1)
    X <- scale(X, scale = FALSE)

  # standardized data
  if(demean == 2)
    X <- scale(X, scale = TRUE)
    
  # standardized ahn and horenstein
  if(demean == 3){
    X <- scale(X, scale = FALSE) 
    X <- t(sapply(1:nrow(X), function(x) X[x,] - rowMeans(X)[x]))
    X <- X + mean(X)
  }

  # cross-section and times
  n <- ncol(X)
  tt <- nrow(X)
  
  # estimating the r factors F'F/T^2=I
  if(n > tt){
    # singular value decomposition
    XXt <- X%*%t(X)
    SVD <- eigen(XXt)
    index <- order(SVD$values[1:kmax], decreasing = TRUE)

    # eigen vectros
    L <- SVD$vectors[,index]
    # factors
    Fhatmax <- L*tt
    # loadings
    lambdamax <- (t(X)%*%Fhatmax)/tt^2
    # normalization-1
    norm <- round(t(Fhatmax)%*%Fhatmax/tt^2)
  }else{
    # estimating the r factors L'L/N=I
    # singular value decomposition
    XtX <- t(X)%*%X
    SVD <- eigen(XtX)
    index <- order(SVD$values[1:kmax], decreasing = TRUE)

    # eigen vectros
    lambdamax <- SVD$vectors[,index]*sqrt(n)
    # loadings
    Fhatmax <- X%*%lambdamax/n
    # normalization-1
    norm <- round(t(lambdamax)%*%lambdamax/n)
  }
  
  # matrices
  Fhatmax <- as.matrix(Fhatmax)

  Fhat <- as.matrix(Fhatmax[,1:r])
  colnames(Fhat) <- paste("f",1:r, sep = "")

  lambda <- lambdamax[,1:r]
    
  # test of number of factors
  k <- r 
  CNT2 <- min(c(n,tt))
  
  # common components
  chat <- Fhat%*%t(lambda)
  chatmax <- Fhatmax%*%t(lambdamax)
  
  # idiosyncratic components
  ehat <- X - chat
  ehatmax <- X - chatmax

  # alpha
  alpha <- tt/(4*log(log(tt)))

  # sigma max
  sigmamax <- mean(colMeans(ehatmax^2))

  # sum squared resids
  Vfk <- mean(colMeans(ehat^2))
  
  # criteria
  IPC1k <- Vfk + ((k*sigmamax*alpha)*((n+tt)/(n*tt)))*log((n*tt)/(n+tt))
  IPC2k <- Vfk + ((k*sigmamax*alpha)*((n+tt)/(n*tt)))*log(CNT2)
  IPC3k <- Vfk + (k*sigmamax*alpha)*((log(CNT2)/CNT2))
  
  # vector
  IPCk <- c(IPC1k, IPC2k, IPC3k)
  names(IPCk) <- c("IPC1k", "IPC2k", "IPC3k")
  
  # result
  result <- list(X, norm, IPCk, Fhat, lambda, ehat, chat)
  names(result) <- c("X", "norm", "IPCk", "Fhat", "lambda", "ehat", "chat")
  class(result) <- "dfm"
  
  return(result)
}


# number of factors onatski 2010 
onatski2010 <- function(X, demean){  

  # raw data
  if(demean == 0)
    X <- X

  # demean data
  if(demean == 1)
    X <- scale(X, scale = FALSE)

  # standardized data
  if(demean == 2)
    X <- scale(X, scale = TRUE)
    
  # standardized ahn and horenstein
  if(demean == 3){
    X <- scale(X, scale = FALSE) 
    X <- t(sapply(1:nrow(X), function(x) X[x,] - rowMeans(X)[x]))
    X <- X + mean(X)
  }
                    
  # transpose
  X <- t(X)
  
  # T and N
  T <- ncol(X)
  N <- nrow(X)
  
  # rnmax (onatski wp)
  rmax <- round(1.55*min(T^(2/5), N^(2/5)))

  # j initial
  j <- rmax + 1

  # singular value decomposition
  XXt <- X%*%t(X)
  SVD <- eigen(XXt/T)

  # lambda
  lambda <- SVD$values
  
  # inital values
  rhat <- Inf
  conver <- Inf
  i <- 1
  Conver <- Inf

  while(length(rhat) < 4 | conver > 0 | conver < 0){
    # iterative
    
    # betahat
    betahat <- 2*abs(coef(lm(lambda[j:(j+4)] ~ c(((j - 1):(j + 3))^(2/3))))[2])
    
    # criterion (page 1008)
    cond <- which(lambda[-length(lambda)] - lambda[-1] >= betahat)
            
    if(length(cond) != 0){
      rHat <- max(cond)
      rhat[i] <- rHat  

      # repeat 2 and 3
      j <- rhat[i] + 1
      # convergence                               
      if(length(rhat) > 4)
        conver <- sum(diff(rhat[(length(rhat)-4):length(rhat)]))
    }else{
      rhat <- rep(0, i)
      conver <- 0
    }    
    
    i <- i + 1
  }    

  # rhat
  r.hat <- c(rhat[length(rhat)], betahat)
  names(r.hat) <- c("ed", "betahat")
  
  # return
  return(r.hat)
}

# principal components + kalman filter  
pckf <- function(Y, q, s = 0, p = 1, demean = 0, 
          type.adf = c("none", "const", "trend"), alpha = 0.05){
  
  # raw data
  if(demean == 0)
    Y <- Y
  
  # demean data
  if(demean == 1)
    Y <- scale(Y, scale = FALSE)
  
  # standardized data
  if(demean == 2)
    Y <- scale(Y, scale = TRUE)
  
  # r static factors
  r <- q*(s+1)
  
  # extract the factors with PC
  dfm <- pcfest(Y, r = r)
  
  # factor
  z <- dfm$Fhat
  
  # adf test specification
  spe <- adf.test <- matrix(0, r, 2)
  colnames(spe) <- colnames(adf.test) <- c("level", "diff")
  
  # select specification
  for(i in 1 : r){  
      bic.i <- matrix(0, ncol(adf.test), length(type.adf))
      colnames(bic.i) <- type.adf
      rownames(bic.i) <- colnames(adf.test)
      
    for(j in 1 : length(type.adf)){
      type.j <- type.adf[j]
      bic.i[1,j] <- adf(z[,i], type.j)$bic  
      bic.i[2,j] <- adf(diff(z[,i]), type.j)$bic  
    
      spe[i,"level"] <- type.adf[which(bic.i[1,] == min(bic.i["level",]))][1] 
      spe[i,"diff"] <- type.adf[which(bic.i[2,] == min(bic.i["diff",]))][1]
    }
  }
        
  # adf test                                        
  adf.test[,"level"] <- sapply(1:ncol(z), function(x) adf(z[,x], 
                        spe[x,"level"])$p.value)
  adf.test[,"diff"] <- sapply(1:ncol(z), function(x) adf(diff(z[,x]), 
                              spe[x,"diff"])$p.value)
  
  # test
  for(i in 1 : r){
    adf.test[i,"level"] <- adf.test(z[,i])$p.value
    adf.test[i,"diff"] <- adf.test(diff(z[,i]))$p.value
  }
  
  # N and T
  N <- ncol(Y)
  T <- nrow(Y)
  
  # stationary?
  ind <- adf.test[,"level"] < alpha
  
  # coeffcient matrix of autoregressive process of factors
  Atemp <- matrix(0, r, r*p)
  initV <- matrix(0, r, r)
  
  # epsilon
  e <- matrix(0, T - p, r)
                                                
  if(sum(ind) > 0){                    
    # coefficients
    for(i in which(ind)){
      arp <- ar.ols(z[,i], FALSE, p, demean = FALSE)
      Atemp[i, i] <- arp$ar[,,1]
      e[,i] <- matrix(arp$resid[-(1:p)])
    }
  }
           
  if(sum(!ind) > 0){
    # coefficients
     for(i in which(!ind)){
       Atemp[i, i] <- 1
       e[,i] <- diff(z[,i])
       initV[i,i] <- 10^7
     }
  }
  
  # coefficient matrix
  A <- Atemp
                  
  # @ variance of idiosyncratic component
  if(r > 1){
    # covariance matrix
    H <- matrix(0,r,r)
    for(i in 1 : r)
      H[i,i] <- var(e[,i])
  }else{
    # variance
    H <- var(e)
  }
  
  # % if s is different from 0 we have q dynamic factors
  Q <- matrix(0, r*p, r*p)
  if(r > q){
    # extract the first q eigenvectors and eigenvalues from cov(e)
    eid <- eigen(H)
    M <- diag(q)
    diag(M) <- eid$values[1:q]
    P <- eid$vectors[,1:q,drop=FALSE]
    # extract the common shocks
    u_orth <- e%*%P%*%M^(-1/2)
    e_pc <- e%*%P%*%t(P)
    # variance of the VAR shock when s>0
    Q[1:r, 1:r] <- P%*%M%*%t(P)
  }else{
    # variance of the VAR shock when s=0 (rename the covariance)
    Q[1:r, 1:r] <- H
  }
  
  # initial covariance matrix for stationary case
  if(sum(ind) > 0){                    
    # coefficients
    for(i in which(ind))
      initV[i,i] <- Q[i,i]/(1-arp$ar[,,1]^2)
  }
  
  # R diagonal
  v <- dfm$lambda
  chi <- z%*%t(v)
  R <- diag(diag(cov(Y - chi)))
  #R <- cov(Y - chi)
  
  # initial state mean
  initx <- rep(0, r*p) 
                    
  # loading matrix
  C <- matrix(0, N, r*p)
  C[,1:r] <- v
  
  # kalman model
  kf <- dlm(FF = C, V = R, GG = A, W = Q, m0 = initx, C0 = initV)
  
  # smoothing
  if(r > 1){
    Fkf <- dlmSmooth(Y, kf)$s[-1,1:r]
  }else{
    Fkf <- as.matrix(dlmSmooth(Y, kf)$s[-1])
  }
  
  # filter
  Ykf <- dlmFilter(Y, kf)$y
  
  # result
  result <- list(C, R, A, Q, initx, initV, z, Fkf, Ykf)
  names(result) <- c("C", "R", "A", "Q", "initx", "initV", "Fpc", "Fkf", "Ykf")

  return(result)
}


# my standard deviation
mysd <- function(x){
  x <- c(x)
  sdu <- sqrt(sum((x-mean(x))^2)/length(x))
  return(sdu)
}

# eigenvalue ratio test
ratio.test <- function(X, kmax = NULL, demean = 3){

  # raw data
  if(demean == 0)
    X <- X

  # demean data
  if(demean == 1)
    X <- scale(X, scale = FALSE)

  # standardized data
  if(demean == 2)
    X <- scale(X, scale = TRUE)
    
  # standardized ahn and horenstein
  if(demean == 3){
    X <- scale(X, scale = FALSE) 
    X <- t(sapply(1:nrow(X), function(x) X[x,] - rowMeans(X)[x]))
    X <- X + mean(X)
  }

  # T and N
  T <- nrow(X)
  N <- ncol(X)                      
  
  # extracting the eigenvalues
  if(T > N){
    mu <- eigen(t(X)%*%X/(N*T))$values
  }else{
    mu <- eigen(X%*%t(X)/(N*T))$values
  }
  
  # Vfk
  Vfk <- c()
  for(i in 1 : c(kmax + 1))        
    Vfk[i] <- sum(mu[i:min(N,T)])
    
  # mock eigenvalue
  mu0 <- sum(mu[1:(min(N,T))])/(log(min(N,T))*N)
  mu <- c(mu0, mu)
         
  # test ER(k)
  ER <- mu[1:c(kmax+2)][-c(kmax+2)]/mu[1:c(kmax+2)][-1]
  if(is.null(Vfk)){
    # test GR(k)
    GR <- NULL
  }else{
    # test GR(k)
    GR <- log(1 + mu[1:c(kmax+2)][-c(kmax+2)]/c(Vfk))/
          log(1 + mu[1:c(kmax+2)][-1]/c(Vfk))
  }
  # estimators
  ker <- which(ER == max(ER)) - 1
  kgr <- which(GR == max(GR)) - 1
  
  # result
  result <- c(ker, kgr, mu)
  names(result) <- c("ker", "kgr", paste("mu_",0:min(N,T), sep = ""))
  # return
  return(result)
}

# pooled test
pooled.test <- function(e, type = "const"){
  N <- ncol(e)
  Pest <- (-2*sum(apply(e, 2, function(x) log(adf(x, type = type)$p.value)))
            - 2*N)/sqrt(4*N)
  p.value <- 1 - pnorm(Pest)

  result <- round(c(Pest, p.value), 4)
  names(result) <- c("Pest", "p.value")

  return(result)
}


# adf test
adf <- function(x, type = c("none", "const", "trend"),
                alternative = c("stationary", "explosive"),
                k = trunc((length(x) - 1)^(1/3))){

    if (NCOL(x) > 1)
        stop("x is not a vector or univariate time series")
    if (any(is.na(x)))
        stop("NAs in x")
    if (k < 0)
        stop("k negative")

    type <- match.arg(type)
    alternative <- match.arg(alternative)
    DNAME <- deparse(substitute(x))

    k <- k + 1
    x <- as.vector(x, mode = "double")
    y <- diff(x)
    n <- length(y)
    z <- embed(y, k)
    yt <- z[, 1]
    xt1 <- x[k:n]
    tt <- k:n

    if(type == "none"){
      if(k > 1){
        yt1 <- z[, 2:k]
        res <- lm(yt ~ xt1 + yt1 - 1)
      }else{
        res <- lm(yt ~ xt1 - 1)
      }
      res.sum <- summary(res)
      STAT <- res.sum$coefficients[1, 1]/res.sum$coefficients[1,
                2]
      Table <- adfTable(trend = "nc", statistic = "t")
    }

    if(type == "const"){
      if(k > 1){
        yt1 <- z[, 2:k]
        res <- lm(yt ~ xt1 + 1 + yt1)
      }else{
        res <- lm(yt ~ xt1 + 1)
      }
      res.sum <- summary(res)
      STAT <- res.sum$coefficients[2, 1]/res.sum$coefficients[2,
                2]
      Table <- adfTable(trend = "c", statistic = "t")
    }

    if(type == "trend"){
      if (k > 1) {
        yt1 <- z[, 2:k]
        res <- lm(yt ~ xt1 + 1 + tt + yt1)
      }else{
        res <- lm(yt ~ xt1 + 1 + tt)
      }
      res.sum <- summary(res)
      STAT <- res.sum$coefficients[2, 1]/res.sum$coefficients[2,
                2]
      Table <- adfTable(trend = "ct", statistic = "t")
    }
    bic <- BIC(res)
    table <- Table$z
    tablen <- dim(table)[2]
    tableT <- Table$x
    tablep <- Table$y
    tableipl <- numeric(tablen)
    for (i in (1:tablen)) tableipl[i] <- approx(tableT, table[,
        i], n, rule = 2)$y
    interpol <- approx(tableipl, tablep, STAT, rule = 2)$y
    if (is.na(approx(tableipl, tablep, STAT, rule = 1)$y))
        if (interpol == min(tablep))
            warning("p-value smaller than printed p-value")
        else warning("p-value greater than printed p-value")
    if (alternative == "stationary")
        PVAL <- interpol
    else if (alternative == "explosive")
        PVAL <- 1 - interpol
    else stop("irregular alternative")
    PARAMETER <- k - 1
    METHOD <- "Augmented Dickey-Fuller Test"
    names(STAT) <- "Dickey-Fuller"
    names(PARAMETER) <- "Lag order"
    structure(list(bic = bic, statistic = STAT, parameter = PARAMETER,
        alternative = alternative, p.value = PVAL, method = METHOD,
        data.name = DNAME), class = "htest")
}

# adf tables
adfTable <- function(trend = c("nc", "c", "ct"), statistic = c("t", "n"),
           includeInf = TRUE)
  {
    # A function implemented by Diethelm Wuertz

    # Description:
    #   Tables critical values for augmented Dickey-Fuller test.

    # Note:
    #   x=-3:0; y=0:3; z=outer(x,y,"*"); rownames(z)=x; colnames(z)=y; z

    # Examples:
    #   adfTable()

    # FUNCTION:

    # Match Arguments:
    type = trend = match.arg(trend)
    statistic = match.arg(statistic)

    # Tables:
    if (statistic == "t") {
      # Hamilton Table B.6 - OLS t-Statistic
      if (type == "nc") {
        table = cbind(
          c(-2.66, -2.26, -1.95, -1.60, +0.92, +1.33, +1.70, +2.16),
          c(-2.62, -2.25, -1.95, -1.61, +0.91, +1.31, +1.66, +2.08),
          c(-2.60, -2.24, -1.95, -1.61, +0.90, +1.29, +1.64, +2.03),
          c(-2.58, -2.23, -1.95, -1.62, +0.89, +1.29, +1.63, +2.01),
          c(-2.58, -2.23, -1.95, -1.62, +0.89, +1.28, +1.62, +2.00),
          c(-2.58, -2.23, -1.95, -1.62, +0.89, +1.28, +1.62, +2.00))
      } else if (type == "c") {
        table = cbind(
          c(-3.75, -3.33, -3.00, -2.63, -0.37, +0.00, +0.34, +0.72),
          c(-3.58, -3.22, -2.93, -2.60, -0.40, -0.03, +0.29, +0.66),
          c(-3.51, -3.17, -2.89, -2.58, -0.42, -0.05, +0.26, +0.63),
          c(-3.46, -3.14, -2.88, -2.57, -0.42, -0.06, +0.24, +0.62),
          c(-3.44, -3.13, -2.87, -2.57, -0.43, -0.07, +0.24, +0.61),
          c(-3.43, -3.12, -2.86, -2.57, -0.44, -0.07, +0.23, +0.60))
      } else if (type == "ct") {
        table = cbind(
          c(-4.38, -3.95, -3.60, -3.24, -1.14, -0.80, -0.50, -0.15),
          c(-4.15, -3.80, -3.50, -3.18, -1.19, -0.87, -0.58, -0.24),
          c(-4.04, -3.73, -3.45, -3.15, -1.22, -0.90, -0.62, -0.28),
          c(-3.99, -3.69, -3.43, -3.13, -1.23, -0.92, -0.64, -0.31),
          c(-3.98, -3.68, -3.42, -3.13, -1.24, -0.93, -0.65, -0.32),
          c(-3.96, -3.66, -3.41, -3.12, -1.25, -0.94, -0.66, -0.33))
      } else {
        stop("Invalid type specified")
      }
    } else if (statistic == "z" || statistic == "n") {
      # Hamilton Table B.5 - Based on OLS Autoregressive Coefficient
      if (type == "nc") {
        table = cbind(
          c(-11.9,  -9.3,  -7.3,  -5.3, +1.01, +1.40, +1.79, +2.28),
          c(-12.9,  -9.9,  -7.7,  -5.5, +0.97, +1.35, +1.70, +2.16),
          c(-13.3, -10.2,  -7.9,  -5.6, +0.95, +1.31, +1.65, +2.09),
          c(-13.6, -10.3,  -8.0,  -5.7, +0.93, +1.28, +1.62, +2.04),
          c(-13.7, -10.4,  -8.0,  -5.7, +0.93, +1.28, +1.61, +2.04),
          c(-13.8, -10.5,  -8.1,  -5.7, +0.93, +1.28, +1.60, +2.03))
      } else if (type == "c") {
        table = cbind(
          c(-17.2, -14.6, -12.5, -10.2, -0.76, +0.01, +0.65, +1.40),
          c(-18.9, -15.7, -13.3, -10.7, -0.81, -0.07, +0.53, +1.22),
          c(-19.8, -16.3, -13.7, -11.0, -0.83, -0.10, +0.47, +1.14),
          c(-20.3, -16.6, -14.0, -11.2, -0.84, -0.12, +0.43, +1.09),
          c(-20.5, -16.8, -14.0, -11.2, -0.84, -0.13, +0.42, +1.06),
          c(-20.7, -16.9, -14.1, -11.3, -0.85, -0.13, +0.41, +1.04))
      } else if (type == "ct") {
        table = cbind(
          c(-22.5, -19.9, -17.9, -15.6, -3.66, -2.51, -1.53, -0.43),
          c(-25.7, -22.4, -19.8, -16.8, -3.71, -2.60, -1.66, -0.65),
          c(-27.4, -23.6, -20.7, -17.5, -3.74, -2.62, -1.73, -0.75),
          c(-28.4, -24.4, -21.3, -18.0, -3.75, -2.64, -1.78, -0.82),
          c(-28.9, -24.8, -21.5, -18.1, -3.76, -2.65, -1.78, -0.84),
          c(-29.5, -25.1, -21.8, -18.3, -3.77, -2.66, -1.79, -0.87))
      } else {
        stop("Invalid type specified")
      }
    } else {
      stop("Invalid statistic specified")
    }

    # Transpose:
    Table = t(table)
    colnames(Table) = c("0.010", "0.025", "0.050", "0.100", "0.900",
                        "0.950", "0.975", "0.990")
    rownames(Table) = c(" 25", " 50", "100", "250", "500", "Inf")
    ans = list(
      x = as.numeric(rownames(Table)),
      y = as.numeric(colnames(Table)),
      z = Table)
    class(ans) = "gridData"

    # Exclude Inf:
    if (!includeInf) {
      nX = length(ans$x)
      ans$x = ans$x[-nX]
      ans$z = ans$z[-nX, ]
    }

    # Add Control:
    attr(ans, "control") <-
      c(table = "adf", trend = trend, statistic = statistic)

    # Return Value:
    ans
  }
  
# m number of non-stationary common factors
panic.f <- function(Fpc, lambda, type = "const"){
  # number of factors
  r <- ncol(Fpc)
  # if r = 1
  if(r == 1){
    test <- adf(Fpc, type = type)
  }else{
    # statistic
    if(type == "const"){
      mt <- c(-Inf, -13.73, -23.535, -32.296, -40.442, -48.442, -57.040)
    }
    if(type == "trend"){
      mt <- c(-Inf, -21.313, -31.356, -40.180, -48.421, - 55.818, - 64.393)
    }
    MQcc <- Inf

    # demean estimator
    if(type == "const"){
      Fpc <- scale(Fpc, TRUE, FALSE)
    }
    if(type == "trend"){
      Fpc <- resid(lm(Fpc ~ c(1:nrow(Fpc))))
    }

    # starting
    m <- r

    while(m != 0 & abs(mt[m+1]) < abs(MQcc)){
      # beta c
      betac <- eigen((t(Fpc)%*%Fpc)/(T^2))$vectors[,1:m,drop=FALSE]
      # Yc
      Tn <- nrow(Fpc)
      Yc <- matrix(0, Tn, m)
      for(t in 1: Tn)
        Yc[t,] <- t(betac)%*%t(Fpc[t,,drop=FALSE])

      # procedure a
      J <- 4*ceiling(min(N,T)/100)^(1/4)
      # residuals i
      if(m > 1){
        epsilonc <- residuals(VAR(Yc, p = 1, type = "none"))
      }else{
        epsilonc <- as.matrix(residuals(arima(Yc, c(1, 0, 0)),
                      include.mean = FALSE))
      }
      Vepsilonc <- matrix(0, m, m)
      tt <- nrow(epsilonc)
      Kj <- 0
      for(j in 1 : J){
        # j
        Kj <- Kj + 1-j/(J+1)
        Vepsilonc <- Vepsilonc + Kj*((t(epsilonc[(j+1):tt,,drop=FALSE])%*%
                      epsilonc[1:(tt-j),,drop = FALSE])/T)
      }

      # ii
      A <- matrix(0, m, m)
      for(t in 2 : Tn){
        A <- A + t(Yc[t-1,,drop = FALSE])%*%Yc[t,,drop = FALSE] +
              t(Yc[t,,drop=FALSE])%*%Yc[t-1,,drop=FALSE]
      }
      # theta c
      thetac <- .5*(A - T*(Vepsilonc + t(Vepsilonc)))%*%solve(t(Yc[-T,,drop=FALSE
                  ])%*%Yc[-T,,drop=FALSE])
      # small eigenvalue
      vmc <- eigen(thetac)$values[m]
      # test
      MQcc <- T*(vmc - 1)
      # test
      if(abs(mt[m+1]) > abs(MQcc)){
        m <- m
      }else{
        m <- m - 1
      }
    }
  }
  return(m)
}

