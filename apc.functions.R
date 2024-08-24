# apc_functions.R
# functions for age, period, cohort analysis
# using the intristic estimator
#
require(MASS)
#
#average.apc.R
average.apc <- function(df1,df2) {
  b <- (df1$estimate + df2$estimate)/2
  bL <- (df1$lower + df2$lower)/2
  bU <- (df1$upper + df2$upper)/2
  x <- df1$x
  apcfact <- df1$apcfact
  return(data.frame(y=b,x=x, apcfact=apcfact, lower=bL, upper=bU ))
}

checknames.apc <- function(df) {
  check_names <- c("N", "Y", "period", "age", "cohort")
  errcount <- 0
  for (name in check_names) {
    if (!name %in% colnames(df)) {
      errcount <- errcount + 1
      cat("variable ", name, " not found in data", "\n")
    }
  }
  if (errcount > 0) {
    print("required variable names are: Y, N, age, period, cohort")
    stop("check variable names")
  }
}

# design.apc.R
design.apc <- function(df, ref="last") {
  # produce design matrix for age, period, cohort analysis
  # assume df with colnames: Y, N, age, period, cohort 
  # process to construct effect coded design matrix
  # options: ref = "last" (default)
  #              = "first" 
  # handles 1st and last factor levels as reference
  
  n <- length(df$Y)
  o <- rep(1,n)
  age    <- df$age
  period <- df$period
  cohort <- df$cohort
  #
  # construct centered design matrix
  #
  # convert apc factors to (0/1) dummy-variable coding
  #
  plab <- unique(period)
  P <- array(NA, c(n,length(plab)))
  for (i in 1:length(plab)) {
    P[,i] <- as.numeric(period==plab[i])
  }
  
  alab <- unique(age)
  A <- array(NA, c(n,length(alab)))
  for (i in 1:length(alab)) {
    A[,i] <- as.numeric(age==alab[i])
  }
  
  clab <- sort(unique(cohort))
  C <- array(NA,c(n, length(clab)))
  for (i in 1:length(clab)) {   
    C[,i] <- as.numeric(cohort==clab[i])
  }
  
  # A,P,C,alab,plab,clab used below
  
  # normalization switch default to LAST if empty
  
  if (ref=="last") {
    # 
    # convert dummies to centered effects (ANOVA) [LAST category as reference] Default
    #
    dimA <- length(alab)
    dimP <- length(plab)
    dimC <- length(clab)
    
    rA <- A[,dimA]
    rP <- P[,dimP]
    rC <- C[,dimC]
    
    cA.L <- A[,1:(dimA-1)] - rA
    cP.L <- P[,1:(dimP-1)] - rP
    cC.L <- C[,1:(dimC-1)] - rC
    
    fname.last <- function(xlab,name) {
      vname <- rep(NA, length(xlab) - 1)
      for(i in 1:length(vname)) {
        #function name for each column
        vname[i] <- paste(name,"_",unique(xlab[i]),sep='')
      }
      return(vname)
    }
    
    # apply names function
    
    fname.last(alab, "age")    -> agename
    fname.last(plab, "period") -> periodname
    fname.last(clab, "cohort") -> cohortname
    # X.str <- as.matrix(data.frame(o, A.str, P.str, C.str))
    X.L   <- as.matrix(data.frame(o,cA.L,cP.L,cC.L))
    colnames(X.L) <- c("Intercept",agename, periodname, cohortname)
    # 
    return(list(X=X.L, alab=alab, plab=plab, clab=clab, ref="last"))
  }
  
  #
  else if (ref=="first") {  
    #
    # centered effects (ANOVA) [using FIRST category as reference]
    #
    dimA <- length(alab)
    dimP <- length(plab)
    dimC <- length(clab)
    
    rA <- A[,1]
    rP <- P[,1] 
    rC <- C[,1] 
    
    cA.F <- A[,2:dimA] - rA
    cP.F <- P[,2:dimP] - rP
    cC.F <- C[,2:dimC] - rC
    
    # apc coef names function
    fname.first <- function(xlab,name) {
      vname <- rep(NA, length(xlab)-1)
      for(i in 1:length(vname)) {
        #function name for each column
        vname[i] <- paste(name,"_",xlab[i+1],sep='')
      }
      return(vname)
    }
    fname.first(alab,"age")    -> agename
    fname.first(plab,"period") -> periodname
    fname.first(clab,"cohort") -> cohortname
    
    # 
    X.F   <- as.matrix(data.frame(o,cA.F,cP.F,cC.F))
    colnames(X.F) <- c("Intercept",agename, periodname, cohortname)
    return(list(X=X.F, alab=alab, plab=plab, clab=clab, ref="first"))
  }
}  


# IE_apc.R
IE_linear <- function(df, ref="last", out="raw") {
  # ref is one of "first" or "last" (default)
  # out NULL=normalized output raw=unormalized output
  # df must contain variables named: age, period, cohort, Y, N
  checknames.apc(df)
  dX <- design.apc(df, ref)
  X     <- dX$X
  if (!is.matrix(X)) X <- as.matrix(X)
  ref   <- dX$ref
  alab  <- dX$alab
  plab  <- dX$plab
  clab  <- dX$clab
  Y   <- df$Y
  N   <- df$N
  XpX <- t(X)%*%X
  
  # empirical log rates
  y <- log(Y/N)
  n <- length(y)
  b <- ginv(XpX)%*%t(X)%*%y
  s2.e  <- sum( (y - X%*%b)^2)/(n - ncol(X) + 1)
  v.b   <- s2.e * ginv(XpX)
  s.b  <- sqrt(diag(v.b))
  rownames(b)     <- colnames(X)  
  rownames(v.b)   <- colnames(X)
  colnames(v.b)   <- colnames(X)
  results <- data.frame(b = b, 
                        se = s.b, 
                        Z = b/s.b,
                        p.val = 2*pnorm(-abs(b/s.b)))
  
  rlist <- list(estimate = b, std.error = s.b,
                vcov = v.b, 
                deviance = s2.e, results  = results, alab = alab, 
                plab = plab, clab = clab, ref = ref)
  
  olist <-  list(estimate = b, std.error = s.b,
                 vcov = v.b, 
                 deviance = s2.e, results=results)
  if (is.null(out)) {
    outp <- normparam.apc(rlist)
    return(outp) }
  else {
    return(olist)
    stop()
  }
}



###################################################
# glm APC uses iterative algorithm (Newton Raphson) 
###################################################

# Poisson Regression: mu = Nexp(XB)

newt.raphP <- function(b,D,R,X) { #ML
  mu <- R*exp(X%*%b)
  g  <- t(X)%*%(D - mu)
  H  <- t(X)%*%(c(mu)*X) # this is negative of hessian (which is neg def) so this is pos def
  b.new <- b + ginv(H)%*%g
  v.b   <- ginv(H)
  s.b   <- sqrt(diag(v.b))
  out   <- list(b=b.new, se=s.b, var.b=v.b)
  return(out)
}      

# Logit Regression: mu = Nexp(XB)/(1 + exp(XB))

newt.raphL <- function(b,D,R,X) { #ML
  mu <- R*exp(X%*%b)/(1 + exp(X%*%b))
  v <- mu/(1 + exp(X%*%b))
  g  <- t(X)%*%(D - mu)
  H  <- t(X)%*%(c(v)*X) # this is negative of hessian (which is neg def) so this is pos def
  b.new <- b + ginv(H)%*%g
  v.b   <- ginv(H)
  s.b   <- sqrt(diag(v.b))
  out   <- list(b=b.new, se=s.b, var.b=v.b)
  return(out) 
}      


# fit rate model
# df must contain variables named: age, period, cohort, Y, N
# ref is one of "first" or "last" (default)
# out NULL=normalized output raw=unormalized output
# bstart is optional arg for start values
IE_rate <- function(df,
                    bstart=NULL, 
                    family="Poisson",
                    ref="last",
                    out=NULL) {
  checknames.apc(df)
  dX <- design.apc(df, ref)
  X   <- dX$X
  ref   <- dX$ref
  alab  <- dX$alab
  plab  <- dX$plab
  clab  <- dX$clab
  Y   <- df$Y
  N   <- df$N
  # other useful returns
  
  b.old <- bstart
  if (is.null(bstart)) b.old  <- rep(0,dim(X)[2])
  if (is.null(family)) family <- "Poisson"  # default
  
  iter  <- 1
  db  <- 1
  tol <- 1.e-8
  if (family == "Poisson") {
    while (db > tol) {
      cat("Iteration = ",iter,"\n")
      out.A   <- newt.raphP(b.old,Y,N,X)
      b.new  <- out.A$b
      s.b    <- out.A$se
      v.b    <- out.A$var.b
      db     <- mean(abs((b.old - b.new)/b.old))
      mu.hat <- N*exp(X%*%b.new)
      dev    <- 2 * sum(Y * log(Y/mu.hat) - Y + mu.hat)
      cat("Deviance = ",dev, "\n")
      b.old <- b.new
      iter <- iter + 1
    }
    b.new <- as.matrix(b.new)
    rownames(b.new) <- colnames(X)  
    rownames(v.b)   <- colnames(X)
    colnames(v.b)   <- colnames(X)
    results <- data.frame(b = b.new, 
                          se = s.b, 
                          Z = b.new/s.b,
                          p.val = 2*pnorm(-abs(b.new/s.b)))
    
    rlist <- list(estimate = b.new, std.error = s.b,
                  vcov = v.b, 
                  deviance = dev, results = results, alab = alab, 
                  plab = plab, clab = clab, ref = ref)
    olist <-  list(estimate = b.new, std.error = s.b,
                   vcov = v.b, 
                   deviance = dev, results = results)
    
    if (is.null(out)) {
      outp <- normparam.apc(rlist)
      return(outp) }
    else {
      return(olist)
      stop()
    }
  }
  else if (family=="Binomial") {
    while (db > tol) {
      cat("Iteration=",iter,"\n")
      out.A  <- newt.raphL(b.old,Y,N,X)
      b.new  <- out.A$b
      s.b    <- out.A$se
      v.b    <- out.A$var.b
      db     <- mean(abs((b.old - b.new)/b.old))
      mu.hat <- N*exp(X%*%b.new)
      dev    <- 2 * sum(Y * log(Y/mu.hat) - Y + mu.hat)
      cat("Deviance = ", dev, "\n")
      b.old <- b.new
      iter  <- iter + 1
    }
    b.new <- as.matrix(b.new)
    rownames(b.new) <- colnames(X)
    rownames(v.b)   <- colnames(X)
    colnames(v.b)   <- colnames(X)
    
    results <- data.frame(b = b.new, 
                          se = s.b, 
                          Z = b.new/s.b,
                          p.val = 2*pnorm(-abs(b.new/s.b)))
    
    rlist <- list(estimate = b.new, std.error = s.b,
                  vcov = v.b, 
                  deviance = dev, results  = results, alab = alab, 
                  plab = plab, clab = clab, ref = ref)
    olist <-  list(estimate = b.new, std.error = s.b,
                   vcov = v.b, 
                   deviance = dev, results=results)
    
    if (is.null(out)) {
      outp <- normparam.apc(rlist)
      return(outp) }
    else {
      return(olist)
      stop()
    }
  }
}    

# # done with ML

#normparam.apc.R
# normalize coefficiens for all
# factor levels
normparam.apc <- function(mod) {
  
  alab <- mod$alab
  plab <- mod$plab
  clab <- mod$clab
  ref  <- mod$ref
  b    <- mod$estimate
  v    <- mod$vcov
  # get ref cat estimates and variances
  if (is.na(ref) | ref=="last") {
    b.refA    <- -sum(b[2:length(alab)])
    v.refA    <-  sum(as.vector(v[2:length(alab),2:length(alab)]))
    b.A       <-  c(b[2:length(alab)],b.refA)
    se.bA     <-  c(sqrt(diag(v[2:length(alab),2:length(alab)])), sqrt(v.refA)) 
    p.idx     <-  length(alab) + 1            # starting index (location) for period effects 
    p.maxi    <-  p.idx + length(plab) - 2 # ending index (location) for period effects
    b.refP    <- -sum(b[p.idx: p.maxi])
    v.refP    <-  sum(as.vector(v[p.idx:p.maxi,p.idx:p.maxi]))
    b.P       <-  c(b[p.idx:p.maxi],b.refP)
    se.bP     <-  c(sqrt(diag(v[p.idx:p.maxi,p.idx:p.maxi])), sqrt(v.refP))
    c.idx     <-  length(alab) + length(plab)  # starting index (location) for cohort effects 
    c.maxi    <-  length(b)                    # ending index (location) for cohort effects
    b.refC    <- -sum(b[c.idx: c.maxi])
    v.refC    <-  sum(as.vector(v[c.idx:c.maxi,c.idx:c.maxi]))
    b.C       <-  c(b[c.idx:c.maxi],b.refC)
    se.bC     <-  c(sqrt(diag(v[c.idx:c.maxi,c.idx:c.maxi])), sqrt(v.refC))
    # assemble new coef and se vector 
    b.norm    <-  as.vector(c(b[1], b.A, b.P, b.C))
    seb.norm  <-  c(sqrt(v[1,1]), se.bA, se.bP, se.bC)
    # revised names:
    # function to generate coef names
    fname.L  <- function(xlab,name) {
      vname <- NULL
      for(i in 1:length(xlab)) {
        #function name for each coef
        vname[i] <- paste(name,"_",xlab[i],sep='')
      }
      return(vname)
    }
    
    fname.L(alab,"age")    -> agename
    fname.L(plab,"period") -> periodname
    fname.L(clab,"cohort") -> cohortname
    
    names(b.norm) <- names(seb.norm) <- c("Intercept",agename, periodname, cohortname)
    
    Z <- b.norm/seb.norm
    p.Z <- 2 * (1-pnorm(abs(Z)))
    b.025 <- b.norm - 1.96 * (seb.norm)
    b.975 <- b.norm + 1.96 * (seb.norm)
    #print(data.frame(b.norm,seb.norm,Z,p.Z, b.025, b.975), digits=3)
    results <- data.frame(estimate=b.norm,
                          std.error=seb.norm,
                          Z=Z,
                          p.val=p.Z, 
                          lower = b.025, 
                          upper = b.975)
    return(results)
  }
  else if (ref=="first") {
    dimA      <- length(alab)
    a.idx     <- 2
    a.maxi    <- dimA  
    b.refA    <- -sum(b[a.idx:a.maxi]) # 
    v.refA    <-  sum(as.vector(v[a.idx:a.maxi,a.idx:a.maxi]))
    b.A       <-  c(b.refA, b[a.idx:a.maxi])
    se.bA     <-  c(sqrt(v.refA), sqrt(diag(v[a.idx:a.maxi,a.idx:a.maxi])))
    p.idx     <-  length(alab) + 1         # starting index (location) for period effects 
    p.maxi    <-  p.idx + length(plab) - 2 # ending index (location) for period effects
    b.refP    <- -sum(b[p.idx: p.maxi])
    v.refP    <-  sum(as.vector(v[p.idx:p.maxi,p.idx:p.maxi]))
    b.P       <-  c(b.refP, b[p.idx:p.maxi])
    se.bP     <-  c(sqrt(v.refP), sqrt(diag(v[p.idx:p.maxi, p.idx:p.maxi])))
    c.idx     <-  length(alab) + length(plab)  # starting index (location) for cohort effects 
    c.maxi    <-  length(b)                                    # ending index (location) for cohort effects
    b.refC    <- -sum(b[c.idx: c.maxi])
    v.refC    <-  sum(as.vector(v[c.idx:c.maxi,c.idx:c.maxi]))
    b.C       <-  c(b.refC, b[c.idx:c.maxi])
    se.bC     <-  c(sqrt(v.refC), sqrt(diag(v[c.idx:c.maxi,c.idx:c.maxi])))
    # assemble new coef and se vector 
    b.norm    <-  as.vector(c(b[1], b.A, b.P, b.C))
    seb.norm  <-  c(sqrt(v[1,1]), se.bA, se.bP, se.bC)
    # revised names:
    # function to generate coef names (not used)
    fname.L  <- function(xlab,name) {
      vname <- NULL
      for(i in 1:length(xlab)) {
        #function name for each coef
        vname[i] <- paste(name,"_", xlab[i],sep='')
      }
      return(vname)
    }
    
    fname.L(alab,"age")    -> agename
    fname.L(plab,"period") -> periodname
    fname.L(clab,"cohort") -> cohortname
    
    names(b.norm) <- names(seb.norm) <- c("Intercept",agename, periodname, cohortname)
    
    Z     <- b.norm/seb.norm
    p.Z   <- 2 * (1-pnorm(abs(Z)))
    b.025 <- b.norm - 1.96 * (seb.norm)
    b.975 <- b.norm + 1.96 * (seb.norm)
    
    results <- data.frame(estimate=b.norm,
                          std.error=seb.norm,
                          Z=Z,
                          p.val=p.Z, 
                          lower = b.025, 
                          upper = b.975)
    return(results)
  }
  
}

