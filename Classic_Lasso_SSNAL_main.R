Classic_Lasso_SSNAL_main <- function(A, b, lambda, parmain, y, xi, x){
  dscale <- parmain$dscale
  Ascaleyes <- parmain$Ascaleyes
  m <- parmain$m
  n <- parmain$n
  
  if ("A" %in% names(parmain)) { A <- parmain$A }
  tstart <- parmain$tstart
  #tic()
  existA <- parmain$existA
  Lip <- parmain$Lip
  scale <- parmain$scale
  maxiter <- parmain$maxiter
  printyes <- parmain$printyes
  rescale <- parmain$rescale
  stoptol <- parmain$stoptol
  orgojbconst <- parmain$orgojbconst
  
  sigmaLip <- 1/Lip
  if (norm(dscale -1, "2")<eps){ld <- lambda}
  else {ld <- lambda * dscale}
  lambdaorg <- lambda
  borg <- b
  normborg <- 1 + norm(borg,"2")
  Atxi <- eigenMatMult(t(A), xi, n_cores=4) #t(A)%*%xi
  Ax <- eigenMatMult(A, x, n_cores=4) #A%*%x
  
  obj1 <- 0.5 * norm(Ax - borg, "2")^2 + lambdaorg*norm(x) + orgojbconst
  obj2 <- -(0.5*(norm(xi,"2")^2)+t(borg) %*% xi) + orgojbconst
  
  # return(rescale)
  # return(c(obj1,obj2, rescale, scale))
  # return(scale)
  bscale <- 1
  cscale <- 1
  if (scale == 1){
    b <- b / sqrt(bscale * cscale)
    xi <- xi / sqrt(bscale * cscale)
    #Amap <- function(x) Amap0(x*sqrt(bscale/cscale))
    #ATmap <- function(x) ATmap0(x*sqrt(bscale/cscale))
    A<-A*sqrt(bscale/cscale)
    if (existA){ A<-A*sqrt(bscale/cscale)}
    lambda <- lambda / cscale
    ld <- ld/cscale
    x <- x/bscale
    y <- y/cscale
    Ax <- Ax / sqrt(bscale*cscale)
    Atxi <- Atxi / cscale
  }
  else {
    #Amap <- Amap0
    #ATmap <- ATmap0
  }
  #return(Ascaleyes)
  if (Ascaleyes){
    #Amap <- function(x) Amap(dscale * x)
    #ATmap <- function(x) dscale * ATmap(x)
    Atxi <- dscale * Atxi
    y <- dscale * y
  }
  
  normb <- 1 + norm(b,"2")
  #Ainput_nal = list(Amap = Amap,
  #                  ATmap = ATmap,
  #                  
  # #)
  # return(normb)
  #if (existA){Ainput_nal <- c(Ainput_nal, A=A)}
  sigma = max(1/sqrt(Lip), min(c(1,sigmaLip,lambdaorg)))
  #return(sigma)
  if (Ascaleyes){sigma <- 3}
  if ("sigma" %in% names(parmain)){sigma <- parmain$sigma}
  #return(c(obj1,obj2))
  Rp1 <- Ax- b
  Rp <- Rp1 + xi
  Rd <- Atxi + y
  primfeas <- norm(Rp,"2") / normb
  dualfeas <- norm(Rd,"2")/(1+norm(y,"2"))
  primfeasorg <- sqrt(bscale * cscale)*norm(Rp,"2")/normborg
  dualfeasorg <- norm(Rd / dscale, "2")*cscale / (1+norm(y/dscale,"2")*cscale)
  maxfeas <- max(primfeas, dualfeas)
  maxfeasorg <- max(primfeasorg, dualfeasorg)
  relgap <- (obj1 - obj2)/(1 + obj1 + obj2)
  runhist <- list(dualfeasorg1 = dualfeasorg,
                  primfeasorg1 = primfeasorg)
  
<<<<<<< Updated upstream
  return(c(primfeas,dualfeas,primfeasorg,
          dualfeasorg,
           maxfeas,maxfeasorg,relgap))
=======
  # return(c(primfeas,dualfeas,primfeasorg,
  #         dualfeasorg,
  #          maxfeas,maxfeasorg,relgap))
>>>>>>> Stashed changes
  
  #if (printyes) {
  #  printf('\n \t\t   Classic Lasso: SSNAL      ')
  #  printf('\n n = %3.0d, m = %3.0d',n, m)
  #  printf('\n bscale = %3.2d, cscale = %3.2d', bscale, cscale)
  #  printf('\n ---------------------------------------------------')
  #  printf('\n  iter|  [pinfeas  dinfeas]  [pinforg  dinforg]    relgaporg|    pobj          dobj    |')
  #  printf(' time | sigma |rankS |')
  #  printf('\n*****************************************************')
  #  printf('**************************************************************')
  #  printf('\n #%3.1d|  %3.2d %3.2d %3.2d %3.2d %- 3.2d %- 8.7d %- 8.7d  %5.1d', 0,primfeas,dualfeas,primfeasorg,dualfeasorg,relgap,obj1,obj2,
  #         toc())#etime(clock,tstart)); 
  #  printf('  %3.2d ',sigma)
  #}
  
  #SSNCG
  
  SSNCG <- 1
  if (SSNCG) {
    parNCG <- list(sigma = sigma,
                   tolconst = 0.5,
                   n = n)
  }
  maxitersub <- 10
  breakyes <- 0
  prim_win <- 0
  dual_win <- 0
  ssncgop <- list(existA=existA,
                  tol = stoptol,
                  precond = 0,
                  bscale = bscale,
                  cscale = cscale,
                  printsub = printyes
  )
  if (Ascaleyes) {
    ssncgop <- c(ssncgop, Ascaleyes=1)
    ssncgop$dscale = dscale
  }
  else {ssncgop <- c(ssncgop, Ascaleyes = 0)}
  
  sigmamax <- 1e7
  sigmamin <- 1e-4
  if (Ascaleyes){
    sigmamax <- sigmamax * mean(dscale)
  }
  
  # ADDED THESE LINES IN TO COMPLY WITH R NEED
  # FOR VARIABLES TO ACTUALLY EXIST BEFORE IF
  # STATEMENTS
  # NORMALLY COMPUTED IN THE COURSE OF THE FUNCTION
  # LATER FOR ITERATIONS mod3==1
  normy <- norm(y, "2")
  normAtxi <- norm(Atxi, "2")
  normx <- norm(x,"2")
  normyxi <- max(normx, normAtxi)
  
  for (iter in c(1:maxiter)) {
    if ((rescale==1 & maxfeas < 5e2 & iter%%3 ==1 & iter > 1)|(!is.null(A) & rescale >= 2 & maxfeas < 1e-1 & abs(relgap)<0.05 & iter >= 5 & max(normx/normyxi, normyxi/normx) > 1.7 & iter%%5 ==1)){
      normy <- norm(y, "2")
      normAtxi <- norm(Atxi, "2")
      normx <- norm(x,"2")
      normyxi <- max(normx, normAtxi)
      #return(c(sigma, normx, normyxi, bscale, cscale))
      c(sigma, bscale2, cscale2, sbc, sboc, bscale, cscale) <- mexscale(sigma, normx, normyxi, bscale, cscale)
      
      Amap <- function(x) Amap(x*sboc)
      ATmap <- function(x) ATmap(x*sboc)
      Ainput_nal <- list(Amap = Amap,
                         ATmap = ATmap)
      if (exist(A)){Ainput_nal <- c(Ainput_nal, A= A*sboc)}
      b <- b/sbc
      ld <- ld / cscale2
      x <- x/bscale2
      xi <- xi/sbc
      Atxi <- Atxi/cscale2
      Ax <- Ax/sbc
      ssncgop$bscale <- bscale
      ssncgop$cscale <- cscale
      normb <- 1 + norm(b,"2")
      if (printyes) {
        printf('\n    ');
        printf('[rescale=%1.0d: %2.0d| %3.2d %3.2d %3.2d | %3.2d %3.2d| %3.2d]',rescale,iter,normx,normAtxi, normy,bscale,cscale,sigma)
      }
      rescale <- rescale + 1
      prim_win <- 0
      dual_win <- 0
    }
    #return("Got here")
    xold <- x
    parNCG$sigma <- sigma
    if (dualfeas < 1e-5){maxitersub <- max(maxitersub, 30)}
    else if(dualfeas > 1e-3){maxitersub <- max(maxitersub, 30)}
    else if(dualfeas > 1e-1){maxitersub <- max(maxitersub, 20)}
    ssncgop <- c(ssncgop, maxitersub = maxitersub)
    
    return(c(n,mean(b),mean(x),
            mean(Ax),mean(Atxi),mean(xi),mean(ld)))
    
    #return(ssncgop)
    
    #c(y,Atxi, xi, parNCG, runhist_NCG, info_NCG) <- classic_Lasso_SSNCG(n,b,Ainput_nal, x, Ax, Atxi, xi, ld, parNCG, ssncgop)
    
    # return(Classic_Lasso_SSNCG(n,b,Ainput_nal, x, Ax, Atxi, xi, ld, parNCG, ssncgop))
    
    if (info_ncg$breakyes < 0){
      parNCG$tolconst <- max(parNCG$tolconst/1.06, 1e-3)
    }
    Rd <- Atxi + y
    X <- -sigma * info_NCG$ytmp
    Ax <- info_NCG$Ax
    
    normRd <- norm(Rd, "2")
    normy <- norm(y, "2")
    dualfeas <- normRd / (1+normy)
    if (Ascaleyes){
      dualfeasorg <- norm(Rd / dscale, "2")*cscale / (1+norm(y/dscale,"2")*cscale)
    }
    else {
      dualfeasorg <- normRd * cscale / (1+normy * cscale)
    }
    Rp1 <- Ax - b
    Rp <- Rp1 + xi
    normRp <- norm(Rp, "2")
    primfeas <- normRp / normb
    primfeasorg <- sqrt(bscale * cscale) * normRp / normborg
    maxfeas <- max(primfeas, dualfeas)
    maxfeasorg <- max(primfeasorg, dualfeasorg)
    runhist <- list( dualfeas = dualfeas,
                     primfeas = primfeas,
                     maxfeas = maxfeas,
                     primfeasorg = primfeasorg,
                     dualfeasorg = dualfeasorg,
                     maxfeasorg = maxfeasorg,
                     sigma = sigma,
                     rankS = sum(1 - parNCG$rr),
                     cnt_Amap = info_NCG$cnt_Amap,
                     cnt_ATmap = info_NCG$cnt_ATmap,
                     cnt_pAATmap = info_NCG$cnt_pAATmap,
                     cnt_fAATmap = info_NCG$cnt_fAATmap)
    if (max(primfeasorg, dualfeasorg) < 500 * max(1e-6, stoptol)) {
      grad <- ATmap0(Rp1 * sqrt(bscale * cscale))
      if (Ascaleyes) {
        etaorg <- norm(grad + projinf(dscale * x %*% bscale - grad, lambdaorg), "2")
        eta <- etaord / (1 + norm(grad,"2")+norm(x%*%bscale))
      }
      else {
        etaorg <- norm(grad + projinf(x %*% bscale -grad, lambdaorg))
        eta <- etaorg / (1 + norm(grad, "2") + norm(x%*% bscale))
      }
      if (eta < stoptol) {
        breakyes <- 1
        msg <- 'converged'
      }
    }
    if (printyes) {
      objscale <- bscale * cscale
      primobj <- objscale * (0.5 * norm(Rp1, "2")^2 + norm(ld*x) + orgojbconst)
      dualobj <- objscale * (-0.5 * norm(xi)^2 + t(b)%*% xi) + orgojbconst
      relgap <- (primobj - dualobj)/(1+abs(primobj) + abs(dualobj))
      # ttime <- measure time
    }
    printf('\n %5.0d| [%3.2e %3.2e] [%3.2e %3.2e]  %- 3.2e| %- 5.4e %- 5.4e |',iter,primfeas,dualfeas,primfeasorg,dualfeasorg,relgap,primobj,dualobj)
    printf(' %5.1f| %3.2e|sigamorg = %3.2e  |',ttime, sigma, sigma*bscale/cscale)
    if (iter >= 1){
      printf('%3.0d|',sum(1- parNCG$rr))
    }
    if (exists(eta)){
      printf('\n [eta = %3.2e, etaorg = %3.2e]',eta, etaorg)
    }
    if (iter %% 3 == 1){
      normx <- norm(x,"2")
      normAtxi <- norm(Atxi,"2")
      normy <- norm(y,"2")
      if (printyes) {
        printf('\n        [normx,Atxi,y =%3.2e %3.2e %3.2e]',normx,normAtxi,normy)
      }
    }
    runhist <- c(runhist, primobj = primobj, dualobj = dualobj, time = ttime, relgap = relgap)
    if (breakyes > 0) {
      if (printyes){
        printf('\n  breakyes = %3.1f, %s',breakyes,msg)
      }
    }
    if(primfeasorg < dualfeasorg) {
      prim_win <- prim_win + 1
    }
    else {
      dual_win <- dual_win + 1
    }
    c(sigma, prim_win, dual_win) <- mexsigma_update_classic_Lasso_SSNAL(sigma, sigmamax, sigmamin, prim_win, dual_win, iter, info_NCG$breakyes)
  }
  if (!printyes) {
    ttime <- toc() # end time
  }
  if ((iter == maxiter) & (breakyes == 0)){
    msg <- "maximum number of iterations reached"
  }
  
  info <- list(iter = iter,
               bscale = bscale,
               cscale = cscale,
               Ax = Ax,
               ttime = ttime,
               msg = msg)
}