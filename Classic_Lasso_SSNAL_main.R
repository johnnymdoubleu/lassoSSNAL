Classic_Lasso_SSNAL_main <- function(A, orig_A, b, lambda, parmain, y, xi, x){
  dscale <- parmain$dscale
  Ascaleyes <- parmain$Ascaleyes
  m <- parmain$m
  n <- parmain$n
  #orig_A <- A
  
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
  # Atxi <- t(A)%*%xi
  Atxi <- eigenMapMatMult(t(A), xi, 4)
  # Ax <- A%*%x
  Ax <- eigenMapMatMult(A, x, 4)
  
  obj1 <- 0.5 * norm(Ax - borg, "2")^2 + lambdaorg*norm(x) + orgojbconst
  obj2 <- -(0.5*(norm(xi,"2")^2)+t(borg) %*% xi) + orgojbconst
  
  #return(rescale)
  #return(c(obj1,obj2))
  #return(scale)
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
  #)
  
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
  
  #return(c(primfeas,dualfeas,primfeasorg,
  #         dualfeasorg,
  #          maxfeas,maxfeasorg,relgap))
  
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
    if ((rescale==1 & maxfeas < 5e2 & iter%%3 ==1 & iter > 1)
        |(is.null(A) & (rescale >= 2 & maxfeas < 1e-1 & abs(relgap)<0.05 & iter >= 5 & max(normx/normyxi, normyxi/normx) > 1.7 & iter%%5 ==1))){
      normy <- norm(y, "2")
      normAtxi <- norm(Atxi, "2")
      normx <- norm(x,"2")
      normyxi <- max(normx, normAtxi)
      #return(c(sigma, normx, normyxi, bscale, cscale))
      
      cat("normy=",normy,"\n")
      cat("normAtxi=",normAtxi,"\n")
      cat("normx=",normx,"\n")
      cat("normyxi=",normyxi,"\n")
      
      mso <- mexscale(sigma, normx, normyxi, bscale, cscale)
      sigma <- mso[1]
      bscale2 <- mso[2]
      cscale2 <- mso[3]
      sbc <- mso[4]
      sboc <- mso[5]
      bscale <- mso[6]
      cscale <- mso[7]
      
      #c(sigma, bscale2, cscale2, sbc, sboc, bscale, cscale) <- mexscale(sigma, normx, normyxi, bscale, cscale)
      
      #Amap <- function(x) Amap(x*sboc)
      #ATmap <- function(x) ATmap(x*sboc)
      #Ainput_nal <- list(Amap = Amap,
      #                   ATmap = ATmap)
      #if (exist(A)){Ainput_nal <- c(Ainput_nal, A= A*sboc)}
      
      A <- A*sboc
      
      b <- b/sbc
      ld <- ld / cscale2
      x <- x/bscale2
      xi <- xi/sbc
      Atxi <- Atxi/cscale2
      Ax <- Ax/sbc
      ssncgop$bscale <- bscale
      ssncgop$cscale <- cscale
      normb <- 1 + norm(b,"2")
      #if (printyes) {
      #  printf('\n    ');
      #  printf('[rescale=%1.0d: %2.0d| %3.2d %3.2d %3.2d | %3.2d %3.2d| %3.2d]',rescale,iter,normx,normAtxi, normy,bscale,cscale,sigma)
      #}
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
    
    #return(c(n,mean(b),mean(x),
    #         mean(Ax),mean(Atxi),mean(xi),mean(ld)))
    
    #return(ssncgop)
    
    #c(y,Atxi, xi, parNCG, runhist_NCG, info_NCG) <- classic_Lasso_SSNCG(n,b,Ainput_nal, x, Ax, Atxi, xi, ld, parNCG, ssncgop)
    
    #return(Classic_Lasso_SSNCG(n,b,A, x, Ax, Atxi, xi, ld, parNCG, ssncgop))
    cat("number of iteration=",iter,"\n")
    clssncgout <- Classic_Lasso_SSNCG(n,b,A, x, Ax, Atxi, xi, ld, parNCG, ssncgop)
    y <- clssncgout$y
    Atxi <- clssncgout$Atxi
    xi <- clssncgout$xi
    parNCG <- clssncgout$par
    runhist_NCG <- clssncgout$runhist
    info_NCG <- clssncgout$info
    
    print(sum(y))
    print(sum(Atxi))
    print(sum(xi))
    print(info_NCG$breakyes)
    
    if (info_NCG$breakyes < 0){
      parNCG$tolconst <- max(parNCG$tolconst/1.06, 1e-3)
    }
    
    cat("parncg_tolconst=",parNCG$tolconst,"\n")
    #return()
    
    Rd <- Atxi + y
    x <- -sigma * info_NCG$ytmp
    Ax <- info_NCG$Ax
    
    normRd <- norm(Rd, "2")
    normy <- norm(y, "2")
    
    cat("normRd=",normRd,"\n")
    cat("normy=",normy,"\n")
    #return()
    
    dualfeas <- normRd / (1+normy)
    
    cat("dualfeas=",dualfeas,"\n")
    #return()
    
    if (Ascaleyes){
      dualfeasorg <- norm(Rd / dscale, "2")*cscale / (1+norm(y/dscale,"2")*cscale)
    }
    else {
      dualfeasorg <- normRd * cscale / (1+normy * cscale)
    }
    
    cat("dualfeasorg=",dualfeasorg,"\n")
    #return()
    
    Rp1 <- Ax - b
    Rp <- Rp1 + xi
    normRp <- norm(Rp, "2")
    primfeas <- normRp / normb
    primfeasorg <- sqrt(bscale * cscale) * normRp / normborg
    maxfeas <- max(primfeas, dualfeas)
    maxfeasorg <- max(primfeasorg, dualfeasorg)
    
    cat("normRp=",normRp,"\n")
    cat("primfeas=",primfeas,"\n")
    cat("primfeasorg=",primfeasorg,"\n")
    cat("maxfeas=",maxfeas,"\n")
    cat("maxfeasorg=",maxfeasorg,"\n")
    cat("printyes=",printyes,"\n")
    #return()
    
    
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
    
    
    
    # IMPORTANT!
    # NOTE:
    #   This is the code for which differences lead to differences in
    #   relgap/eta in the R implementation.
    #   I believe cascading numerical imprecisions throughout the algorithm
    #   prohibit doing any better, and appears intractable.
    
    
    
    
    
    
    
    
    
    if (max(primfeasorg, dualfeasorg) < 500 * max(1e-6, stoptol)) {
      #grad <- ATmap0(Rp1 * sqrt(bscale * cscale))
      # grad <- t(orig_A) %*% (Rp1 * sqrt(bscale*cscale))
      grad <- eigenMapMatMult(t(orig_A), (Rp1 * sqrt(bscale * cscale)), 4)
      
      # cat("rtbc=",sqrt(bscale*cscale),"\n")
      cat("Rp1norm=",norm(Rp1,"2"),"\n")
      cat("dscale=",norm(dscale,"2"),"\n")
      cat("lambdaorg=",lambdaorg,"\n")
      # cat("innorm=",norm(dscale * x * bscale - grad,"2"),"\n")
      # cat("outnorm_y=",norm(proj_inf(dscale * x * bscale - grad, lambdaorg)$y,"2"),"\n")
      # cat("outnorm_rr=",norm(proj_inf(dscale * x * bscale - grad, lambdaorg)$rr,"2"),"\n")
      print(norm(grad,"2"))
      
      if (Ascaleyes) {
        # etaorg <- norm(grad + proj_inf((dscale * x) %*% bscale - grad, lambdaorg)$y, "2")
        # eta <- etaorg / (1 + norm(grad, "2")+norm((dscale * x) %*% bscale, "2"))
        etaorg <- norm(grad + proj_inf(eigenMapMatMult((dscale*x), bscale, 4) - grad, lambdaorg)$y, "2")
        eta <- etaorg / (1 + norm(grad, "2")+norm(eigenMapMatMult((dscale*x), bscale, 4), "2"))
        print("######")
        print("Nearing the end")
        print("######")
        # cat("etaorg=",etaorg,"\n")
        # cat("etaorg=",eta,"\n")
        
      }
      else {
        # etaorg <- norm(grad + projinf(x %*% bscale - grad, lambdaorg, "2"))
        # eta <- etaorg / (1 + norm(grad, "2") + norm(x %*% bscale, "2"))
        etaorg <- norm(grad + projinf(eigenMapMatMult(x, bscale,4) - grad, lambdaorg, "2"))
        eta <- etaorg / (1 + norm(grad, "2") + norm(eigenMapMatMult(x, bscale,4), "2"))
      }
      # cat("number of iteration = ", iter, "\n")
      cat("etaorg=",etaorg,"\n")
      cat("eta=",eta,"\n")
      
      if (eta < stoptol) {
        breakyes <- 1
        msg <- 'converged'
        cat(msg, "\n")
      }
    }
    
    
    if (printyes) {
      objscale <- bscale * cscale
      primobj <- objscale * (0.5 * norm(Rp1,"2")^2 + norm(ld * x,"1")) + orgojbconst
      dualobj <- objscale * (-0.5 * norm(xi,"2")^2 + t(b) %*% xi) + orgojbconst
      # dualobj <- objscale * (-0.5 * norm(xi,"2")^2 + eigenMapMatMult(t(b),xi,4)) + orgojbconst
      relgap <- (primobj - dualobj)/(1+ abs(primobj) + abs(dualobj))
      # ttime <- measure time
    }
    
    cat("primobj=",primobj,"\n")
    cat("dualobj=",dualobj,"\n")
    cat("relgap=",relgap,"\n")
    #return()
    
    #printf('\n %5.0d| [%3.2e %3.2e] [%3.2e %3.2e]  %- 3.2e| %- 5.4e %- 5.4e |',iter,primfeas,dualfeas,primfeasorg,dualfeasorg,relgap,primobj,dualobj)
    #printf(' %5.1f| %3.2e|sigamorg = %3.2e  |',ttime, sigma, sigma*bscale/cscale)
    if (iter >= 1){
      #printf('%3.0d|',sum(1- parNCG$rr))
    }
    #if (exists(eta)){
    #  printf('\n [eta = %3.2e, etaorg = %3.2e]',eta, etaorg)
    #}
    if (iter %% 3 == 1){
      normx <- norm(x,"2")
      normAtxi <- norm(Atxi,"2")
      normy <- norm(y,"2")
      #if (printyes) {
      #  printf('\n        [normx,Atxi,y =%3.2e %3.2e %3.2e]',normx,normAtxi,normy)
      #}
    }
    
    #runhist <- c(runhist, 
    #             primobj[iter] = primobj, 
    #             dualobj[iter] = dualobj,
    #             time[iter] = ttime, 
    #             relgap[iter] = relgap)
    
    #return(runhist)
    
    cat("breakyes=",breakyes,"\n")
    #return()
    
    if (breakyes > 0) {
      #if (printyes){
      #  printf('\n  breakyes = %3.1f, %s',breakyes,msg)
      #}
      break
    }
    
    if(primfeasorg < dualfeasorg) {
      prim_win <- prim_win + 1
    }
    else {
      dual_win <- dual_win + 1
    }
    
    cat("prim_win=",prim_win,"\n")
    cat("dual_win=",dual_win,"\n")
    
    msuout <- mexsigma_update(sigma, sigmamax, sigmamin, prim_win, dual_win, iter, info_NCG$breakyes)
    
    sigma <- msuout[1]
    prim_win <- msuout[2]
    dual_win <- msuout[3]
    #       !!!! IMPORTANT !!!!!
      
    #  Comment/uncomment the following "Stop" input command
    #   in both R + MATLAB if you want to manually debug
    #   or supervise the algorithm.
    
    #dbgtest <- readline("Stop")
    #if(dbgtest=="die") break
    
    #c(sigma, prim_win, dual_win) <- mexsigma_update_classic_Lasso_SSNAL(sigma, sigmamax, sigmamin, prim_win, dual_win, iter, info_NCG$breakyes)
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
               msg = msg)
  
  output <- list(obj=c(obj1,obj2),
                 y=y,
                 xi=xi,
                 x=x,
                 info=info,
                 runhist=runhist)
  return(output)
}