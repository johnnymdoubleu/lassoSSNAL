simplified_SSNCG <- function(n, b, A, x0, Ax0, Atxi0, xi0, ld, par, options) {
  
  printsub <- 0
  breakyes <- 0
  maxitersub <- 50
  tiny <- 1e-10
  tol <- 1e-6
  maxitpsqmr <- 500
  precond <- 0
  Ascaleyes <- 0
  
  if ("printsub" %in% names(options)) printsub <- options$printsub
  if ("maxitersub" %in% names(options)) maxitersub <- options$maxitersub
  if ("tiny" %in% names(options)) tiny <- options$tiny
  if ("tol" %in% names(options)) tol <- options$tol
  if ("precond" %in% names(options)) precond <- options$precond
  if ("Ascaleyes" %in% names(options)) Ascaleyes <- options$Ascaleyes
  
  existA <- options$existA
  sig <- par$sigma
  bscale <- options$bscale
  cscale <- options$cscale
  normborg <- 1 + norm(b, "2") * sqrt(bscale * cscale)
  
  # CHANGE ALL OCCURRENCES OF Amap and ATmap to raw calls
  
  yinput <- -Atxi0 - x0/sig
  
  pi_out <- proj_inf(yinput, ld)
  
  # PROJ_INF IS IN 'solvers' FOLDER!
  par$rr <- pi_out$rr
  y <- pi_out$y
  
  #return(c(mean(par$rr),mean(y)))
  
  Rpb <- Ax0 - b + xi0
  normRp <- norm(Rpb, "2")
  ytmp <- yinput - y
  Atxi <- Atxi0
  xi <- xi0
  Ly <- t(b) %*% xi - 0.5 * norm(xi,"2")^2 - 0.5 * sig * norm(ytmp, "2")^2
  
  runhist <- list()
  runhist$psqmr[1] <- 0
  runhist$findstep[1] <- 0
  cnt_Amap <- 0
  cnt_ATmap <- 0
  cnt_pAATmap <- 0
  cnt_fAATmap <- 0
  
  for (itersub in seq(1, maxitersub)) {    
    yold <- y
    xiold <- xi
    Atxiold <- Atxi
    Rdz <- Atxi + y;
    normRd <- norm(Rdz,"2")
    msigAytmp <- -sig * A %*% ytmp
    GradLxi <- -(xi - b + msigAytmp)
    cnt_Amap <- cnt_Amap + 1
    normGradLxi <- norm(GradLxi, "2") * sqrt(bscale * cscale) / normborg
    priminf_sub <- normGradLxi
    
    if(Ascaleyes == 1) {
      dualinf_sub <- norm(Rdz / options$dscale, "2") * cscale / 
        (1 + norm(y / options$dscale, "2") * cscale)
    } else {
      dualinf_sub <- normRd * cscale/ (1 + norm(y, "2") * cscale)
    }
    
    
    if(max(priminf_sub,dualinf_sub) < tol) {
      tolsubconst <- 0.9
    } else {
      tolsubconst <- 0.05
    }
    
    tolsub <- max(min(1, par$tolconst * dualinf_sub), tolsubconst * tol)
    runhist$priminf[itersub] <- priminf_sub
    runhist$dualinf[itersub] <- dualinf_sub
    runhist$Ly[itersub] <- Ly
    
    if (printsub){
      # cat("iteration for subproblem:", itersub, "\n")
      # cat("Ly:", Ly, "\n")
      # cat("priminf_sub:", priminf_sub, "\n")
      # cat("dualinf_sub:", dualinf_sub, "\n")
      # cat("par.tolconst:", par$tolconst, "\n")
    }
    #if (printsub)
    #  fprintf('\n      %2.0d  %- 11.10e %3.2e %3.2e %3.2e',...
    #          itersub,Ly,priminf_sub,dualinf_sub,par.tolconst);
    #end
    
    if (max(normGradLxi) < max(tolsub) && itersub > 1) {
      if (printsub) {
        cat("good termination in subproblem:", "\n")
        cat("dualinf in subproblem:", dualinf_sub, "\n")
        cat("tolsub:", tolsub, "\n")
        cat("normGradLxi:", normGradLxi, "\n")
      }
      #msg = 'good termination in subproblem:';
      #if printsub
      #fprintf('\n       %s  ',msg);
      #fprintf(' dualinfes = %3.2e, gradLyxi = %3.2e, tolsub = %3.2e',...
      #        dualinf_sub,normGradLxi,tolsub);
      #end
      breakyes <- -1
      break
    }
    #%% Compute Newton direction
    #%% precond = 0, 
    
    par$epsilon <- min(1e-3, 0.1 * normGradLxi) #%% good to add
    par$precond <- precond
    if(precond == 1) {
      par$invdiagM <- 1 / (1 + sig)
    }
    if( (dualinf_sub > 1e-3) | (itersub <= 5) ) {
      maxitpsqmr <- max(maxitpsqmr, 200)
    } else if (dualinf_sub > 1e-4) {	 
      maxitpsqmr <- max(maxitpsqmr, 300) 
    } else if (dualinf_sub > 1e-5) {	 
      maxitpsqmr <- max(maxitpsqmr, 400) 
    } else if (dualinf_sub > 5e-6) {
      maxitpsqmr <- max(maxitpsqmr, 500) 
    }
    
    if (itersub > 1) {
      prim_ratio <- priminf_sub / runhist$priminf[itersub - 1]
      dual_ratio <- dualinf_sub / runhist$dualinf[itersub - 1]
      if (is.nan(dual_ratio)) {
        prim_ratio <- 0
        dual_ratio <- 0
      }
    } else {
      prim_ratio <- 0
      dual_ratio <- 0
    }
    
    rhs <- GradLxi
    
    #if (Ascaleyes == 1 && false
    #tolpsqmr = min(5e-3, 0.1*norm(rhs));
    #else
    #end
    
    tolpsqmr <- min(5e-3, 0.1*norm(rhs,"2"))
    
    const2 <- 1
    if (itersub > 1 &
        (prim_ratio > 0.5 | priminf_sub > 0.1 * runhist$priminf[1])){
      const2 <- 0.5 * const2
    }
    
    if (dual_ratio > 1.1){
      const2 <- 0.5 * const2
    }
    
    tolpsqmr <- const2 * tolpsqmr
    par$tol <- tolpsqmr
    par$maxit <- maxitpsqmr
    
    
    lsout <- linsyssolve(A,rhs,par) #replacing Classic_Lasso_linsys_solver
    dxi <- lsout$xi
    resnrm <- lsout$resnrm
    solve_ok <- lsout$solve_ok
    
    # Atdxi = t(A) %*% dxi
    Atdxi <- t(t(dxi) %*% A)
    # Atdxi <- eigenMapMatMult(t(A), dxi, 4)
    cnt_ATmap <- cnt_ATmap + 1
    iterpsqmr <- length(resnrm) - 1
    
    if (iterpsqmr == 0) {
      cnt_pAATmap <- cnt_pAATmap + 1
    } else {
      if (existA == TRUE) {
        cnt_pAATmap <- cnt_pAATmap + iterpsqmr
      } else {
        cnt_fAATmap <- cnt_fAATmap + iterpsqmr
      }
    }
    
    if (printsub) {
      cat("tolerance:", par$tol, "\n")
      cat("resnrm:", tail(resnrm, n=1), "\n")
      cat("PSQMR iteration no.:", iterpsqmr, "\n")
    }
    #if (printsub)
    #  fprintf('| %3.1e %3.1e %3.0d %-3d',par.tol,resnrm(end),iterpsqmr);
    #fprintf(' %2.1f %2.0d',const2,sum(1-par.rr));
    #end
    
    par$iter <- itersub;
    if ((itersub <= 3) & (dualinf_sub > 1e-4) | (par$iter <3)) {
      stepop <- 1
    } else {
      stepop <- 2
    }
    steptol <- 1e-5
    step_op <- list()
    step_op$stepop <- stepop
    
    
    fsret <- findstep(par, b, ld, Ly, xi, Atxi, y, ytmp, dxi, Atdxi,
                      steptol, step_op)
    par <- fsret$par
    Ly <- fsret$Ly
    xi <- fsret$xi
    Atxi <- fsret$Atxi
    y <- fsret$y
    ytmp <- fsret$ytmp
    alp <- fsret$alp
    iterstep <- fsret$iter
    
    runhist$solve_ok[itersub] <- solve_ok
    runhist$psqmr[itersub] <- iterpsqmr 
    runhist$findstep[itersub] <- iterstep 
    if (alp < tiny) breakyes <- 1
    Ly_ratio <- 1
    if (itersub > 1) {
      Ly_ratio <- (Ly - runhist$Ly[itersub - 1]) / (abs(Ly) + eps)
    }
    
    if (printsub){
      cat("alp:", alp, "\n")
      cat("iterstep:", iterstep, "\n")
      if (Ly_ratio < 0){
        cat("-----------------------------", "\n")
      }
    }
    #if (printsub)
    #  fprintf(' %3.2e %2.0f',alp,iterstep);
    #if (Ly_ratio < 0); fprintf('-'); end
    #end
    
    ### CHECK FOR STAGNATION LINE 161 GitHub code
    #%% check for stagnation
    if (itersub > 4) {
      idx <- seq(max(1, itersub - 3) : itersub)
      tmp <- runhist$priminf[idx]
      ratio <- min(tmp) / max(tmp)
      if ((all(runhist$solve_ok[idx] <= -1)) & (ratio > 0.9)  
          & (min(runhist$psqmr[idx]) == max(runhist$psqmr[idx])) 
          & (max(tmp) < 5 * tol)) {
        breakyes <- 1
      }
      const3 <- 0.7
      priminf_1half <- min(runhist$priminf[1:ceil(itersub * const3)])
      priminf_2half <- min(runhist$priminf[(ceil(itersub * const3)+1):itersub])
      priminf_best <- min(runhist$priminf[1:itersub-1])
      priminf_ratio <- runhist$priminf[itersub] / runhist$priminf[itersub-1]
      dualinf_ratio <- runhist$dualinf[itersub] / runhist$dualinf[itersub-1]
      stagnate_idx <- which(runhist$solve_ok[1:itersub] <= -1)
      stagnate_count <- length(stagnate_idx)
      idx2 <- seq(max(1, itersub - 7): itersub)
      
      
      if ((itersub >= 10) & all(runhist$solve_ok[idx2] == -1)  
          & (priminf_best < 1e-2) & (dualinf_sub < 1e-3) ) {                   
        tmp <- runhist$priminf[idx2] 
        ratio <- min(tmp) / max(tmp)
        if (ratio > 0.5) {
          #if (printsub) fprintf('##')
          breakyes <- 2 
        }
      }
      
      if ((itersub >= 15) & (priminf_1half < min(2e-3, priminf_2half)) & 
          (dualinf_sub < 0.8 * runhist$dualinf[1]) & (dualinf_sub < 1e-3) & 
          (stagnate_count >= 3)){
        #if (printsub); fprintf('###'); end
        breakyes <- 3
      }
      
      if ((itersub >= 15) & (priminf_ratio < 0.1) 
          & (priminf_sub < 0.8 * priminf_1half) 
          & (dualinf_sub < min(1e-3, 2 * priminf_sub)) 
          & ((priminf_sub < 2e-3) | 
             (dualinf_sub < 1e-5 & priminf_sub < 5e-3)) 
          & (stagnate_count >= 3) ) {
        #if (printsub); fprintf(' $$'); end
        breakyes <- 4
      }
      
      
      if ((itersub >=10) & (dualinf_sub > 5 * min(runhist$dualinf)) 
          & (priminf_sub > 2 * min(runhist$priminf))) {
        #if (printsub); fprintf('$$$'); end
        breakyes <- 5
      }
      
      if (itersub >= 20) {
        dualinf_ratioall <- runhist$dualinf[2:itersub] / 
          runhist$dualinf[1:(itersub-1)]
        idx <- which(dualinf_ratioall > 1) 
        if (length(idx) >= 3) {
          dualinf_increment <- mean(dualinf_ratioall[idx])
          if (dualinf_increment > 1.25) {
            breakyes <- 6             
          }                    
        }              
      }
      
      if (breakyes > 0) {
        Rdz <- Atxi + y
        # msigAytmp = -sig*A %*% ytmp
        msigAytmp <- eigenMapMatMult(-sig*A, ytmp, 4)
        cnt_Amap <- cnt_Amap + 1
        if (printsub){
          dualfeasorg <- norm(Rdz/options$dscale, "2")*cscale/(1+norm(y/options$dscale, "2")*cscale)
          cat("dualfeasorg:", dualfeasorg, "\n")
        }
        #%normRd = norm(Rdz);
        #if printsub
        #if Ascaleyes
        #fprintf('\n new dualfeasorg = %3.2e', norm(Rdz./options.dscale)*cscale/(1+norm(y./options.dscale)*cscale));
        #else
        #  fprintf('\n new dualfeasorg = %3.2e', norm(Rdz)*cscale/(1+norm(y)*cscale));
        #end
        #end
        break
      }
    }
  }
  info <- list()
  info$maxCG <- max(runhist$psqmr)
  info$avgCG <- sum(runhist$psqmr) / itersub
  info$breakyes <- breakyes
  info$itersub <- itersub
  info$tolconst <- par$tolconst
  info$RpGradratio <- normRp * sqrt(bscale * cscale)/(normGradLxi * normborg)
  info$rankX <- par$rr
  info$ytmp <- ytmp
  info$cnt_Amap <- cnt_Amap
  info$cnt_ATmap <- cnt_ATmap
  info$Ax <- msigAytmp
  info$cnt_pAATmap <- cnt_pAATmap
  info$cnt_fAATmap <- cnt_fAATmap
  
  
  output <- list(y = y,
                 Atxi = Atxi,
                 xi = xi,
                 par = par,
                 runhist = runhist,
                 info = info)
  return(output)
}