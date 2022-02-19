Classic_Lasso_SSNAL <- function(Ainput, b, n, lambda, options, y=NULL, xi=NULL, x=NULL){
  maxiter <- 5000
  stoptol <- 1e-6
  printyes <- 1
  scale <- 0
  dscale <- matrix(1,n,1)
  rescale <- 1
  Lip <- 1
  Ascale <- 0
  Ascaleyes <- 0
  orgojbconst <- 0
  
  
  #replace the assigned values if the values in options were defined
  if ("maxiter" %in% names(options)){maxiter <- options$maxiter}
  if ("stoptol" %in% names(options)){stoptol <- options$stoptol}
  if ("printyes" %in% names(options)){printyes <- options$printyes}
  if ("rescale" %in% names(options)){rescale <- options$rescale}
  if ("Lip" %in% names(options)){Lip <- options$Lip}
  if ("Ascale" %in% names(options)){Ascale <- options$Ascale}
  if ("orgojbconst" %in% names(options)){orgojbconst <- options$orgojbconst}
  
  #start time / clock time     tstart
  tstart <- as.numeric(strsplit(format(Sys.time(), "%Y %m %d %H %M %S")," ")[[1]])/rep(1000,6)
  m <- length(b)
  #if (is.list(Ainput)){
  #  A <- Ainput$A
  #  Amap0 <- Ainput.Amap0
  #  ATmap0 <- Ainput.ATmap0
  #}
  #else {
  #  A <- Ainput
  #  Amap <- A%*%x
  #  ATmap <- t(A) %*% x
  #}
  
  A <- Ainput
  existA <- exists("A")
  if (Ascale!=0 & existA) {
    #start 2nd timer
    if (Ascale == 1) {
      dscale <- 1/pmax(1, t(sqrt(colSums(A*A))))
    }
    else if (Ascale == 2){
      dscale <- 1 / sqrt(colSums(A*A))
    }
    
    #A <- A %*% diag(dscale,n,n)
    A <- as.matrix(A%*%Diagonal(x=dscale))
    Ascaleyes <- 1
    #print out the time taken
  }

  if(is.null(x) | is.null(xi) | is.null(y)){
    x <- matrix(0,n,1)
    xi <- matrix(0,m,1)
    y <- x
  }
  
  
  parmain <- list(dscale = dscale,
                  Ascaleyes = Ascaleyes,
                  m = m,
                  n = n,
                  scale = scale,
                  existA = existA,
                  orgojbconst = orgojbconst,
                  A = A,
                  tstart = tstart,
                  Lip = Lip,
                  maxiter = maxiter,
                  printyes = printyes,
                  rescale = rescale,
                  stoptol = stoptol
  )
  if("Sigma" %in% names(options)) parmain$Sigma = options$Sigma
  
  clsmo <- Classic_Lasso_SSNAL_main(A,Ainput,b,lambda,parmain,y,xi,x)
  
  obj_main <- clsmo$obj
  y <- clsmo$y
  xi <- clsmo$xi
  x <- clsmo$x
  info_main <- clsmo$info
  runhist <- clsmo$runhist
  
  iter <- info_main$iter
  bscale <- info_main$bscale
  cscale <- info_main$cscale
  Ax <- info_main$Ax
  #ttime <- info_main$ttime
  msg <- info_main$msg
  
  if (iter == maxiter){
    msg <- "Maximum iteration reached"
  }
  
  xi <- xi * sqrt(bscale * cscale)
  # At <- t(Ainput)
  # Atxi <- t(Ainput) %*% xi
  Atxi <- t(t(xi) %*% Ainput)
  # Atxi <- eigenMapMatMult(At, xi, 4)
  y <- y * cscale
  x <- x * bscale
  if (Ascaleyes) {
    x <- dscale * x
    y <- dscale / y
  }
  Rd <- Atxi + y
  dualfeasorg <- norm(Rd, "2")/ (1 + norm(y,"2"))
  Ax <- Ax * sqrt(bscale * cscale)
  Rp <- Ax - b + xi
  primfeasorg <- norm(Rp, "2")/ (1+norm(b,"2"))
  primobj <- 0.5 * norm(Ax-b,"2")^2 + lambda*norm(x) + orgojbconst
  dualobj <- -0.5 * norm(xi,"2")^2 + t(b) %*% xi + orgojbconst
  
  relgap <- (primobj-dualobj)/(1+abs(primobj)+abs(dualobj))
  obj <- c(primobj, dualobj)
  grad <- t(t(Ax-b)%*%Ainput)
  # grad <- eigenMapMatMult(At, (Ax-b), 4)
  # grad <- t(Ainput) %*% (Ax-b)
  etaorg <- norm(grad + proj_inf(x - grad, lambda)$y, "2")
  eta <- etaorg / (1 + norm(grad,"2") + norm(x, "2"))
  
  #runhist <- list(m = m,
  #                n=n,
  #                iter = iter,
  #                totaltime = ttime,
  #                primobjorg = primobj,
  #                dualobjorg = dualobj,
  #                maxfeas = max(c(dualfeasorg, primfeasorg)),
  #                eta = eta,
  #                etaorg = etaorg)
  #runhist <- list()
  
  info <- list(m = m,
               n = n,
               minx = min(min(x)),
               max = max(max(x)),
               relgap = relgap,
               iter = iter,
   #            totaltime = ttime,
               eta = eta,
               etaorg = etaorg,
               obj = obj,
               maxfeas = max(c(dualfeasorg, primfeasorg)),
    #           cnt_Amap = sum(runhist$cnt_Amap),
    #           cnt_ATap = sum(runhist$cnt_ATmap),
    #           cnt_pAATmap = sum(runhist$cnt_pAATmap),
    #           cnt_fAATmap = sum(runhist$cnt_fAATmap),
    #           nnz = findnnz(x, 0.999),
               x = x
  )
  #if (printyes){
  #  if (is.null(msg)){
  #    printf('\n %d', msg)
  #  }
  #  printf('\n  number iter = %2.0d',iter)
  #  printf('\n  time = %3.2d',ttime)
  #  printf('\n  time per iter = %5.4d', ttime/iter)
  #  printf('\n  primobj = %9.8d, dualobj = %9.8d, relgap = %3.2d',primobj,dualobj, relgap)
  #  printf('\n  primfeasorg = %3.2d, dualfeasorg = %3.2d',primfeasorg, dualfeasorg)
  #  printf('\n  eta = %3.2, etaorg = %3.2', eta, etaorg)
  #  printf('\n  min(X) = %3.2d, max(X) = %3.2d', info$minx,info$maxx)
  #  printf('\n  Amap cnt = %3d, ATmap cnt = %3d, partial AATmap cnt = %3d, full AATmap cnt = %3d', info$cnt_Amap, info$cnt_ATmap, info$cnt_pAATmap, info$cnt_fAATmap)
  #  printf('\n  number of nonzeros in x (0.999) = %3.0d',findnnz(x,0.999))
  #}
  
  
  output <- list(obj = obj,
                 y = y,
                 xi = xi,
                 x = x,
                 info = info,
                 runhist = runhist)
}