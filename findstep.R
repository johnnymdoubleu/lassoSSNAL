#par,Ly,xi,Atxi,y,ytmp,alp,iterstep


#output <- list(par = par,
#               Ly = Ly,
#               xi = xi,
#               Atxi = Atxi,
#               y = y,
#               ytmp = ytmp,
#               alp = alp,
#               iter = iter)
#return(output)

findstep <- function(par,b,ld,Ly0,xi0,Atxi0,y0,ytmp0,
                     dxi,Atdxi,tol,options) {
  if("stepop" %in% names(options)) {
    stepop = options$stepop
  }
  printlevel=0
  maxit = ceil(log(1/(tol+eps))/log(2))
  c1 = 1e-4
  c2 = 0.9
  sig = par$sigma
  #%%
  
  tmp1 = t(dxi) %*% (b-xi0)
  tmp2 = norm(dxi,"2")^2 
  g0  = tmp1 + sig * t(Atdxi) %*% ytmp0
  #Ly = []; 
  
  if (g0 <= 0) {
    alp = 0
    iter = 0
    #if (printlevel) 
    #  fprintf('\n Need an ascent direction, %2.1e  ',g0); 
    #end
    xi = xi0
    Atxi = Atxi0
    y = y0
    ytmp = ytmp0
    Ly = Ly0
    
    ## FINALIZE THIS!!
    output <- list(par = par,
                   Ly = Ly,
                   xi = xi,
                   Atxi = Atxi,
                   y = y,
                   ytmp = ytmp,
                   alp = alp,
                   iter = iter)
    return(output)
  }
  
  alp = 1 
  alpconst = 0.5
  for (iter in seq(1,maxit)) {
    if (iter==1) {          
      alp = 1
      LB = 0
      UB = 1
    } else {
      alp = alpconst*(LB+UB)
    }
    xi = xi0 + alp*dxi
    yinput = ytmp0 + y0 - alp*Atdxi
    
    ## DEAL WITH THIS!
    
    #[y,rr] = proj_inf(yinput,ld);
    #par.rr = rr;
    
    pi_out <- proj_inf(yinput,ld)
    
    # PROJ_INF IS IN 'solvers' FOLDER!
    par$rr = pi_out$rr
    y = pi_out$y
    
    ytmp = yinput - y;
    galp = tmp1 - alp*tmp2 + sig*t(Atdxi) %*% ytmp;
    if (iter==1) {
      gLB = g0
      gUB = galp
      if (sign(gLB)*sign(gUB) > 0) {
        #if (printlevel); fprintf('|'); end
        Atxi = Atxi0+alp*Atdxi;
        Ly = t(b) %*% xi - 0.5*norm(xi,"2")^2 - 0.5*sig*norm(ytmp,"2")^2;             
        
        ## DEAL WITH THIS RETURN VALUE LATER
        output <- list(par = par,
                       Ly = Ly,
                       xi = xi,
                       Atxi = Atxi,
                       y = y,
                       ytmp = ytmp,
                       alp = alp,
                       iter = iter)
        return(output)
      }
    }
    
    if (abs(galp) < c2*abs(g0)) { 
      Ly = t(b) %*% xi - 0.5*norm(xi,"2")^2 - 0.5*sig*norm(ytmp,"2")^2
      if (Ly-Ly0-c1*alp*g0 > -1e-8/max(1,abs(Ly0))
          && ((stepop==1) || (stepop==2 && 
                              abs(galp)<tol))) {
        #if (printlevel); fprintf(':'); end
        Atxi = Atxi0+alp*Atdxi;
        
        ### HERE!!!
        output <- list(par = par,
                       Ly = Ly,
                       xi = xi,
                       Atxi = Atxi,
                       y = y,
                       ytmp = ytmp,
                       alp = alp,
                       iter = iter)
        return(output)         
      }
    }
    
    if (sign(galp)*sign(gUB) < 0) {
      LB = alp
      gLB = galp
    } else if (sign(galp)*sign(gLB) < 0) { 
      UB = alp
      gUB = galp
    }
  }
  if (iter == maxit) {
    Atxi = Atxi0+alp*Atdxi
  }
  #if (printlevel); fprintf('m'); end
  if (!exists("Ly")) {
    Ly = t(b) %*% xi - 0.5*norm(xi,"2")^2 - 0.5*sig*norm(ytmp,"2")^2            
  }
  
  output <- list(par = par,
                 Ly = Ly,
                 xi = xi,
                 Atxi = Atxi,
                 y = y,
                 ytmp = ytmp,
                 alp = alp,
                 iter = iter)
  return(output)
}