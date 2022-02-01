
function [obj,y,xi,x,info,runhist] = Classic_Lasso_SSNAL_main(Amap0,ATmap0,b,lambda,parmain,y,xi,x) 

   dscale    = parmain.dscale;
   Ascaleyes = parmain.Ascaleyes;
   m = parmain.m; n = parmain.n;
   
   if isfield(parmain,'A'); A = parmain.A; end   
   tstart = parmain.tstart;
   existA = parmain.existA;
   
   Lip = parmain.Lip;
   scale = parmain.scale;
   
   maxiter = parmain.maxiter;
   printyes = parmain.printyes;
   rescale = parmain.rescale;
   stoptol = parmain.stoptol;
   orgojbconst = parmain.orgojbconst;

   sigmaLip = 1/Lip;
   if norm(dscale - 1) < eps
      ld = lambda;
   else
      ld = lambda*dscale;
   end
   lambdaorg = lambda;
   borg = b;
   normborg = 1 + norm(borg);
   Atxi = ATmap0(xi); Ax = Amap0(x);
   obj(1) = 0.5*norm(Ax - borg)^2 + lambdaorg*norm(x,1) + orgojbconst;
   obj(2) = -(0.5*norm(xi)^2 + borg'*xi) + orgojbconst;

   %DBG
   %fprintf('\n  obj1 = %3.2f',obj(1)); %DBGGG
   %fprintf('\n  obj2 = %3.2f',obj(2));
   %fprintf('\n  dscale = %3.2f',parmain.dscale); %DBGGG
   %input("Stop") %DBG

   bscale = 1; cscale = 1; 
   if (scale == 1)
      b = b/sqrt(bscale*cscale);
      xi = xi/sqrt(bscale*cscale);
      Amap = @(x) Amap0(x*sqrt(bscale/cscale));
      ATmap = @(x) ATmap0(x*sqrt(bscale/cscale));
      if existA; A = A*sqrt(bscale/cscale); end
      lambda = lambda/cscale; ld = ld/cscale;
      x = x/bscale; y = y/cscale; 
      Ax = Ax/sqrt(bscale*cscale);
      Atxi = Atxi/cscale;
   else
      Amap = Amap0; ATmap = ATmap0;
   end

   %fprintf('\n  ascaleyes = %3.2f',Ascaleyes); %DBGGG
   %input("Stop") %DBG

   if Ascaleyes
      Amap = @(x) Amap(dscale.*x);
      ATmap = @(x) dscale.*ATmap(x);
      Atxi = dscale.*Atxi; y = dscale.*y; 
   end
   normb = 1+norm(b);
   Ainput_nal.Amap = Amap;
   Ainput_nal.ATmap = ATmap;
   if existA; Ainput_nal.A = A; end
   sigma = max(1/sqrt(Lip),min([1,sigmaLip,lambdaorg]));
   

   if Ascaleyes; 
      sigma = 3; 
   end
   if isfield(parmain,'sigma'); sigma = parmain.sigma; end

   %fprintf('\n  sigma = %3.7f',sigma); %DBGGG
   %input("Stop") %DBG

   Rp1 = Ax - b;
   Rp = Rp1 + xi;
   Rd = Atxi + y;
   primfeas = norm(Rp)/normb;
   dualfeas = norm(Rd)/(1+norm(y));
   primfeasorg = sqrt(bscale*cscale)*norm(Rp)/normborg;
   dualfeasorg = norm(Rd./dscale)*cscale/(1+norm(y./dscale)*cscale);
   maxfeas = max(primfeas, dualfeas);
   maxfeasorg = max(primfeasorg, dualfeasorg);
   relgap = (obj(1) - obj(2))/(1+obj(1)+obj(2));
   runhist.dualfeasorg(1) = dualfeasorg;
   runhist.primfeasorg(1) = primfeasorg;
  
   %{
   fprintf('\n  primfeas = %3.7f',primfeas); %DBGGG
   fprintf('\n  dualfeas = %3.7f',dualfeas); %DBGGG
   fprintf('\n  primfeasorg = %3.7f',primfeasorg); %DBGGG
   fprintf('\n  dualfeasorg = %3.7f',dualfeasorg); %DBGGG
   fprintf('\n  maxfeas = %3.7f',maxfeas); %DBGGG
   fprintf('\n  maxfeasorg = %3.7f',maxfeasorg); %DBGGG
   fprintf('\n  relgap = %3.7f',relgap); %DBGGG
   input("Stop") %DBG
   %}

   if printyes
        fprintf('\n*******************************************************');
        fprintf('******************************');
        fprintf('\n \t\t   Classic Lasso: SSNAL      ');
        fprintf('\n******************************************************');
        fprintf('*******************************');
        fprintf('\n n = %3.0f, m = %3.0f',n, m);
        fprintf('\n bscale = %3.2e, cscale = %3.2e', bscale, cscale);
        fprintf('\n ---------------------------------------------------');
        fprintf('\n  iter|  [pinfeas  dinfeas]  [pinforg  dinforg]    relgaporg|    pobj          dobj    |');
        fprintf(' time | sigma |rankS |');
        fprintf('\n*****************************************************');
        fprintf('**************************************************************');
        fprintf('\n #%3.1d|  %3.2e %3.2e %3.2e %3.2e %- 3.2e %- 8.7e %- 8.7e  %5.1f',...
           0,primfeas,dualfeas,primfeasorg,dualfeasorg,relgap,obj(1),obj(2),etime(clock,tstart)); 
        fprintf('  %3.2e ',sigma);
   end
   %% ssncg
   SSNCG = 1;
   if SSNCG
      parNCG.sigma = sigma;
      parNCG.tolconst = 0.5;
      parNCG.n = n;
   end
   maxitersub = 10;
   breakyes = 0;
   prim_win = 0;
   dual_win = 0;
   ssncgop.existA = existA;
   ssncgop.tol = stoptol;
   ssncgop.precond = 0;
   ssncgop.bscale = bscale;
   ssncgop.cscale = cscale;
   ssncgop.printsub = printyes;
   if Ascaleyes
      ssncgop.Ascaleyes = 1;
      ssncgop.dscale = dscale;
   else
      ssncgop.Ascaleyes = 0;
   end
   sigmamax = 1e7; %1e5;
   sigmamin = 1e-4; 
   if Ascaleyes; sigmamax = sigmamax*mean(dscale); end
   
   %fprintf('\n  normyxi = %3.7f',(max(normx/normyxi,normyxi/normx) > 1.7)); %DBGGG
   %input("Stop") %DBG

   for iter = 1:maxiter
      if ((rescale == 1) && (maxfeas < 5e2) && (rem(iter,3) == 1) && (iter > 1) ) ...
          || (~exist('A','var') && ((rescale >= 2) && maxfeas < 1e-1 && (abs(relgap) < 0.05) ...
             && (iter >= 5) && (max(normx/normyxi,normyxi/normx) > 1.7) && rem(iter,5)==1))
         normy = norm(y); normAtxi = norm(Atxi);
         normx = norm(x);         
         normyxi = max(normx,normAtxi); 
         [sigma,bscale2,cscale2,sbc,sboc,bscale,cscale] = ...
             mexscale(sigma,normx,normyxi,bscale,cscale);
         Amap = @(x) Amap(x*sboc);
         ATmap = @(x) ATmap(x*sboc);
         Ainput_nal.Amap = Amap;
         Ainput_nal.ATmap = ATmap;
         if exist('A','var'); Ainput_nal.A = A*sboc; end
         b = b/sbc; 
         ld = ld/cscale2; 
         x = x/bscale2;
         xi = xi/sbc; 
         Atxi = Atxi/cscale2; Ax = Ax/sbc;         
         ssncgop.bscale = bscale;
         ssncgop.cscale = cscale;
         normb = 1+norm(b);
         if printyes
            fprintf('\n    ');
            fprintf('[rescale=%1.0f: %2.0f| %3.2e %3.2e %3.2e | %3.2e %3.2e| %3.2e]',...
            rescale,iter,normx,normAtxi, normy,bscale,cscale,sigma);
         end         
         rescale = rescale+1; 
         prim_win = 0; dual_win = 0; 
      end 
      xold = x;      
      parNCG.sigma = sigma;
      if dualfeas < 1e-5
         maxitersub = max(maxitersub,30);
      elseif dualfeas < 1e-3
         maxitersub = max(maxitersub,30);
      elseif dualfeas < 1e-1
         maxitersub = max(maxitersub,20);
      end
      ssncgop.maxitersub = maxitersub; 
%{
      fprintf('\n  n = %3.7f',n); %DBGGG
      fprintf('\n  b = %3.7f',mean(b)); %DBGGG
      fprintf('\n  x = %3.7f',mean(x)); %DBGGG
      fprintf('\n  Ax = %3.7f',mean(Ax)); %DBGGG
      fprintf('\n  Atxi = %3.7f',mean(Atxi)); %DBGGG
      fprintf('\n  xi = %3.7f',mean(xi)); %DBGGG
      fprintf('\n  ld = %3.7f',mean(ld)); %DBGGG
      fprintf('\n ssncgop_precond = %3.7f',ssncgop.precond);
      fprintf('\n ssncgop_maxitersub = %3.7f',ssncgop.maxitersub);
      input("Stop") %DBG
%}
      [y,Atxi,xi,parNCG,runhist_NCG,info_NCG] = ...
          Classic_Lasso_SSNCG(n,b,Ainput_nal,x,Ax,Atxi,xi,ld,parNCG,ssncgop);
      
      %{
      fprintf('\n sy = %3.7f',sum(y));
      fprintf('\n satxi = %3.7f',sum(Atxi));
      fprintf('\n sxi = %3.7f',sum(xi));
      fprintf('\n breakyes = %3.7f',info_NCG.breakyes);
      input("Stop") %DBG
      %}

      if info_NCG.breakyes < 0
         parNCG.tolconst = max(parNCG.tolconst/1.06,1e-3);
      end

      fprintf('\n parncg tolconst = %3.7f',parNCG.tolconst);
      %input("Stop") %DBG

      Rd = Atxi + y;
      x = -sigma*info_NCG.ytmp;
      Ax = info_NCG.Ax;
      normRd = norm(Rd); normy = norm(y);

      fprintf('\n normrd = %3.7f',normRd);
      fprintf('\n normy = %3.7f',normy);
      %input("Stop") %DBG

      dualfeas = normRd/(1+normy);

      fprintf('\n dualfeas = %3.7f',dualfeas);
      %input("Stop") %DBG

      if Ascaleyes
         dualfeasorg = norm(Rd./dscale)*cscale/(1+norm(y./dscale)*cscale);
      else
         dualfeasorg = normRd*cscale/(1+normy*cscale);
      end

      fprintf('\n dualfeasorg = %3.7f',dualfeasorg);
      %input("Stop") %DBG

      Rp1 = Ax - b;
      Rp = Rp1 + xi;
      normRp = norm(Rp);
      primfeas = normRp/normb;
      primfeasorg = sqrt(bscale*cscale)*normRp/normborg;
      maxfeas = max(primfeas,dualfeas);
      maxfeasorg = max(primfeasorg, dualfeasorg);

      fprintf('\n normRp = %3.7f',normRp);
      fprintf('\n primfeas = %3.7f',primfeas);
      fprintf('\n primfeasorg = %3.7f',primfeasorg);
      fprintf('\n maxfeas = %3.7f',maxfeas);
      fprintf('\n maxfeasorg = %3.7f',maxfeasorg);
      %input("Stop") %DBG


      runhist.dualfeas(iter+1) = dualfeas;
      runhist.primfeas(iter+1) = primfeas;
      runhist.maxfeas(iter+1)  = maxfeas;
      runhist.primfeasorg(iter+1) = primfeasorg;
      runhist.dualfeasorg(iter+1) = dualfeasorg;
      runhist.maxfeasorg(iter)  = maxfeasorg;
      runhist.sigma(iter) = sigma;
      runhist.rankS(iter) = sum(1 - parNCG.rr);
      runhist.cnt_Amap(iter) = info_NCG.cnt_Amap;
      runhist.cnt_ATmap(iter) = info_NCG.cnt_ATmap;
      runhist.cnt_pAATmap(iter) = info_NCG.cnt_pAATmap;
      runhist.cnt_fAATmap(iter) = info_NCG.cnt_fAATmap;
      %%---------------------------------------------------------
%% check for termination
%%---------------------------------------------------------    
      if (max([primfeasorg,dualfeasorg]) < 500*max(1e-6, stoptol)) 
          %fprintf('\n Stopping here!');
          %input("Stop") %DBG

          grad = ATmap0(Rp1*sqrt(bscale*cscale));
          if Ascaleyes
             etaorg = norm(grad + proj_inf(dscale.*x*bscale - grad,lambdaorg));
             eta = etaorg/(1+norm(grad)+norm(dscale.*x*bscale));
          else
             etaorg = norm(grad + proj_inf(x*bscale - grad,lambdaorg));
             eta = etaorg/(1+norm(grad)+norm(x*bscale));
          end
          if eta < stoptol
             breakyes = 1;
             msg = 'converged';
          end
      end
      %input("Stop") %DBG
%%--------------------------------------------------------    
%% print results
%%--------------------------------------------------------
      
      if printyes
         objscale = bscale*cscale;
         primobj = objscale*(0.5*norm(Rp1)^2 + norm(ld.*x,1)) + orgojbconst;
         dualobj = objscale*(-0.5*norm(xi)^2 + b'*xi) + orgojbconst;
         relgap = (primobj-dualobj)/( 1+abs(primobj)+abs(dualobj));
         ttime = etime(clock,tstart);
         if (printyes)
            fprintf('\n %5.0d| [%3.2e %3.2e] [%3.2e %3.2e]  %- 3.2e| %- 5.4e %- 5.4e |',...
               iter,primfeas,dualfeas,primfeasorg, dualfeasorg,relgap,primobj,dualobj); 
            fprintf(' %5.1f| %3.2e|sigamorg = %3.2e  |',ttime, sigma, sigma*bscale/cscale); 
            if iter >= 1
               fprintf('%3.0d|',sum(1- parNCG.rr));
            end
            if exist('eta'); fprintf('\n        [eta = %3.2e, etaorg = %3.2e]',eta, etaorg);end
         end
	     if (rem(iter,3*1)==1) 
            normx = norm(x); normAtxi = norm(Atxi); normy = norm(y);
            if (printyes)
               fprintf('\n        [normx,Atxi,y =%3.2e %3.2e %3.2e]',...
               normx,normAtxi,normy);
            end
         end
         runhist.primobj(iter)   = primobj;
         runhist.dualobj(iter)   = dualobj;
         runhist.time(iter)      = ttime; 
         runhist.relgap(iter)    = relgap;
      end

        fprintf('\n primobj=%3.7f',primobj);
        fprintf('\n dualobj=%3.7f',dualobj);
        fprintf('\n relgap=%3.7f',relgap);
        %input("Stop") %DBG

      if (breakyes > 0) 
         if printyes
            fprintf('\n  breakyes = %3.1f, %s',breakyes,msg); 
         end
         break; 
      end
      if (primfeasorg < dualfeasorg)
         prim_win = prim_win+1; 
      else
         dual_win = dual_win+1; 
      end   

      fprintf('\n prim_win = %3.7f',prim_win);
      fprintf('\n dual_win = %3.7f',dual_win);
      %input("Stop") %DBG

      [sigma,prim_win,dual_win] = ...
          mexsigma_update_Classic_Lasso_SSNAL(sigma,sigmamax,sigmamin,prim_win,dual_win,iter,info_NCG.breakyes);
   
      fprintf('\n sigma=%3.7f',sigma);
      fprintf('\n prim_win = %3.7f',prim_win);
      fprintf('\n dual_win = %3.7f',dual_win);
      %input("Stop") %DBG

   end
   
   if ~printyes 
      ttime = etime(clock,tstart);
   end
   
   if (iter == maxiter) && (breakyes == 0) %% added by X. D. Li @2018/10/23/15:14
      msg = 'maximum number of iterations reached';
   end
   
   info.iter = iter;
   info.bscale = bscale;
   info.cscale = cscale;
   info.Ax = Ax;
   info.ttime = ttime;
   info.msg = msg;
%%********************************************************************