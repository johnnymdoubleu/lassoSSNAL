function Ay = matvec_ClassicLasso(y,par,AP)
tmp = AP'*y;
%fprintf('\n  SAPall = %3.8f',sum(AP,'all')); %DBGGG
%fprintf('\n  APdim = %3.8f',size(AP)); %DBGGG
%fprintf('\n  SP = %3.8f',iterpsqmr); %DBGGG
%input("Stop") %DBG
Ay = y + par.sigma*(AP*tmp);