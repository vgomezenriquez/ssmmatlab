function oparm=tradid(y,Y,infm,parm,ser,ols,a,tradval,fid,fmarqdt)
%
% this function automatically estimates the number of TD variables for an 
% ARIMA model with TD correction
%
% Input arguments:
% y: vector containing the data
% Y: matrix containing regression variables
%   infm     : structure containing function names and optimization options
%   .f  :   a function to evaluate the vector ff of individual functions
%           such that ff'*ff is minimized
%   .tr :   >0 x is passed from marqdt to f but not passed from f to marqdt  
%           =0 x is passed from marqdt to f and passed from f to marqdt 
%   .tol:   a parameter used for stopping
%   .jac:   =1 evaluation of jacobian and gradient at the solution is performed
%           =0 no evaluation of jacobian and gradient at the solution is performed
% .maxit:   maximum number of iterations
%   .nu0:   initial value of the nu parameter
%   .prt:   =1 printing of results
%           =0 no printing of results
% parm: astructure containing model information, where
% .s:  seasonality
% .S:  second seasonality
% .p:  AR order
% .ps: order of the AR of order s
% .q:  order of the regular MA
% .qs: order of the MA of order s (1 at most)
% .qS: order of the MA of order S (1 at most)
% .dr: order of regular differencing
% .ds: order of differencing of order s
% .dS: order of differencing of order S
% .pvar:  array containing the indices of variable parameters
% .pfix:  array containing the indices of fixed parameters
% ser  : a structure containing the series parameters (the ones specified
%        by the user in the spec file and the default ones)
% ols  : = 1, perform OLS, = 0, use the Durbin Levinson algorithm in the HR
%        method
% a    : an integer, the degree of the AR approximation in the first step
%        of the Hanna-Rissanen method.
% tradval: an integer array containing the possible numbers of TD variables 
%         (0 is also a value)
% fid    : the number of the external output file 
% fmarqdt: a parameter for the estimation method
%          = 1 Levenberg-Marquardt method
%          = 0 Lsqnonlin (Matlab)
%
% Output arguments:
% oparm: the input parm structure plus the field
% .trad : the estimated number of TD variables
%
% Copyright (c) 21 July 2003 by Victor Gomez
% Ministerio de Hacienda, Direccion Gral. de Presupuestos, 
% Subdireccion Gral. de Analisis y P.E., 
% Alberto Alcocer 2, 1-P, D-34, 28046, Madrid, SPAIN.
% Phone : +34-915835439
% E-mail: VGomez@sepg.minhap.es

%
% The author assumes no responsibility for errors or damage resulting from the use
% of this code. Usage of this code in applications and/or alterations of it should 
% be referenced. This code may be redistributed if nothing has been added or
% removed and no money is charged. Positive or negative feedback would be appreciated.
%

if ~isstruct(parm)
   error('tradid: requires a parameter structure');
end;
oparm=parm; ny=length(y);
s=parm.s; S=parm.S; dr=parm.dr; ds=parm.ds; dS=parm.dS; pvar=parm.pvar; pfix=parm.pfix;
p=oparm.p; q=oparm.q; ps=oparm.ps; qs=oparm.qs; qS=oparm.qS;
bg_year=ser.bg_year; bg_per=ser.bg_per; freq=ser.freq;
bicm=1.d10; ntrad=length(tradval);
for i=1:ntrad
 trad=tradval(i);
 if trad > 0 
  Ytr=trade(bg_year,bg_per,ny,trad,0,[],freq); YY=[Y Ytr]; 
 else 
  YY=Y; 
 end
 est=1; x0 = cinest(y,YY,parm,est,ols,a,0,fid); xv=x0(pvar); xf=x0(pfix); 
 if ~isempty(xv)
  xv=arimaopt(fmarqdt,fid,x0,xv,xf,y,YY,parm,infm,0);
  x0(pvar)=xv;
 %  [F,e,beta,M]=residual2(x0,y,YY,s,dr,ds,p,ps,q,qs,1); nyd=ny-dr-ds*s; 
 end
 [F,e,beta,M] = residual2x(x0,y,YY,s,S,dr,ds,dS,p,ps,q,qs,qS);
 nyd=ny-dr-ds*s-dS*S;
 nbeta=length(beta);
 dn=double(nyd); var=e'*e/dn; ldn=log(dn);  bic=log(var)+double(nbeta)*ldn/dn;
	if (bic < bicm) 
   bicm=bic;   oparm.trad=trad;  
	end
end
