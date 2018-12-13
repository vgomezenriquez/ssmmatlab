function [order,kro] = varmaxgenid(y,x,seas,maxorder,hr3,ct,prt)
% PURPOSE: identifies a VARMAX model for the series y with inputs x.
% According to the sign of maxorder, it identifies: 
% a) if maxorder > 0, a VARMAX model using the generic neighborhood and ct 
%    = 'AIC', 'BIC' or 'HQ'.
% b) if maxorder = 0, a VARMAX(p,p,p) model using sequential LR tests. 
% c) if maxorder < 0, a VARMAX(p,p,p) model using ct = 'AIC', 'BIC' or
%    'HQ'.
%---------------------------------------------------
% USAGE: [order,kro] = varmaxgenid(y,x,seas,maxorder,hr3,ct,prt)
% where:    y      = an (nobs x neqs) matrix of y-vectors
%           x      = matrix of input variables (nobs x nx)
%                   (NOTE: constant vector automatically included)
%        seas      = seasonality
%    maxorder      = >0, maximum order to estimate the McMillan degree 
%                    using the generic neighborhood and AIC, BIC, HQ.
%                  =  0, estimate the order of the optimum VARMAX(p,p,p)
%                    model by LR tests with maximum order obtained by the 
%                    program. In this case, input ct is ignored (it can be
%                    set to empty for example).
%                    <0, -maxorder is the maximum order to estimate the 
%                    McMillan degree using equal K.i. and AIC, BIC, HQ.
%          hr3      = 1 perform only the first two stages of the HR method
%                    0 perform the three stages of the HR method
%          ct      = 'AIC','BIC', 'HQ'
%         prt      = 1 print results of the VARX, VARMAX(p,p,p) tests
%                    and different criteria
%---------------------------------------------------
% RETURNS: order = the system order, that is, the McMillan degree
%            kro = the Kronecker indices
%---------------------------------------------------
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

[ny,s] = size(y);
if ~isempty(x)
 [nx,m] = size(x);
 if (nx ~= ny)
  error('varmaxgenid: nobs in x-matrix not the same as y-matrix');
 end
else
 m=0;
end
order=[]; kro=[];
if nargin < 6
 prt=0;
end

if maxorder == 0
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 %Estimate p in VARMAX(p,p,p) by LR using the first two steps of the HR
 %method.
 %First, determination of optimum VARX length
 if (m > 0)  
  minp=5; maxp=8; 
 else
  minp=3; maxp=8; 
 end
 if seas > 1
  a=2.; 
 else
  a=1.25; 
 end
 pt=max(maxp,seas+minp); 
 minlags=0; 
 maxlags=ceil(max(log(ny)^a,pt)); %VARX length is given by this formula
 if m == 0
  [lagsopt,initres] = lratiocr(y,maxlags,minlags,prt);
  res = var_res(y,lagsopt);
 else
  [lagsopt,initres] = lratiocrx(y,maxlags,minlags,prt,x);
  res = varx_res(y,lagsopt,x);
 end
 residv=[initres(1:lagsopt,:); res]; %residuals of the VARX approximation
 if prt == 1
  fprintf(1,'lagsopt in VARX = %2d\n\n',lagsopt);
 end
 %Then, determination of p
 %%%%%%%%%%%%%%%%%%%%%%
 if (seas > 1)
  maxlagsax=min(max(lagsopt,seas),minp+seas); 
  minlags=maxlagsax-seas;   
  incr=seas;
  [lagsopts] = lratiocrax(y,maxlagsax,minlags,incr,prt,residv,x);
  if prt == 1
   fprintf(1,'lagsopts = %2d\n\n',lagsopts);
  end
  minlags=seas;  
 else
  minlags=0; maxlagsax=min(lagsopt,maxp); 
 end
 incr=1; resida=residv;
 [lagsopt] = lratiocrax(y,maxlagsax,minlags,incr,prt,resida,x);
 if prt == 1
  fprintf(1,'lagsopt in VARMAX(p,p,p) = %2d\n',lagsopt);
 end
 kro=repmat(lagsopt,1,s); order=sum(kro);
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
else 
 if maxorder > 0
  incrm=-1;
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 %Estimate McMillan degree using the generic neighborhood and the HR
 %method. The models are selected using  AIC, BIC or HQ criteria.
 else
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 %Estimate McMillan degree using equal K.i. and the HR
 %method. The models are selected using  AIC, BIC or HQ criteria.
  q=floor(-maxorder/s); maxorder=s*q;
  incrm=-s;   
 end
 
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 % First step of the HR method. This step is common to all models.
 % parameters for the long VARX
 if (m > 0)  
  minp=5; 
 else
  minp=3;  
 end
 a=2.0; 
 pt=max(8,seas+minp); 
 minlags=0; 
 maxlags=ceil(max(log(ny)^a,pt)); %VARX length is given by this formula
 if prt == 1
  fprintf(1,'VARX length = %2d\n',maxlags);
 end
 initres=zeros(maxlags,s);
 for i=minlags:maxlags-1
  if m == 0
   resid = var_res(y,i); initres(i+1,:)=resid(1,:);
  else
   resid = varx_res(y,i,x); initres(i+1,:)=resid(1,:);
  end
 end
 if m == 0
  res = var_res(y,maxlags);  
 else
  res = varx_res(y,maxlags,x);  
 end
 residv=[initres(1:maxlags,:); res];  % VARX residuals
%  %alternative to fixed VARX order
%  if m == 0
%   [lagsopt,initres] = lratiocr(y,maxlags,minlags,0);
%   res = var_res(y,lagsopt);
%  else
%   [lagsopt,initres] = lratiocrx(y,maxlags,minlags,0,x);
%   res = varx_res(y,lagsopt,x);
%  end
%  residv=[initres(1:lagsopt,:); res]; % VARX residuals
%  %end of alternative 
 % str.sigmarv=cov(res,1);  %felix uses this, not all of the residuals
 sigmarv=cov(residv,1); %this covariance matrix will be used in hanris2
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


if strcmp(ct,'bic') == 1 %ct == 'bic'
 cmin0=0.01; cmin1=0.;% crit=dlny;  
elseif strcmp(ct,'aic') == 1 %ct == 'aic'
 cmin0=0.; cmin1=0.; %  crit=2.d0;  
elseif strcmp(ct,'hq ') == 1 %ct == 'hq '
 cmin0=0.; cmin1=0.; %  crit=dllny; 
end
 bicm=1.d10;
 %loop over all McMillan degrees
 for n = maxorder:incrm:0
  q=floor(n/s); r=n-s*q; gkro=zeros(1,s); %generic neighborhood
  for i=1:s
   if i <= r
    gkro(i)=q+1;
   else
    gkro(i)=q;
   end
  end
  %Second and, possibly, third steps of the Hannan-Rissanen method
  str = matechelon(gkro,s,m); maxgkro=max(gkro);
  if (seas > 1) && (maxgkro >= seas)
   for j=minp+2:seas
    for ii=1:s
     for jj=1:s
      str.phi(ii,jj,j)=0; str.theta(ii,jj,j)=0;
      str.nparm=str.nparm-2;
     end
    end
   end
  end
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 % Second step of HR method:
  str.sigmarv=sigmarv;  str = vecparwr(str); 
  str = mhanris2(y,residv,x,str);
%   %invert model if noninvertible
%   noninv2=str.noninv2;  
%   if (noninv2 > 0 )  
%    str = inv2(str);
%   end
%   %end of invert model
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 % Third step of HR method: only if model is invertible and hr3 = 0
  if (str.noninv2 == 0) && (hr3 == 0)
   str = mhanris3(y,x,str); 
  end
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  noninv=0; 
  if (str.noninv2 > 0) || (isfield(str,'noninv3') && str.noninv3 == 1)
   noninv=1;
  end
  nonst=0; 
  if (str.nonst2 > 0) || (isfield(str,'nonst3') && str.nonst3 == 1)
   nonst=1;
  end
  if (noninv + nonst == 0) %only invertible and stationary models are considered
   %compute criterion statistic
   bic = mcbic(y,x,str,ct);
   rbic=abs(1-bic/bicm);
   if (bic < bicm)
    if ((n == 0) && (rbic > cmin0)) ||((n== 1) && (rbic > cmin1)) || ( n > 1)
      bicm=bic; order=n; kro=gkro;
    end
   end
   if prt > 0
    out=[n,bic,bicm,rbic,order]; 
    fprintf(1,['n = %2d, bic = %6.4g, bicm = %6.4g, '....
               'rbic = %10.4g, order = %2d,'],out);
    fprintf(1,' kro = '),fprintf(1,'%2d ',kro),fprintf(1,'\n');   end
  end
 end
end
