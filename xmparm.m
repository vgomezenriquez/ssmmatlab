function [xx,pvar,pfix,parm]= xmparm(ninput,x0,xm,xf,pvar,pfix,nr,nlagtf,g,...
                                     tford,parm,nreg)
%
% this function automatically identifies the input filters in a transfer
% function model using the LTF method. See Liu, L. M., and Hanssens, D. M. 
% (1982), "Identification of Multiple{Input Transfer Function Models", 
% Communications in Statistics, Theory and Methods, 11, 297-314. 
%
% Input arguments:
% x0    : array containing the initial model parameters
% xm    : array containing the ARMA variable parameters
% xf    : array containing the ARMA fixed parameters
% pvar  : array containing the indices of ARMA variable parameters
% pfix  : array containing the indices of ARMA fixed parameters
% nr    : the number of ARMA variable parameters
% nlagtf: the number of lags for the polynomial approximations to the 
%         rational input filter expansions (see LTF method)
% g     : array containing the estimated weigths of the polynomial
%         approximations to the rational input filter expansions (see LTF
%         method)
% tford : (ninput x 3) array containing for each input variable the delay,
%         the degree of the ma part, and the degree of the ar part
% parm : a structure where
%  s:  seasonality
%  S:  second seasonality
% .p:  AR order
% .ps: order of the AR of order s
% .q:  order of the regular MA
% .qs: order of the MA of order s (1 at most)
% .qS: order of the MA of order S (1 at most)
% .dr: order of regular differencing
% .ds: order of differencing of order s
% .dS: order of differencing of order S
% nreg: number of regression variables
%
% Output arguments:
% xx  : an array containing the transfer function model parameters
% pvar: an array containing the transfer function model variable parameters
% pfix: an array containing the transfer function model fixed parameters
% parm: a structure containing the transfer function model information. In
%       addition to the input fields, it contains
% .pvar:  array containing the indices of variable parameters
% .pfix:  array containing the indices of fixed parameters
% .ninput: number of inputs
% .delay: array with the delays of the input filters
% .ma: array with the ma parameters of the input filters  
% .ar: array with the ar parameters of the input filters  
%
%
% Copyright (c) 21 July 2015 by Victor Gomez
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

parm.ninput = ninput;
parm.pvar = []; %indices should start after nr (number of arima parameters)
parm.pfix = []; 
xx=x0; xx(pvar)=xm; xx(pfix)=xf; org=nr;  
if ninput > 0
 parm.delay = zeros(ninput,1);
 parm.ma = zeros(ninput,1);
 parm.ar = zeros(ninput,1);  
 for i=1:ninput
  [bb,aa,err]=shank(g(nreg+(nlagtf+1)*(i-1)+tford(i,1)+1:nreg+(nlagtf+1)*i),...
                tford(i,2),tford(i,3));
  parm.delay(i)=tford(i,1);
  parm.ma(i)=tford(i,2);
  parm.ar(i)=tford(i,3);
  parm.pvar=[parm.pvar org+1:org+tford(i,2)+1+tford(i,3)];
  org=org+tford(i,2)+1+tford(i,3);
  xx=[xx,bb'];
  if tford(i,3) > 0
   xx=[xx,aa(2:end)'];
  end
 end
 pvar=[pvar parm.pvar]; parm.pvar=pvar;
 pfix=[pfix parm.pfix]; parm.pfix=pfix; 
end
