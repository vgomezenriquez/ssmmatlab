function [str, ierror] = stamodel(str)
% PURPOSE: given a structure with a nonstationary model, it obtains a
% stationary model.
% The model is assumed to be in echelon form. The stationary model is also
% in echelon form, but it imposes less constraints. The constraints are
% equal to those of the MA part.
%---------------------------------------------------
% USAGE: str=stamodel(str)
% where:    str    = a structure containing the model information
%---------------------------------------------------
% RETURNS: str  = the structure with the invertible model
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

ierror = 0;

str.nonsta = 1;


% obtain the stationary model solving the DARE for the inverse process
kro = str.kro;
s = str.s;
stre.F = str.Fs;
stre.H = str.Hs;
stre.K = str.Ks;
stre.Sigma = str.sigmar2;
stre.kro = kro;
stre.phi = str.phis;
stre.theta = str.thetas;
stre.npar = str.npar; %stre.gamma=str.gamma; stre.gammas=str.gammas;
stre.gamma = str.gammas;
stre.s = s;
stre.m = str.m;
stre.kro = str.kro; %stre.vgam=str.vgam;
%  stre.vgams=str.vgams;

%  %solving the DARE for the inverse process
%  stre = enfsta(stre);

% for i=1:1
%  n=sum(kro); %McMillan degree
%  q=floor(n/s); r=n-s*q; gkro=zeros(size(kro)); %generic neighborhood
%  for i=1:s
%   if i <= r
%    gkro(i)=q+1;
%   else
%    gkro(i)=q;
%   end
%  end
%  gkro
%  if kro == gkro
%  %use the inverse process. It may change the Kronecker indices, but it
%  %works at least for the generic neighborhood
%   [stre,ierror] = enfstap(stre);
%  else
[stre, ierror] = enfstap(stre);
if ierror > 0
    disp('Enforcing stationarity in stamodel did not work')
    ierror = 1;
    return
end
%   stre = enfstap2(stre); %use phi as an MA model
%  end
%  iar=chkstainv(stre.Fb);  %if iar >1, the model is not stationary
%  if iar == 0
%   break
%  else
% %   i
% %   [H,F,G,J,ierror] = qarmax2ss2(stre.phib,eye(str.s));
% %   abs(eig(F))
%   stre.phi=stre.phib;
%  end
% end


%  maxmodarc=max(abs(eig(stre.Fb))); str.maxmodarc=maxmodarc;

% change estimates
str.sigmar2 = stre.Sigmab;
str.phis = stre.phib;
str.thetas = stre.thetab;

stre.phi = str.phis;
stre.theta = str.thetas;
stre.gamma = str.gammas;
stre.s = str.s;
stre.m = str.m;
stre.kro = str.kro;
stre = vecparwr(stre);

neqs = str.s;
vgamsc = stre.vgam;
vgams = str.vgams;
vgams(1:end-neqs) = vgamsc(1:end-neqs);

bind = str.bind;
[nbind, mbind] = size(bind);
for i = 1:nbind
    beta(i) = vgams(bind(i));
end

str.vgams = vgams;
str = param2sse(str);
