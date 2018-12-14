function str = sta2(str)
% PURPOSE: given a structure passed after executing the second step of HR method such
% that the model is not stationary, this function makes it stationary.
%---------------------------------------------------
% USAGE: str = sta2(str)
% where:    str    = a structure containing the structure of the VARMAX
%---------------------------------------------------
% RETURNS: str = a structure containing the inverted model
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

stre = str;
stre.thetas = str.thetas;
stre.sigmar2 = str.sigmar2;
stre.vgams = str.vgams;
stre.phis = str.phis;
stre.gammas = str.gammas;
%make the model stationary
[stre, ierror] = stamodel(stre);
if ierror > 0
    disp('Polynomial P(z) is transformed into P(lambda*z)')
    stre = sta2r(stre);
end
str.phis = stre.phis;
str.sigmar2 = stre.sigmar2;
str.nonst2 = 0;
str.Fs = stre.Fs;
str.Hs = stre.Hs;
str.Ks = stre.Ks;
% insert new parameters in beta and vgams
vgams = stre.vgams;
bind = str.bind;
[nbind, mbind] = size(bind);
beta = zeros(nbind, 1);
for i = 1:nbind
    beta(i) = vgams(bind(i));
end
str.beta = beta;
str.vgams = vgams;
%obtain new armax parameters (phi(1)=I, theta(1)=I).
[phist, thetast, gammast, ierror] = armaxe2armax(str.phis, ...
    str.thetas, str.gammas);
str.phist = phist;
str.thetast = thetast;
str.gammast = gammast;
