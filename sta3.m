function [str, ierror] = sta3(str)
% PURPOSE: given a structure passed after executing the third step of HR method such
% that the model is not stationary, this function makes it stationary.
%---------------------------------------------------
% USAGE: str = sta3(str)
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
stre.thetas = str.thetas3;
stre.sigmar2 = str.sigmar3;
stre.vgams = str.vgams3;
stre.phis = str.phis3;
stre.gammas = str.gammas3;
%make the model stationary
[stre, ierror] = stamodel(stre);
if ierror > 0
    disp('Polynomial P(z) is transformed into P(lambda*z)')
    stre = sta3r(stre);
end
str.phis3 = stre.phis;
str.sigmar3 = stre.sigmar2;
str.nonst3 = 0;
str.Fs3 = stre.Fs;
str.Hs3 = stre.Hs;
str.Ks3 = stre.Ks;
% insert new parameters in beta3 and vgams3
vgams = stre.vgams;
bind = str.bind;
[nbind, mbind] = size(bind);
beta = zeros(nbind, 1);
for i = 1:nbind
    beta(i) = vgams(bind(i));
end
str.beta3 = beta;
str.vgams3 = vgams;
%obtain new armax parameters (phi(1)=I, theta(1)=I).
[phist3, thetast3, gammast3, ierror] = armaxe2armax(str.phis3, ...
    str.thetas3, str.gammas3);
str.phist3 = phist3;
str.thetast3 = thetast3;
str.gammast3 = gammast3;
