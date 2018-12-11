function str = inv2(str)
% PURPOSE: given a structure passed after executing the second step of HR method such
% that the model is not invertible, this function inverts the model using the DARE.
%---------------------------------------------------
% USAGE: str = inv2(str)
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
%invert ma part; the additional constraints are not imposed
[stre, ierror] = invmodel(stre);
if ierror > 0
    disp('Polynomial P(z) is transformed into P(lambda*z)')
    stre = str;
    stre = inv2r(stre);
    %  disp('enforcing invertibility did not work in inv3')
    %  return
end
str.thetas = stre.thetas;
str.sigmar2 = stre.sigmar2;
str.noninv2 = 0;
str.Ks = stre.Ks;
% insert new parameters in beta and vgams
str.beta = stre.beta;
str.vgams = stre.vgams;
%obtain new armax paramters (phi(1)=I, theta(1)=I).
[phist, thetast, gammast, ierror] = armaxe2armax(str.phis, ...
    str.thetas, str.gammas);
str.phist = phist;
str.thetast = thetast;
str.gammast = gammast;
