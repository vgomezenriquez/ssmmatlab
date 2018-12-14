function [str, ierror] = invmodel(str)
% PURPOSE: given a structure with a noninvertible model, it obtains the invertible model
% with the same covariance generating function than the original model.
% The model is assumed to be in echelon form. The invertible model is also
% in echelon form. However, if there are additional constraints in the
% original model, the invertible model does not impose those constraints.
% It only imposes the constraints of the echelon form.
%---------------------------------------------------
% USAGE: str=invmodel(str)
% where:    str    = a structure containing the model information
%---------------------------------------------------
% RETURNS: str  = the structure with the invertible model
%---------------------------------------------------
% The method is based on factorizing the moving average part. That is,
% given the covariance generating function
%      G(z) = Theta(z)*Sigma*Theta(z^{-1})',
% where Theta(z) is noninvertible, an invertible matrix polynomial Thetab(z) and a
% new positive definite matrix Sigmab are computed such that
%      G(z) = Thetab(z)*Sigmab*Thetab(z^{-1})'.
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

str.noninv = 1;
ierror = 0;


% obtain the invertible model using the DARE
stre.F = str.Fs;
stre.H = str.Hs;
stre.K = str.Ks;
stre.Sigma = str.sigmar2;
stre.kro = str.kro;
stre.npar = str.npar;
stre.s = str.s;
stre.m = str.m;


stre.phi = str.phis;
stre.theta = str.thetas;
stre.gamma = str.gammas;


%obtain the invertible process using polynomial methods
%use polynomial methods to solve the covariance factorization problem
[stre, ferror] = enfinvp(stre);
if (ferror > 0)
    disp('enforcing invertibility did not work in invmodel')
    ierror = 1;
    return
end

%  maxmodmac=max(abs(eig(str.Fs-stre.Kb*str.Hs))); str.maxmodmac=maxmodmac;

% change estimates
str.sigmar2 = stre.Sigmab;
str.thetas = stre.thetab;
str.Ks = stre.Kb;


stre.phi = str.phis;
stre.theta = str.thetas;
stre.gamma = str.gammas;
stre.s = str.s;
stre.m = str.m;
stre.kro = str.kro;
stre = vecparwr(stre); %fill in vgam with new model parameters
neqs = str.s;
vgamsc = stre.vgam;
vgams = str.vgams;
vgams(1:end-neqs) = vgamsc(1:end-neqs);
%there may be more nonzero paramters in vgams than in vgam due to the
%invertibility of the model. This happens when we first make zero some
%echelon form parameters and then we invert the model.

bind = str.bind;
[nbind, mbind] = size(bind);
beta = zeros(size(bind));
for i = 1:nbind
    beta(i) = vgams(bind(i));
end

str.vgams = vgams;
str.beta = beta;
