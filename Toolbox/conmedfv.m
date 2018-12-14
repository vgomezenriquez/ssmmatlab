function [residv, beta, str] = conmedfv(beta, y, x, str)
% PURPOSE: given a structure, it computes the model residuals and their
% covariance matrix using the third step parameters of HR method.
%---------------------------------------------------
% USAGE: [residv,beta,str]=conmedfv(beta,y,x,str)
% where:   beta    = an (1 x nparm) vector of parameters
%           y      = an (nobs x neqs) matrix of y-vectors
%           x      = matrix of input variables (nobs x nx)
%           str    = a structure containing the model information
%---------------------------------------------------
% RETURNS: residv  = the residuals
%           beta   = an [(nparm + s) x 1] vector containing the parameters,
%                    where
%                    nparm : number of parameters
%                        s : number of ouputs
%            str   = updated structure containing the model information
%                    More specifically, if the model corresponding to str
%                    that is input to conmedfv is nonstationary or
%                    noninvertivble, an appropriate transformation is made
%                    to the AR or the MA part of the model to make them
%                    stable.
%                    The transformed parameters are in the updated field
%                    vgams.
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

kro = str.kro;
nlag = max(kro);
vgams = str.vgam;
bind = str.bind;
[nbind, mbind] = size(bind);
for i = 1:nbind
    vgams(bind(i)) = beta(i);
end
str.vgams = vgams;
str = param2sse(str);


%check stationarity and invertibility. If necessary, change parameters.
iar = chkstainv(str.Fs); %if iar >1, the model is not stationary
if iar > 1
    %  fprintf(1,'model nonstationary, iar = %2d\n',iar);
    %convert Phi(z) into Phi(lambda*z) for an appropriate lambda
    vgam = enfstab(str, 'phi  ');
    str.vgams = vgam;
    str = param2sse(str);
    for i = 1:nbind
        beta(i) = vgam(bind(i));
    end
end
ima = chkstainv(str.Fs-str.Ks*str.Hs); %if ima >1, the model is not invertible
if ima > 1
    %  fprintf(1,'model noninvertible, ima = %2d\n',ima);
    %convert Th(z) into Th(lambda*z) for an appropriate lambda
    vgam = enfstab(str, 'theta');
    str.vgams = vgam;
    str = param2sse(str);
    for i = 1:nbind
        beta(i) = vgam(bind(i));
    end
end

[resid, sigmar] = compresde0(y, x, str); % constant is passed in str.mu
% str.sigmar23=sigmar;  %store covariance matrix of residuals
[R, p] = chol(sigmar);
L = R';
if (p > 0)
    error('covariance matrix of residuals2 singular in conmedfv')
end
residv = vec(L \ resid(nlag+1:end, :)');
