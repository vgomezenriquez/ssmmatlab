function [str, ierror] = enfinvp(str)
% PURPOSE: enforces invertibility in a VARMAX model in echelon form
%          using polynomial methods
%---------------------------------------------------
% USAGE:  str = enfinv(str)
% where:    str is a structure containing all the information about the
%           model
%---------------------------------------------------
% RETURNS: the same structure with the invertible model
%---------------------------------------------------
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
phi = str.phi;
theta = str.theta;
gamma = str.gamma;
Sigma = str.Sigma;
kro = str.kro;
s = str.s;

p = max(kro);
pp1 = p + 1;
thetab = zeros(size(theta));


% ths=theta;
% if (sum(sum(Sigma - eye(s))) ~= 0)
%  for i=1:pp1
%   ths(:,:,i)=theta(:,:,i)*Sigma;
%  end
% end
% G1=pmmulbf(ths,theta); %generating function of MA part in echelon form


% left divide theta(z) by theta(0)
A = theta(:, :, 1);
if any(any(A-eye(s)))
    theta(:, :, 1) = eye(s);
    for i = 2:pp1
        theta(:, :, i) = A \ theta(:, :, i);
    end
end

ths = theta;
if (sum(sum(Sigma-eye(s))) ~= 0)
    for i = 1:pp1
        ths(:, :, i) = theta(:, :, i) * Sigma;
    end
end
% ths
% theta
G = pmmulbf(ths, theta); %generating function of MA part
Gm = G(:, :, pp1:end);
[Sigmab, Theta, ierror, iter, normdif] = pmspectfac(Gm, 20, 1e-3);
if (ierror > 0)
    disp('enforcing invertibility did not work in enfinvp')
    return
end
% S1=Sigmab


% % solve the DARE for comparison
% H=str.H; F=str.F; K=str.K;
% P=msdare(F,K,H,Sigma,Sigma,Sigma);
% Sigmab=H*P*H'+Sigma;
% Kb=(F*P*H'+K*Sigma)/Sigmab;
% S2=Sigmab  %it should be equal to S1


%adjustment for the echelon form because Theta(0)=I
thetab(:, :, 1) = A;
for i = 2:pp1
    thetab(:, :, i) = A * Theta(:, :, i-1);
    %enforce echelon form constraints. Roots may change
    %  for j=1:s
    %   for k=1:s
    %    if (str.theta(j,k,i) == 0)
    %     thetab(j,k,i)=0;
    %    end
    %   end
    %  end
end


% ths=thetab;
% for i=1:pp1
%  ths(:,:,i)=thetab(:,:,i)*Sigmab;
% end
% G2=pmmulbf(ths,thetab); %new generating function of MA part. It should be
% %equal to G1. Thus, the following should be zero
% G1-G2
% % Theta
% % thetab


stre.phis = phi;
stre.thetas = thetab;
stre.gammas = gamma;
stre.npar = str.npar;
stre.s = str.s;
stre.kro = str.kro;
stre.m = str.m;
stre = armaxe2sse(stre);
ima = chkstainv(stre.Fs-stre.Ks*stre.Hs); %if ima > 0 the model is noninvertible
if ima > 0
    disp('enforcing invertibility did not work in enfinvp')
    ierror = 1;
    return
end

Kb = stre.Ks;

str.Kb = Kb;
str.Sigmab = Sigmab;
str.thetab = thetab;
% Kb
% Sigmab
% thetab
