function [Omega, Theta, ierror, iter, normdif] = pmspectfac(Lp, niter, lim)
%
% This function computes the spectral factorization
%   Lp(z,z^{-1})=Theta(z)*Omega*Theta'(z^{-1}),
% where Theta(0)=I.
%This is achieved by solving the similar problem
%   Lp(z,z^{-1})=phi(z)*phi'(z^{-1}),
% where phi is a polynomial matrix,
%   phi(z)=phi_0+phi_1*z+ ... + phi_p*z^p,
% such that phi_0 is a lower triangular matrix, Lp is a symmetric Laurent
% polynomial matrix of the form
%   Lp(z,z^{-1})= L'_pz^{-p}+.....+L'_1z^{-1}+L_0+L_1z+...L_pz^p,
% and L_0 is a symmetric, positive definite, matrix.
% The solution is found using Newton's method, iterating in the
% symmetric polynomial matrix equation
%    X(z)phi'(z^{-1}) + phi(z)X'(z^{-1}) = 2*Lp(z,z^{-1}),
% where p=degree(phi(z))= degree(X(z)), and the starting value for phi
% is phi_0=C*C', where C is a lower triangular matrix such that C*C'=L_0,
% and phi_i=0, i=1,...,p.
%
% Input:  Lp  = [n,n,p+1] matrix containing L_0+L_1z+...L_pz^p
%       niter = maximum number of iterations (default 10)
%        lim  = coefficient error limit for convergence
%              max(|phi(z)*phi'(z^{-1})-Lp|) < lim (default 1e-6)
% Output: Theta = [n,n,p] matrix containing the solution Theta(z) without
%                 Theta(0)=I
%         Omega = n x n symmetric, positive definite
%        ierror = 0,1  a flag for errors in dimensions
%          iter = number of needed iterations
%       normdif = norm of the difference upon convergence
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
%
[nl, ml, pl] = size(Lp);
Omega = [];
Theta = zeros(nl, nl, pl-1);
ierror = 0;
if (nl ~= ml)
    disp('wrong dimensions of Lp in pmspectfac');
    ierror = 1;
    return
end
if nargin <= 2
    lim = 1e-6;
    if nargin == 1
        niter = 10;
    end
end;

s = nl;
p = pl - 1;
X = zeros(size(Lp));
phi = zeros(size(Lp));
phi1 = phi;
Lp0 = Lp(:, :, 1);
Gme = 2 * Lp;
[R, r] = chol(Lp0);
C = R';
phi(:, :, 1) = C;

k = 0;
da = 1;
while (da > lim)
    k = k + 1;
    X = sympmeq(phi, Gme); %X(:,:,1) is constrained to by symmetric.
    phi1 = .5 * (phi + X);
    dphi = phi1 - phi;
    da = 0;
    for i = 1:pl
        da = da + norm(dphi(:, :, i), 1);
    end
    %  k
    %  da
    phi = phi1;
    if (k == niter)
        ierror = 2;
        break
    end
end
% k,phi

iter = k;
normdif = da;
C = phi(:, :, 1);
Omega = C * C';
for i = 2:pl
    Theta(:, :, i-1) = phi(:, :, i) / C;
end
