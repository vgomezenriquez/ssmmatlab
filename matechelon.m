function [str, ferror] = matechelon(kro, s, m)
% PURPOSE: obtains the matrix structure of a VARMAX model in echelon form,
%          both in polynomial and state space form. The state space echelon
%          form is:
%
%  alpha_{t+1} = F*alpha_{t} + B*x_t{t} + K*a_{t}
%      y_{t}   = H*alpha_{t} + D*x_{t}  + a_{t}
%
% The VARMAX echelon form is

%  phi(B)*y_{t} = gamma(B)*x_{t} + theta(B)a_{t},
%
% where B is the backshift operator, B*y_{t} = y_{t-1}.
%--------------------------------------------------------------------------
% USAGE:  str = matechelon(kro,s,m)
% where:  kro = a (1 x s) vector containing the Kronecker indices
%           s = number of outputs
%           m = number of inputs
%--------------------------------------------------------------------------
% RETURNS: a structure with the following fields
%         s: number of outputs
%         m: number of inputs
%       kro: a (1 x s) vector containing the Kronecker indices
%       phi: an (s x s x maxkro) array
%     theta: an (s x s x maxkro) array
%     gamma: an (s x m x maxkro) array
%     nparm: number of parameters
%      npar: an (s x s) array
%         F: an (n x n) matrix
%         H: an (s x n) matrix
%         K: an (n x s) matrix
%         B: an (n x m) matrix
%         D: an (s x m) matrix
% where maxkro = max(Kronecker indices), n=McMillan degree = sum(Kronecker
% indices), npar is an array of indices necessary to define the Kronecker
% indices, phi,theta and gamma are the VARMAX matrices in echelon form and
% F, H, K, B and D are the VARMAX matrices in state space echelon form.
% The parameters in the echelon form are indicated with NaN.
%--------------------------------------------------------------------------
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
[s1, s2] = size(kro);
ferror = 0;
if (s ~= s1 && s ~= s2)
    ferror = 1;
    disp('wrong kro vector in matechelon');
    return
end
npar = zeros(s, s);
for i = 1:s
    for j = 1:i - 1
        npar(i, j) = min(kro(i)+1, kro(j));
    end
    for j = i:s
        npar(i, j) = min(kro(i), kro(j));
    end
end
% npar
mkro = max(kro);
phi = zeros(s, s, mkro+1);
theta = zeros(s, s, mkro+1);
gamma = [];
phi(:, :, 1) = eye(s);
nparm = 0;
for i = 1:s
    for p = 1:i - 1
        for k = kro(i) + 2 - npar(i, p):kro(i) + 1
            phi(i, p, k) = NaN;
            nparm = nparm + 1;
        end
    end
    for k = 2:kro(i) + 1
        phi(i, i, k) = NaN;
        nparm = nparm + 1;
    end
    for p = i + 1:s
        for k = kro(i) + 2 - npar(i, p):kro(i) + 1
            phi(i, p, k) = NaN;
            nparm = nparm + 1;
        end
    end
    for p = 1:s
        for k = 2:kro(i) + 1
            theta(i, p, k) = NaN;
            nparm = nparm + 1;
        end
    end
end
theta(:, :, 1) = phi(:, :, 1);
if (m > 0)
    gamma = zeros(s, m, mkro+1);
    for i = 1:s
        for p = 1:m
            for k = 1:kro(i) + 1
                gamma(i, p, k) = NaN;
                nparm = nparm + 1;
            end
        end
    end
end
str.s = s;
str.m = m;
str.kro = kro;
str.phi = phi;
str.theta = theta;
str.gamma = gamma;
str.nparm = nparm;
str.npar = npar;

% state space echelon form
n = sum(kro); % McMillan degree
F = zeros(n, n);
H = zeros(s, n);
K = NaN(n, s); % system matrices
B = NaN(n, m);
D = NaN(s, m);
ch = 0;
for i = 1:s
    if kro(i) > 0
        H(i, ch+1) = 1;
        F(ch+1:ch+kro(i)-1, ch+2:ch+kro(i)) = eye(kro(i)-1);
    else
        % if kro(i)=0, we should put in matrix H the coefficients of
        % the \Phi_0 matrix corresponding to the i-th variable that
        % are different from zero.
        chb = 0;
        for j = 1:s
            if npar(i, j) > 0
                H(i, chb+1) = NaN;
            end
            chb = chb + kro(j);
        end
    end
    chb = 0;
    for j = 1:s
        if npar(i, j) > 0
            F(ch+kro(i), chb+1:chb+npar(i, j)) = NaN(1, npar(i, j));
        end
        chb = chb + kro(j);
    end
    ch = ch + kro(i);
end
str.F = F;
str.H = H;
str.K = K;
str.B = B;
str.D = D;
