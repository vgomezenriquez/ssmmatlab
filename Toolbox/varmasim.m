function [z, ferror] = varmasim(l, N, phi, th, stda, seed)
%
%
%        This function generates a VARMA model
%        It uses Akaike's state space representation
%
%        Input parameters:
%        l: number of observations discarded at the beginning of the series
%        N: number of observations of the simulated series
%        phi: AR matrix polynomial
%        th:  MA matrix polynomial
%        stda:  covariance matrix of the innovations  (default: I)
%        seed: a number to start random normal generation
%
%        Output parameters:
%        z: the simulated series
%        ferror: a flag for errors
%
% Copyright (c) 21 December 2010 by Victor Gomez
% Ministerio de Hacienda, Direccion Gral. de Presupuestos,
% Subdireccion Gral. de Analisis y P.E.,
% Alberto Alcocer 2, 1-P, D-34, 28046, Madrid, SPAIN.
% Phone : +34-915835439
% E-mail: vgomez@sgpg.minhac.es
%
% The author assumes no responsibility for errors or damage resulting from the use
% of this code. Usage of this code in applications and/or alterations of it should
% be referenced. This code may be redistributed if nothing has been added or
% removed and no money is charged. Positive or negative feedback would be appreciated.
%

z = [];
ferror = 0;

[np, mp, kp] = size(phi);
[nt, mt, kt] = size(th);
if (np ~= mp) | (nt ~= mt) | (np ~= nt)
    disp('dimension mismatch in varmasim')
    ferror = 1;
    return
end

farg = 0;
if nargin < 4
    disp('number of arguments must be at least 4 in varmasim')
    ferror = 2;
    return
elseif nargin == 4
    stda = eye(np);
    farg = 4;
elseif nargin == 5
    R = chol(stda);
    stda = R';
    farg = 5;
elseif nargin == 6
    R = chol(stda);
    stda = R';
end


% Function qarmax2ss1 puts the armax model into Akaike's state space form
%   x(t+1) = F*x(t) + G*u(t)
%    y(t)  = H*x(t) + J*u(t),
% where
%            [0 I 0  ... ....   0]        [ Psi_1    ]
%            [0 0 I  ... ....   0]        [ Psi_2    ]
%     F =    [ ...   ...  ...    ],   G = [ ...      ],
%            [0 0 0  ... ....   I]        [ Psi_{r-1}]
%            [-bphi_r ... -bphi_1]        [ Psi_r    ]
%     H =    [I 0 0  ...  ... 0],   J =  Psi_0,
% bphi_i = phi_0^{-1}*phi_i and phi^{1}(z)*theta(z) = Psi_0 + Psi_1*z
%  + Psi_2*z^2+ ..., and r = max{degree(phip), degree(theta)}.
%
[H, F, K, J, ierror] = qarmax2ss1(phi, th);

% compute initial state vector
[A, Sigma, Xi] = vincovma(F, K, stda);

if (farg == 0)
    z = varimass(l, N, H, F, K, A, Sigma, Xi, stda, seed);
elseif (farg >= 4)
    z = varimass(l, N, H, F, K, A, Sigma, Xi, stda);
end