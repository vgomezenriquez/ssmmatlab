function [H, F, G, J, ierror] = qarmax2ss1(phi, theta)
%
%---------------------------------------------------
% This function puts the armax model into Akaike's state space form
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

% USAGE: [H,F,G,J,ierror] = qarmax2ss(phi,theta)
% where:    phi   = a matrix polynomial with phi(0) nonsingular
%           theta = a matrix polynomial not necessarily square
%---------------------------------------------------
% RETURNS:
%           H = a k x n matrix
%           F = an n x n matrix
%           G = a n x m matrix
%           J = a k x m matrix
%        ierror =1, dimension mismatch in phi and theta
%               =0, there are no errors on input
%---------------------------------------------------
%
% Copyright (c) 21 July 2003 by Victor Gomez
% Ministerio de Hacienda, Direccion Gral. de Presupuestos,
% Subdireccion Gral. de Analisis y P.E.,
% Alberto Alcocer 2, 1-P, D-34, 28046, Madrid, SPAIN.
% Phone : +34-915835439
% E-mail: vgomez@sgpg.meh.es
%
% The author assumes no responsibility for errors or damage resulting from the use
% of this code. Usage of this code in applications and/or alterations of it should
% be referenced. This code may be redistributed if nothing has been added or
% removed and no money is charged. Positive or negative feedback would be appreciated.
%
H = [];
F = [];
G = [];
J = [];
ierror = 0;
[s1, s2, np] = size(phi);
[s3, s4, nq] = size(theta);
if (s1 ~= s2) | (s1 ~= s3)
    disp('wrong dimensions of phi or theta in qarmax2ss1');
    ierror = 1;
    return
end
p = np - 1;
q = nq - 1;
r = max(p, q);
[Psi, ierror] = ptransfer(phi, theta, r+1);
P0 = phi(:, :, 1);
if any(any(P0-eye(s1)))
    for i = 2:np
        phi(:, :, i) = P0 \ phi(:, :, i);
    end
    phi(:, :, 1) = eye(s1);
end
n = s1 * r;
m = s4;
k = s1;
J = Psi(:, :, 1);
F = zeros(n);
H = zeros(k, n);
G = zeros(n, m);
H(1:k, 1:k) = eye(k);
bphi = zeros(s1, s2, r+1);
bphi(:, :, 1:np) = phi;
for i = 1:r
    G((i - 1)*k+1:i*k, :) = Psi(:, :, i+1);
    F(end-k+1:end, (i - 1)*k+1:i*k) = -bphi(:, :, r-i+2);
end
F(1:end-k, k+1:end) = eye(n-k);
