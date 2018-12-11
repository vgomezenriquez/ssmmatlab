function [H, F, G, J, ierror] = qarmax2ss2(phi, theta)
%
%---------------------------------------------------
% USAGE: [H,F,G,J,ierror] = qarmax2ss(phi,theta)
% where:    phi   = a matrix polynomial with phi(0) nonsingular
%           theta = a matrix polynomial not necessarily square
% This function puts the armax model into a state space form
%   x(t+1) = F*x(t) + G*u(t)
%    y(t)  = H*x(t) + J*u(t),
% where
%            [-bphi_1     I 0  ... .... 0]        [ theta_1-phi_1*Psi0        ]
%            [-bphi_2     0 I  ... .... 0]        [ theta_2-phi_2*Psi0        ]
%     F =    [  ...            ...  ...  ],   G = [       ...                 ],
%            [-bphi_{r-1} 0 0  ...  ... I]        [ theta_{r-1}-phi_{r-1}*Psi0]
%            [-bphi_r     0 0  ...  ... 0]        [ theta_r-phi_r*Psi0        ]
%     H =    [phi_0^{-1} 0 0  ...  ... 0],    J =  Psi_0,
% bphi_i = phi_i*phi_0^{-1} and phi^{-1}(z)*theta(z) = Psi_0 + Psi_1*z + Psi_2*z^2+ ...
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
% E-mail: VGomez@sepg.minhap.es
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
P0 = phi(:, :, 1);
p = np - 1;
q = nq - 1;
Psi0 = theta(:, :, 1);
invp = 0;
if any(any(P0-eye(s1)))
    invp = 1;
    Psi0 = P0 \ Psi0;
end
r = max(p, q);
n = s1 * r;
m = s4;
k = s1;
J = Psi0;
F = zeros(n);
H = zeros(k, n);
G = zeros(n, m);
if invp == 1
    H(1:k, 1:k) = pinv(P0);
else
    H(1:k, 1:k) = eye(k);
end
bphi = zeros(s1, s2, r+1);
bphi(:, :, 1:np) = phi;
btheta = zeros(s3, s4, r+1);
btheta(:, :, 1:nq) = theta;
for i = 1:r
    if invp == 1
        F((i - 1)*k+1:i*k, 1:k) = -bphi(:, :, i+1) / P0;
    else
        F((i - 1)*k+1:i*k, 1:k) = -bphi(:, :, i+1);
    end
    G((i - 1)*k+1:i*k, :) = btheta(:, :, i+1) - bphi(:, :, i+1) * Psi0;
end
F(1:end-k, k+1:end) = eye(n-k);
