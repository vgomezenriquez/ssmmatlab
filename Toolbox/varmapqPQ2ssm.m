function [phirs, thrs, H, F, G, J, ferror] = varmapqPQ2ssm(phi, th, Phi, Th, L, str)
% PURPOSE: given the matrix polynomials of a VARMA(p,q)(P,Q)_s model and
% the Cholesky factor of the innovations covariance, this function puts the
% model into the state space form
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
%
%and computes the AR and MA matrix polynomials
%---------------------------------------------------
% USAGE: [phirs,thrs,H,F,G,J,ferror] = varmapqPQ2ssm(phi,th,Phi,Th,L,str)
% where:   phi      = the regular AR matrix polynomial
%          th       = the regular MA matrix polynomial
%          Phi      = the seasonal AR matrix polynomial
%          Th       = the seasonal MA matrix polynomial
%          L        = the Cholesky factor of the innovations covariance
%                     matrix
%          str      = a structure containing model information
%---------------------------------------------------
%---------------------------------------------------
% RETURNS: phirs = the overall AR matrix polynomial
%          thrs  = the overall MA matrix polynomial
%          H     = a matrix of the state space form
%          F     = a matrix of the state space form
%          G     = a matrix of the state space form
%          J     = a matrix of the state space form
%      ferror    = flag for errors
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

ferror = 0;
H = [];
F = [];
G = [];
J = [];

[np, mp, pm1] = size(phi);
[nt, mt, qm1] = size(th);
[nP, mP, Pm1] = size(Phi);
[nT, mT, Qm1] = size(Th);
[ns, ms] = size(L);

if ((np ~= nt) | (np ~= nP) | (np ~= nT) | (np ~= ns) | ...
        (mp ~= mt) | (mp ~= mP) | (np ~= mT) | (np ~= ms))
    ferror = 1;
    disp('dimensions of phi, th, Phi, Th and L should agree in varmapqPQ2ssm');
    return
end

%transform multiplicative to non-multiplicative model
seas = str.freq;
p = pm1 - 1;
q = qm1 - 1;
P = Pm1 - 1;
Q = Qm1 - 1;
if (seas > 1)
    prsm1 = p + seas * P + 1;
    phirs = zeros(np, mp, prsm1);
    qrsm1 = q + seas * Q + 1;
    thrs = zeros(np, mp, qrsm1);
    phirs(:, :, 1:p+1) = phi;
    thrs(:, :, 1:q+1) = th;
    for i = 1:P
        seasi = seas * i;
        phirs(:, :, seasi+1:seasi+p+1) = pmatmul(phi, Phi(:, :, i+1));
    end
    for i = 1:Q
        seasi = seas * i;
        thrs(:, :, seasi+1:seasi+q+1) = pmatmul(th, Th(:, :, i+1));
    end
else
    phirs = phi;
    thrs = th;
end

%set up state space model
[H, F, G, J, ierror] = qarmax2ss2(phirs, thrs);
G = G * L;
J = L;
