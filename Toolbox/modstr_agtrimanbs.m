function [modelsn, H_p, J_p, Q] = modstr_agtrimanbs(xx, pfix, pvar, xf, models, Q, stcs, p, nt)
%**************************************************************************
% Auxiliary function called in agtrimanssbs_d.m.
% Reference: ''A new State Space Methodology to
% Disaggregate Multivariate Time Series'', Journal of Time Series Analysis,
% 30, 97-124, by Gomez and Aparicio, (2009).
%
%  INPUTS:
%        xx : array with parameters to be estimated
%      pfix : array with fixed parameter indices
%      pvar : array with variable parameter indices
%        xf : array with fixed parameters
%    models : structure with initial univariate model information
%         Q : orthogonal matrix needed for tranformation of the model for
%             stacked observations;
%             required input if stcs is not equal 1;
%             if stcs = 1, Q is calculated by mobsvf.m
%      stcs : if stcs = 1, staircase alhorithm is applied, where matrix Q
%             is computed by mobsvf.m
%         p : number of variables in the data matrix y
%        nt : length of y
%
%    OUTPUTS:
%   modelsn : structure with multivariate model information
%   H_p,J_p : matrices in the observation equation of the model for stacked
%             observations:
%               Y_t = H_p*alpha_t + J_p*U_t
%         Q : orthogonal matrix needed for transformation of the matrices
%             of the model for stacked observations
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
%*************************************************************************


Z = models.Z;
T = models.T;
G = models.G;
H = models.H;
X = models.X;
stord = models.stord;
conc = models.conc;
trend = models.trend;
slope = models.slope;
seas = models.seas;
irreg = models.irreg;
arp = models.arp;
nalphao = models.nalphao;

npar = length(pfix) + length(pvar);
x = zeros(1, npar);
x(pfix) = xf;
x(pvar) = xx;

xp = zeros(1, 12);
xp(stord) = x;
xp(conc) = 1; %concentrated parameter

[nh, mh] = size(H);
H(1, 1) = xp(1);
H(2, 1) = xp(2);
H(2, 2) = xp(3);
H(3, 3) = xp(4);
H(4, 3) = xp(5);
H(4, 4) = xp(6);
j = 5;
for i = 6:2:nh
    H(i, j) = xp(8);
    H(i-1, j) = xp(7);
    H(i, j+1) = xp(9);
    j = j + 2;
end
% for i=6:2:nh
%  H(i,5)=xp(8); H(i-1,5)=xp(7); H(i,6)=xp(9);
% end
G(1, mh-1) = xp(10);
G(2, mh-1) = xp(11);
G(2, mh) = xp(12);

% T
% H
% Z
% G
%
% high frequency model
%
J = G;
G = H;
[nalpha, junk] = size(G);
H = Z;
F = T;
% D=kron([1/3 1/3 1/3],eye(2));              %aggregation pattern
D = models.D;
% F
% G
% H
% J
% D
% error('matrices en modstr')
%
% stacked model
%
Fp = eye(nalpha);
G_p = [];
H_p = [];
[nj, mj] = size(J);
J_p = zeros(p*nj, p*mj);
if ~isempty(X)
    V_p = reshape(X, p, nt); % INEM la agrupamos de 3 en 3 meses matriz de 3x98
else
    V_p = [];
end
aux = J; % auxiliary variable for J_p
for k = 1:p - 1
    aux = [H * F^(k - 1) * G, aux];
end
for k = 1:p
    Fpm1 = Fp;
    Fp = Fp * F;
    G_p = [Fpm1 * G, G_p];
    H_p = [H_p; H * Fpm1];
    J_p((k - 1)*nj+1:k*nj, 1:k*mj) = aux(:, p*mj-k*mj+1:p*mj);
end
% G_p
% H_p
% J_p
% error('matrices en modstr')
H_pbar = D * H_p;
J_pbar = D * J_p;
if ~isempty(V_p)
    V_pbar = (D * V_p)';
    X = V_pbar;
end
% Fp
% G_p
% H_p
% V_pbar
%
% staircase algorithm
%
if stcs == 1
    %  [Fno,Gno,Hno,Q,k] = obsvf(Fp,G_p,H_pbar);
    [Fno, Gno, Hno, Q, k] = mobsvf(Fp, G_p, H_pbar);
    %  error('prueba')
    nalphao = sum(k);
    models.nalphao = nalphao;
end
Fno = Q * Fp * Q';
Gno = Q * G_p;
Hno = H_pbar * Q';
% nalphao=5;
%
% low-frequency model
%
first = nalpha - nalphao + 1;
T = Fno(first:nalpha, first:nalpha);
H = Gno(first:nalpha, :);
G = J_pbar;
Z = Hno(:, first:nalpha);
% T
% H
% G
% Z
W = [];
C = [];
A = eye(nalphao);
[nC, junk] = size(C);
[junk, nA] = size(A);
ins = [C, A];
i = [nC, 0, 0, nA];

modelsn.Z = Z;
modelsn.T = T;
modelsn.X = X;
modelsn.W = W;
modelsn.G = G;
modelsn.H = H;
modelsn.ins = ins;
modelsn.i = i;
modelsn.stord = stord;
modelsn.conc = conc;
modelsn.trend = trend;
modelsn.slope = slope;
modelsn.seas = seas;
modelsn.irreg = irreg;
modelsn.arp = arp;
modelsn.nalphao = nalphao;
