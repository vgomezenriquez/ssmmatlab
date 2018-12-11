function models = modstr_mulcycuswcv(xx, pfix, pvar, xf, stordt, conc, n, r, l)
%**************************************************************************
% Auxiliary function called in mulcycusfunwcv to set up the state space
% model
%
%  INPUTS:
%        xx : array with parameters to be estimated
%      pfix : array with fixed parameter indices
%      pvar : array with variable parameter indices
%        xf : array with fixed parameters
%    stordt : array with indices for standard deviations
%      conc : index for the standard deviation to be concentrated out
%         n : number of variables in the data matrix y
%         r : number of slope factors
%         l : numbef of quarterly series
%
%    OUTPUTS:
%    models : structure with multivariate model information
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

npar = length(pfix) + length(pvar);
x = zeros(1, npar);
x(pfix) = xf;
x(pvar) = xx;

stord = setdiff(stordt, conc); %standard deviations to be estimated

%parameters for K matrix
xk = x(1:n-1);
%parameters for standard deviations
xs = zeros(1, length(stordt));
xs(stord) = x(n:end);
xs(conc) = 1.; %concentrated parameter


% Multivariate monthly structural model
npr = n + r;
np1 = n + 1;
tn = 2 * n;
tnpr = tn + r;
%matrix Tm
Tm = eye(npr);
K = zeros(n, r);
K(:, 1) = [1.; xk'];
Tm(1:n, np1:npr) = K;
%matrix Zm
Zm = [eye(n), zeros(n, r)];
%matrix Hm
Hm = zeros(npr, tnpr);
Hm(1:n, 1:n) = diag(xs(1:n)');
Hm(np1:npr, np1:npr) = diag(xs(tn+1:tnpr)');
%matrix Gm
Gm = zeros(n, tnpr);
Gm(:, npr+1:tnpr) = diag(xs(np1:tn)');
%initial conditions
insm = eye(npr);
im = [0, 0, 0, npr];

% Multivariate monthly structural model for the trend component, mu_t
%matrix Hm
Hmu = zeros(npr, npr);
Hmu(1:n, 1:n) = diag(xs(1:n)');
Hmu(np1:npr, np1:npr) = diag(xs(tn+1:tnpr)');
%matrix Gm
Gmu = zeros(n, npr);
models.Tm = Tm;
models.Hm = Hmu;
models.Zm = Zm;
models.Gm = Gmu;
models.W = [];
models.X = [];

models.T = Tm;
models.H = Hm;
models.Z = Zm;
models.G = Gm;
models.W = [];
models.X = [];
models.i = im;
models.ins = insm;
models.stord = stord;
models.conc = conc;
