function models = modstr_viviusa(y, xx, pfix, pvar, xf, stordt, conc, n, r)
%**************************************************************************
% Auxiliary function called in viviusa_d.m to set up the state space model
%
%  INPUTS:
%        xx : array with parameters to be estimated
%      pfix : array with fixed parameter indices
%      pvar : array with variable parameter indices
%        xf : array with fixed parameters
%    stordt : array index for the standard deviations in the model
%      conc : index for the standard deviation to be concentrated out
%         n : number of variables in the data y
%         r : number of parameters in matrix K
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

%number of parameters in K matrix
nk = 0;
for i = 1:r
    nk = nk + n - i;
end

%parameters for K matrix
xk = x(1:nk);
%parameters for standard deviations
xs = zeros(1, length(stordt));
xs(stord) = x(nk+1:end);
xs(conc) = 1.; %concentrated parameter


% Multivariate monthly structural model
npr = n + r;
np1 = n + 1;
tn = 2 * n;
tnpr = tn + r;
%matrix Tm
Tm = eye(npr);
K = zeros(n, r);
K(1:r, 1:r) = eye(r);
nk = 0;
for i = 1:r
    nkn = nk + n - i;
    K(i+1:n, i) = xk(nk+1:nkn)';
    nk = nkn;
end
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
[insm, im] = incossm(Tm, Hm, npr);


%matrix Hm
Hmu = zeros(npr, npr);
Hmu(1:n, 1:n) = diag(xs(1:n)');
Hmu(np1:npr, np1:npr) = diag(xs(tn+1:tnpr)');
%matrix Gm
Gmu = zeros(n, npr);

models.xv = xx;
models.xf = xf;
models.pfix = pfix;
models.pvar = pvar;
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
models.stordt = stordt;
models.conc = conc;
models.n = n;
models.r = r;

[ny, p] = size(y);
% %interventions for variable VU according to the results given by TRAMO. The
% %interventions are: TC 367, LS 360, AO 37.
% X1=zeros(ny*p,3); fac=.7;
% for j=1:ny
%  ip=(j-1)*p+1;
% %intervention TC 367
%  if (j >= 367)
%   X1(ip,1)=fac;
%   fac=min(1.d-10,fac*fac);
%  end
% %intervention LS 360
%  if (j >= 360)
%   X1(ip,2)=1.;
%  end
% %intervention AO 37
%  if (j == 37)
%   X1(ip,3)=1.;
%  end
% end

X1 = [];

%interventions for variable PVN according to the results given by TRAMO.
%The interventions are: AO 370, LS 89, LS 344, LS 127, LS 289.
X2 = zeros(ny*p, 5);
for j = 1:ny
    ip = (j - 1) * p + 4;
    %intervention AO 370
    if (j == 370)
        X2(ip, 1) = 1.;
    end
    %intervention LS 89
    if (j >= 89)
        X2(ip, 2) = 1.;
    end
    %intervention LS 344
    if (j >= 344)
        X2(ip, 3) = 1.;
    end
    %intervention LS 127
    if (j >= 127)
        X2(ip, 4) = 1.;
    end
    %intervention LS 289
    if (j >= 289)
        X2(ip, 5) = 1.;
    end
end

X = [X1, X2];

models.X = X;
