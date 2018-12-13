function [result, ferror] = varmapqPQestim(y, str, Y)
%
% This function estimates a VARMA(p,q)(P,Q)_s model using the exact maximum
% likelihood method.
%
%
% Inputs: y: matrix containing the input series
%         Y: an (n*neqs x nbeta) matrix containing the regression matrix
%       str: a structure containing the initial model information given by
%            function suvarmapqPQ.m. If an additional field, .nmiss, is
%            present, it is used to control output in the case of missing
%            observations.
%
%  Output: .xvf : estimated parameters
%           .xf : vector of fixed parameters
%      .sigma2c : concentrated parameter estimate
%       .Sigmar : estimated exact covariance matrix of residuals
%           .tv : t-values of the estimated varma parameters
%    .residexct : matrix containing recursive residuals, only if Y is empty
%            .e : vector of standardized residuals at the end of estimation
%                 (Q'_2*y)
%           .ff : vector of nonlinear functions whose sum of squares is
%                 minimized at the end of estimation
%           .h  : vector of estimated regression estimates
%           .H  : matrix of mse of h
%           .A  : estimated state vector, x_{t|t-1}, obtained with the
%                 Kalman filter at the end of the sample
%           .P  : Mse of A
%           .tvr: vector of t-values for h
%       .ferror : flag for errors
%
% Copyright (c) January 2010 by Victor Gomez
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

result = [];
ferror = 0;

if ~isfield(str, 'profile')
    profile = 0;
else
    profile = str.profile;
end
if ~isfield(str, 'prt')
    prt = 2;
else
    prt = str.prt;
end


[ny, my] = size(y); %each row is an observation, each column is a variable
[nY, mY] = size(Y);
nbeta = mY;

%check stationarity and invertibility.
H = str.Z;
F = str.T;
G = str.H;
J = str.G;
nxi = length(str.xi);
nxv = sum(str.xi);
nxf = nxi - nxv;
xv = zeros(1, nxv);
xf = zeros(1, nxf);
xi = str.xi;
%check stationarity
[U, T] = schur(F', 'complex');
T = T';
nonstat = size(find(abs(diag(T)) >= .995), 1);
if nonstat > 0
    if (prt > 0)
        fprintf(1, 'initial model nonstationary in varmapqPQestim\n');
    end
    %  ferror=1; return
    phi = enfstabpol(str.phi);
    Phi = enfstabpol(str.Phi);
    th = str.th;
    Th = str.Th;
    Sigma = str.Sigma;
    freq = str.freq;
    xf = str.xf;
    [str, ferror] = suvarmapqPQ(phi, th, Phi, Th, Sigma, freq);
    cont = 0; %contf=0;
    for i = 1:nxi
        if xi(i) == 1
            cont = cont + 1;
            xv(cont) = str.xv(i);
            %   else
            %    contf=contf+1;
            %    xf(contf)=str.xv(i);
        end
    end
    str.xv = xv;
    str.xf = xf;
    str.xi = xi;
end
%check invertibility
FKH = F - (G / J) * H;
[U, T] = schur(FKH', 'complex');
T = T';
noninv = size(find(abs(diag(T)) >= .995), 1);
if noninv > 0
    if (prt > 0)
        fprintf(1, 'initial model noninvertible in varmapqPQestim\n');
    end
    %  ferror=1; return
    th = enfstabpol(str.th);
    Th = enfstabpol(str.Th);
    phi = str.phi;
    Phi = str.Phi;
    Sigma = str.Sigma;
    freq = str.freq;
    xf = str.xf;
    [str, ferror] = suvarmapqPQ(phi, th, Phi, Th, Sigma, freq);
    cont = 0; %contf=0;
    for i = 1:nxi
        if xi(i) == 1
            cont = cont + 1;
            xv(cont) = str.xv(i);
            %   else
            %    contf=contf+1;
            %    xf(contf)=str.xv(i);
        end
    end
    str.xv = xv;
    str.xf = xf;
    str.xi = xi;
end

%initial parameters
xv = str.xv;
xf = str.xf;
chb = 0;

str.profile = profile;


f = 'lkevarmapqPQ';
tr = 0;
mvx = 1; %parameters for Levenberg-Marquardt method:
tolf = 1e-4; %f  :   a function to evaluate the vector ff of individual
maxit = 100;
nu0 = .01; %       functions such that ff'*ff is minimized
jac = 1; %prt=2;             %tr :   >0 x is passed from marqdt to f but not passed from f
clear infm
infm = minfm(f, tr, mvx, tolf, maxit, nu0, jac, prt, [], []);
%Levenberg-Marquardt method
[xvf, J, ff, g, iter, conf] = marqdt(infm, xv, y, Y, xf, str, chb);
%this part added 21-2-2011
if (abs(sum(sum(J))) <= 1.d-8)
    if (prt > 0)
        disp('Jacobian computation failed in varmapqPQestim')
    end
    ferror = 1;
end
%end of addition

if ~isempty(Y) || (profile == 1)
    chb = 1;
end
[F, xvf, xf, e, f, h, H, A, P] = lkevarmapqPQ(xvf, y, Y, xf, str, chb);
sigma2c = e' * e / double(ny*my-nbeta); %concentrated parameter estimate

Lh = str.Lh;
nlhm1 = length(Lh) - 1;
L = zeros(my);
%inser in xx all parameters
nparm = str.nparm;
xx = zeros(1, nparm);
xi = str.xi;
j = 0;
jj = 0;
for i = 1:nparm
    if xi(i) == 1
        j = j + 1;
        xx(i) = xvf(j);
    else
        jj = jj + 1;
        xx(i) = xf(jj);
    end
end
betam = [1., xx(end-nlhm1+1:end)];

l = 0;
for i = 1:my
    cont = my - i + 1;
    ind = l + 1:l + cont;
    L(i:end, i) = betam(ind);
    l = l + cont;
end
Sigmar = L * L' * sigma2c; %estimated exact covariance matrix of residuals
residexct = []; %estimated exact residuals
if isfield(str, 'nmiss')
    nmiss = str.nmiss;
else
    nmiss = 0;
end
if isempty(Y) && (profile == 0) && (nmiss == 0)
    for i = 1:ny
        ind = (i - 1) * my + 1:i * my;
        residexct(i, 1:my) = e(ind);
    end
end

if ferror == 0
    %t-values
    V = -J; % negative derivatives for the regression
    res2t = ff + V * xvf'; % regression
    [beta, tv] = btval([], [res2t, V]); % tv contains the t-values
else
    tv = [];
end

result.xvf = xvf;
result.xf = xf;
result.sigma2c = sigma2c;
result.Sigmar = Sigmar;
result.tv = tv;
result.residexct = residexct;
result.e = e;
result.ff = ff;
result.h = h;
result.H = H;
result.A = A;
result.P = P;
if ~isempty(Y) || (profile == 1)
    %t-values of regression parameters
    H = H * sigma2c;
    tvr = h ./ sqrt(diag(H));
    result.tvr = tvr;
else
    result.tvr = [];
end
result.ferror = ferror;
