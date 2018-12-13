function [F, xv, xf, e, f, g, M, A, P] = lkevarmapqPQ(xv, y, Y, xf, str, chb)
% PURPOSE: given a structure containing information about a VARMA model, it
% evaluates the likelihood of that model after putting it into state space
% form.
%---------------------------------------------------
% USAGE: [F,xv,xf,e,f,g,M,A,P] = lkevarmapqPQ(xv,y,Y,xf,str,chb)
% where:
%          xv       = a vector containing the parameters to be estimated
%          y        = an (n x neqs) matrix containing the data
%          Y        = an (n x (neqs x nbeta)) matrix containing regression
%                     matrix
%          xf       = a vector containing the fixed parameters
%          str      = a structure containing the initial model information
%          chb      =1, compute regression estimate and covariance matrix
%                   =0, do not compute
%---------------------------------------------------
%---------------------------------------------------
% RETURNS: F    = a vector containing the individual functions at the
%                 solution
%         xv    = a vector containing the estimated parameters
%         xf    = a vector containing the fixed parameters
%          e    = a vector containing the standardized residuals
%          f    = a scalar containing the determinantal term
%          g    = a vector containing the regression estimates
%          M    = a matrix containing the mse of g
%          A    = the estimated augmented state vector
%                 at the end of filtering
%          P    = the mse of A at the end of filtering
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

if ~isfield(str, 'profile')
    profile = 0;
else
    profile = str.profile;
end


%pass parameters to VARMA model
[phi, th, Phi, Th, L] = pr2varmapqPQ(xv, xf, str);

%transform multiplicative to non-multiplicative model
[phirs, thrs, H, F, G, J, ferror] = varmapqPQ2ssm(phi, th, Phi, Th, L, str);

%initial state covariance matrix
nxi = length(str.xi);
nxv = sum(str.xi);
nxf = nxi - nxv;
xxv = zeros(1, nxv);
xxf = str.xf;
xxi = str.xi;
%check stationarity
[U, T] = schur(F', 'complex');
T = T';
nonstat = size(find(abs(diag(T)) >= .995), 1);
%check stationarity
if nonstat > 0
    %  fprintf(1,'model nonstationary, nonstat = %2d\n',nonstat);
    phi = enfstabpol(phi);
    Phi = enfstabpol(Phi);
    [phirs, thrs, H, F, G, J, ferror] = varmapqPQ2ssm(phi, th, Phi, Th, L, str);
    [str, ferror] = suvarmapqPQ(phi, th, Phi, Th, str.Sigma, str.freq);
    cont = 0; %contf=0;
    for i = 1:nxi
        if xxi(i) == 1
            cont = cont + 1;
            xxv(cont) = str.xv(i);
            %   else
            %    contf=contf+1;
            %    xxf(contf)=str.xv(i);
        end
    end
    str.xv = xxv;
    str.xf = xxf;
    str.xi = xxi;
    xv = str.xv;
    xf = str.xf;
    if (profile == 0)
        %  [Sigma,ferror] = mlyapunov(F,G*G');
        Sigma = dlyapsq(F, G);
        Sigma = Sigma' * Sigma;
        %  Sigma = dlyap(F,G*G');      %matlab Lyapunov
    end
end
%check invertibility
FKH = F - (G / L) * H;
[U, T] = schur(FKH', 'complex');
T = T';
noninv = size(find(abs(diag(T)) >= .995), 1);
if noninv > 0
    %  fprintf(1,'model noninvertible, noninv = %2d\n',noninv);
    th = enfstabpol(th);
    Th = enfstabpol(Th);
    [phirs, thrs, H, F, G, J, ferror] = varmapqPQ2ssm(phi, th, Phi, Th, L, str);
    [str, ferror] = suvarmapqPQ(phi, th, Phi, Th, str.Sigma, str.freq);
    cont = 0; %contf=0;
    for i = 1:nxi
        if xxi(i) == 1
            cont = cont + 1;
            xxv(cont) = str.xv(i);
            %   else
            %    contf=contf+1;
            %    xxf(contf)=str.xv(i);
        end
    end
    str.xv = xxv;
    str.xf = xxf;
    str.xi = xxi;
    xv = str.xv;
    xf = str.xf;
    if (profile == 0)
        %  [Sigma,ferror] = mlyapunov(F,G*G');
        Sigma = dlyapsq(F, G);
        Sigma = Sigma' * Sigma;
        %  Sigma = dlyap(F,G*G');   %matlab Lyapunov
    end
end
if (nonstat + noninv == 0) && (profile == 0)
    %  [Sigma,ferror] = mlyapunov(F,G*G');
    Sigma = dlyapsq(F, G);
    Sigma = Sigma' * Sigma;
    % Sigma = dlyap(F,G*G');      %matlab Lyapunov
end
[nalpha, mf] = size(F);

%set up regression matrices
X = Y;
W = [];


%set up initial state
if profile == 1
    [ny, my] = size(y);
    C = zeros(nalpha);
    [nx, mx] = size(X);
    if nx == my
        X = [X, zeros(my, nalpha)];
    else
        X = [X, repmat(zeros(my, nalpha), ny, 1)];
    end
    [junk, nX] = size(X);
    if nX > nalpha
        W0 = [zeros(nalpha, nX-nalpha), eye(nalpha)];
    else
        W0 = eye(nalpha);
    end
    [junk, nW] = size(W0);
    ins = [C, W0];
    i = [nalpha, nW, 0, 0];
else
    ins = Sigma;
    i = [nalpha, 0, 0, 0];
end

%        chb= 1 compute g and M
%             0 do not compute g and M
%set up system matrices
T = F;
Z = H;
GG = J;
HH = G;
%likelihood evaluation
if profile == 0
    %  [e,f,g,M,A,P]=scakfle2(y,X,Z,GG,W,T,HH,ins,i,chb);
    [e, f, g, M, A, P] = scakflesqrt(y, X, Z, GG, W, T, HH, ins, i, chb);
else
    [e, f, g, M, A, P] = scakflepc(y, X, Z, GG, W, T, HH, ins, i, chb);
end
F = e * f;
