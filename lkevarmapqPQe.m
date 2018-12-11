function [F, xv, xf, e, f, g, M, A, P] = lkevarmapqPQe(xv, y, Y, xf, str, chb, constant)
% PURPOSE: given a structure containing information about a cointegrated
% VARMA model, it evaluates the likelihood of that model after putting it
% into state space form.
%---------------------------------------------------
% USAGE: [F,xv,xf,e,f,g,M] = lkevarmapqPQe(xv,y,Y,xf,str,chb,constant)
% where:   y        = an (n x neqs) matrix containing the data
%          Y        = an (n x (neqs x nbeta)) matrix containing regression
%                     matrix
%          xv       = a vector containing the parameters to be estimated
%          xf       = a vector containing the fixed parameters
%          xi       = an index vector, if xi(i)=1, the i-th parameter is to
%                     be estimated, =0, not
%          str      = a structure containing the initial model information
%          chb      =1, compute regression estimate and covariance matrix
%                   =0, do not compute
%     constant      =1 a constant should be included in the model for the
%                      differenced series
%                    0 no constant in the model for the differenced series
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

%obtain differenced series and parameters for the VARMA model followed by
%the differenced series.
[yd, xvv, xff, DA, Dr, Ds, ferror] = pr2varmapqPQd(y, xv, xf, str);
[ny, my] = size(y);
X = [];
if ~isempty(Y)
    Xd = Y;
    [nxd, mxd] = size(Xd);
    if nxd == my
        %we make Y always time varying
        Xd = zeros(ny*my, mxd);
        for i = 1:ny
            ip = (i - 1) * my + 1:i * my;
            Xd(ip, :) = Y;
        end
    end
    if isfield(str, 'nr')
        nr = str.nr;
    else
        nr = 0;
    end
    if isfield(str, 'ns')
        ns = str.ns;
    else
        ns = 0;
    end
    if (ns > 0)
        Xdcold = [];
        for i = 1:mxd
            Xdcol = Xd(:, i);
            Xdcolt = reshape(Xdcol, ny, my);
            [XX, ferror] = dtimesy(Ds, Xdcolt);
            Xdcold = [Xdcold, vec(XX)];
        end
        Xd = Xdcold;
    end
    if (nr > 0)
        Xdcold = [];
        for i = 1:mxd
            Xdcol = Xd(:, i);
            Xdcolt = reshape(Xdcol, ny, my);
            [XX, ferror] = dtimesy(Dr, Xdcolt);
            Xdcold = [Xdcold, vec(XX)];
        end
        Xd = Xdcold;
    end
    if constant == 1
        Cst = eye(my);
        [nxd, mxd] = size(Xd);
        [nyd, myd] = size(yd);
        X = zeros(nxd, my+mxd);
        for i = 1:nyd
            ip = (i - 1) * my + 1:i * my;
            X(ip, :) = [Cst, Xd(ip, :)];
        end
    else
        X = Xd;
    end
else
    if constant == 1
        X = eye(my);
    end
end


%pass parameters to ecf model
[Lambda, alpha, betap, th, Th, L, ferror] = pr2ecf(xvv, xff, DA, str);
%transform error correction model to model in differences
[phi, D, DA, ferror] = mecf2mid(Lambda, alpha, betap);
Phi(:, :, 1) = eye(my);

%transform multiplicative to non-multiplicative model
[phirs, thrs, H, F, G, J, ferror] = varmapqPQ2ssm(phi, th, Phi, Th, L, str);


%initial state covariance matrix
nxi = length(str.xi);
nxv = sum(str.xi);
nxf = nxi - nxv;
xxv = zeros(1, nxv);
xxf = str.xf;
xxi = str.xi;
% nr=str.nr; ns=str.ns;
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
    [Pi, Lambda, alpha, betap, ferror] = mid2mecf(phi, D, DA);
    xvd = str.xvd;
    xfd = str.xfd;
    xid = str.xid;
    [str, ferror] = suvarmapqPQe(Lambda, alpha, betap, th, Th, str.Sigma, str.freq);
    %  [str,ferror] = suvarmapqPQ(phi,th,Phi,Th,str.Sigma,str.freq);
    %  [str,ferror] = aurivarmapqPQ(str,nr,ns,DA);
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
    str.xvd = xvd;
    str.xfd = xfd;
    str.xid = xid;
    xv = [xvd, xv];
    xf = [xfd, xf];
    %  [Sigma,ferror] = mlyapunov(F,G*G');
    Sigma = dlyapsq(F, G);
    Sigma = Sigma' * Sigma;
    %  Sigma = dlyap(F,G*G');      %matlab Lyapunov
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
    [Pi, Lambda, alpha, betap, ferror] = mid2mecf(phi, D, DA);
    xvd = str.xvd;
    xfd = str.xfd;
    xid = str.xid;
    [str, ferror] = suvarmapqPQe(Lambda, alpha, betap, th, Th, str.Sigma, str.freq);
    %  [str,ferror] = suvarmapqPQ(phi,th,Phi,Th,str.Sigma,str.freq);
    %  [str,ferror] = aurivarmapqPQ(str,nr,ns,DA);
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
    str.xvd = xvd;
    str.xfd = xfd;
    str.xid = xid;
    xv = [xvd, xv];
    xf = [xfd, xf];
    %  [Sigma,ferror] = mlyapunov(F,G*G');
    Sigma = dlyapsq(F, G);
    Sigma = Sigma' * Sigma;
    %  Sigma = dlyap(F,G*G');   %matlab Lyapunov
end
if (nonstat + noninv == 0)
    %  [Sigma,ferror] = mlyapunov(F,G*G');
    Sigma = dlyapsq(F, G);
    Sigma = Sigma' * Sigma;
    % Sigma = dlyap(F,G*G');      %matlab Lyapunov
end
[nalpha, mf] = size(F);


%set up initial state
ins = Sigma;
i = [nalpha, 0, 0, 0];
%set up regression matrices
% X=Y;
W = [];
%        chb= 1 compute g and M
%             0 do not compute g and M
%set up system matrices
T = F;
Z = H;
GG = J;
HH = G;
%likelihood evaluation
[e, f, g, M, A, P] = scakflesqrt(yd, X, Z, GG, W, T, HH, ins, i, chb);
F = e * f;
