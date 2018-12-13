function [nrout, nind, tip, Yo] = outlr(y, Y, parm, iout, infm, npr, x0, sp1, sp2, ...
    fmarqdt, ols, aa)
%
% this function performs automatic outlier detection
% three types of outlier are considered:
% AO: defined by a one at t=T and zeros elsewhere
% TC: defined by a one at t=T followed by delta^i at t=T+i, i=1,2,...
%     Usually, delta=.7.
% LS: defined by ones from t=T to the end.
%
% Input arguments:
% y: vector containing the data
% Y: matrix containing regression variables
% parm: astructure containing model information, where
% .s:  seasonality
% .S:  second seasonality
% .p:  AR order
% .ps: order of the AR of order s
% .q:  order of the regular MA
% .qs: order of the MA of order s (1 at most)
% .qS: order of the MA of order S (1 at most)
% .dr: order of regular differencing
% .ds: order of differencing of order s
% .dS: order of differencing of order S
% .pvar:  array containing the indices of variable parameters
% .pfix:  array containing the indices of fixed parameters
% iout:   a structure containing information for outlier detection,
%         where
% .C:     critical value for outlier detection
% .delta: the value for delta in TC outliers
% .mthd: method to compute ARMA parameter estimates (0 Hannan Rissanen,
%                                                    1 Max. Lik.)
% .schr: =0 outliers of type AO and TC are considered (default)
%        =1 outliers of type AO, TC and LS are considered
%   infm     : structure containing function names and optimization
%              options
%   .f  :   a function to evaluate the vector ff of individual functions
%           such that ff'*ff is minimized
%   .tr :   >0 x is passed from marqdt to f but not passed from f to
%           marqdt
%           =0 x is passed from marqdt to f and passed from f to marqdt
%   .tol:   a parameter used for stopping
%   .jac:   =1 evaluation of jacobian and gradient at the solution is
%              performed
%           =0 no evaluation of jacobian and gradient at the solution is
%              performed
% .maxit:   maximum number of iterations
%   .nu0:   initial value of the nu parameter
%   .prt:   =1 printing of results
%           =0 no printing of results
%   .chb:   = 1  compute the beta estimate and its MSE
%             0  do not compute the beta estimate and its MSE
%   .inc:   = 0, the initial states in the filter equations to obtain the
%                filtered variables are equal to zero (not estimated)
%           = 1, the initial states in the filter equations are estimated
% x0: initial parameter vector
% npr: number of forecasts
% Note that AO and LS can be obtained by setting delta=1 and schr=2
% sp1,sp2: the outliers are searched in the time span (sp1,sp2)
% fmarqdt: a parameter for the estimation method
%          = 1 Levenberg-Marquardt method
%          = 0 Lsqnonlin (Matlab)
%
% Output arguments:
% nrout: number of outliers detected
% nind : index numbers of detected outliers
% tip  : array containing the type of the detectec outliers
% Yo: new desing matrix containing the detected outliers in the last
%     columns. That is, Yo=[Y X], where X contains the outlier
%     variables.
%
% Copyright (c) 21 July 2015 by Victor Gomez
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

%

if ~isstruct(parm)
    error('outlr: requires a parameter structure');
end;
if ~isstruct(iout)
    error('outlr: requires an outlier structure');
end;
if ~isstruct(infm)
    error('outlr: requires a structure for optimization');
end;


s = parm.s;
S = parm.S;
dr = parm.dr;
ds = parm.ds;
dS = parm.dS;
p = parm.p;
ps = parm.ps;
q = parm.q;
qs = parm.qs;
qS = parm.qS;
pfix = parm.pfix;
pvar = parm.pvar;
mthd = iout.omet;
C = iout.C;
delta = iout.delta;
schr = iout.schr;

n = length(y);
% maximum number of outliers allowed
nmaxout = floor(.15*n);
% maximum number of eliminated outliers allowed
nmaxelim = 2 * nmaxout;
[mY, nY] = size(Y);
YY = Y;
b = 1;
a = [1, -delta];
nvd = max(n+npr, mY);
% vdelta is a vector containing the delta^i weights
vdelta = poldiv(b, a, nvd);
tol = 1.e-10;
xf = x0(pfix);

% difference series and compute OLS estimator
[yd, beta] = diffest(y, Y, s, S, dr, ds, dS, 1);

%main loop for outlier detection
finish = 0;
nrout = 0;
nind = [];
ntip = [];
eind = [];
etip = [];
search = 1;
elim = 0;
while finish == 0
    out = 0;
    mr = 0;
    % compute initial conditions for parameter vector by Hannan-Rissanen method
    x = inest(yd, beta, s, S, p, ps, q, qs, qS, ols, aa);
    chk = chkroots(x, p, ps, q, qs, qS);
    if chk == 1
        x = invroots(x, p, ps, q, qs, qS);
        if mthd == 0
            if fmarqdt == 0
                %       [x,J,exitflag]=optimize(x,y,Y,s,S,p,ps,dr,ds,dS,q,qs,qS,pfix,pvar);
                [x, J, exitflag] = optim1(x, y, Y, parm, infm);
            else
                [xx, J, exitflag] = marqdt(infm, x(pvar), y, Y, parm, infm, xf);
                x(pvar) = xx;
            end
        end
    elseif chk == -1
        nrout = 0;
        nind = [];
        tip = [];
        Yo = YY;
        return
    end
    x(pfix) = xf;
    if mthd == 1
        if fmarqdt == 0
            %       [x,J,exitflag]=optimize(x,y,Y,s,S,p,ps,dr,ds,dS,q,qs,qS,pfix,pvar);
            [x, J, exitflag] = optim1(x, y, Y, parm, infm);
        else
            [xx, J, exitflag] = marqdt(infm, x(pvar), y, Y, parm, infm, xf);
            x(pvar) = xx;
        end
    end
    [nd, md] = size(yd);
    % compute L matrix using the Durbin-Levinson algorithm
    %  L=lm1(x,nd,s,S,p,ps,q,qs,qS);  %L is sparse matrix (lower triangular)
    % compute L matrix using the fast Kalman filter algorithm
    L = lm1KF(x, nd, s, S, p, ps, q, qs, qS); %L is sparse matrix (lower triangular)
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%
    % filter the differenced series.
    ydf = L * yd;
    if md > 1
        [beta, tv] = btval(x, ydf);
        res = ydf(:, 1) - ydf(:, 2:md) * beta;
    else
        res = ydf;
    end
    if search == 1 && nrout < nmaxout
        % compute the MAD estimator of the standard deviation of the residuals
        mad = 1.483 * median(abs(res-median(res)));
        tmax = -1;
        ind = 0;
        bmax = 0;
        tipomax = 0;
        % detect outliers one by one
        for i = sp1:sp2
            % form outlier variables
            v1 = zeros(n, 1);
            v1(i) = 1;
            v2 = v1;
            v2(i:n) = vdelta(1:n-i+1);
            if schr == 0
                % difference outlier variables
                [v, junk] = diffest(v1, v2, s, S, dr, ds, dS, 0);
            else
                v3 = v1;
                v3(i:n) = ones(n-i+1, 1);
                [v, junk] = diffest(v1, [v2, v3], s, S, dr, ds, dS, 0);
            end
            % filter the outlier variables and compute regression t statistics
            v = L * sparse(v);
            v1r = v(:, 1)' * res;
            v1v1 = v(:, 1)' * v(:, 1);
            v2r = v(:, 2)' * res;
            v2v2 = v(:, 2)' * v(:, 2);
            t1 = v1r / (sqrt(v1v1) * mad);
            t1 = abs(t1);
            t2 = v2r / (sqrt(v2v2) * mad);
            t2 = abs(t2);
            if schr == 1
                v3r = v(:, 3)' * res;
                v3v3 = v(:, 3)' * v(:, 3);
                % v(:,3) can be zero if t=1 and a regular difference is applied
                if v3v3 > tol
                    t3 = v3r / (sqrt(v3v3) * mad);
                    t3 = abs(t3);
                else
                    t3 = -1;
                end
            else
                t3 = -1;
            end
            % compute maximum t for the observation
            if t1 > t2 && t1 > t3
                t = t1;
                b = v1r / v1v1;
                tipo = 1;
            elseif t2 > t1 && t2 > t3
                t = t2;
                b = v2r / v2v2;
                tipo = 2;
            elseif schr == 1
                t = t3;
                b = v3r / v3v3;
                tipo = 3;
            end
            % compare the maximum t with the global maximum t so far computed
            if t > tmax
                tmax = t;
                ind = i;
                bmax = b;
                tipomax = tipo;
            end
        end
        if tmax > C
            if tipomax == 1
                v1 = zeros(nvd, 1);
                v1(ind) = 1;
                Y = [Y, v1];
                nrout = nrout + 1;
                nind = [nind; ind];
                ntip = [ntip; 1];
            elseif tipomax == 2
                v2 = zeros(nvd, 1);
                v2(ind:nvd) = vdelta(1:nvd-ind+1);
                Y = [Y, v2];
                nrout = nrout + 1;
                nind = [nind; ind];
                ntip = [ntip; 2];
            elseif tipomax == 3
                v3 = zeros(nvd, 1);
                v3(ind:nvd) = ones(nvd-ind+1, 1);
                Y = [Y, v3];
                nrout = nrout + 1;
                nind = [nind; ind];
                ntip = [ntip; 3];
            end
            beta = [beta; bmax];
            [yd, beta] = diffest(y, Y, s, S, dr, ds, dS, 1);
            out = 1;
        end
    end
    %  nrout
    %  nind
    %  ntip
    % multiple regression (stepwise elimination)
    if out == 0
        md1 = md - 1;
        imin = 0;
        if md1 > nY
            tmin = 9999;
            nout = md1 - nY;
            for i = 1:nout
                tvm = abs(tv(nY+i));
                if tvm < C && tvm < tmin
                    tmin = tvm;
                    imin = i;
                end
            end
            %    imin
            %    tv
            if imin > 0
                % Eliminate all non significant outliers one by one.
                mr = 1;
                elim = elim + 1;
                ncoin = coincid(nind, eind);
                %     nind
                %     ntip
                %     eind
                %     etip
                if ncoin ~= 0
                    % Do not search for more outliers if the present outlier to be eliminated
                    % coincides with some of the previously eliminated outliers.
                    search = 0;
                end
                % store information on the eliminated outlier
                [eind, etip] = sinfelo(eind, etip, nind, ntip, imin);
                % do not search for more outliers if the number of eliminated outliers
                % so far is greater than nmaxout or greater than nmaxelim
                if length(eind) >= nmaxout || elim >= nmaxelim
                    search = 0;
                end
                % delete the column corresponding to the outlier with minimum t<C
                Y(:, nY+imin) = [];
                nind(imin) = [];
                ntip(imin) = [];
                beta(imin) = [];
                [nn, mm] = size(Y);
                mnm = min(nn, mm);
                if mnm == 0, Y = [];
                end
                [yd, beta] = diffest(y, Y, s, S, dr, ds, dS, 1);
                nrout = nrout - 1;
            end
        end
        if mr == 0
            finish = 1;
        end
    end
end

% search is finished. Now complete Output
Yo = Y;
tip = [];
for i = 1:nrout
    if ntip(i) == 1
        tip = [tip; 'AO'];
    elseif ntip(i) == 2
        tip = [tip; 'TC'];
    else
        tip = [tip; 'LS'];
    end
end
