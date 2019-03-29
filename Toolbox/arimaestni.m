function outa = arimaestni(dbname, ser, fidr, ii)
%
% function to identify, estimate and forecast an ARIMA model for one
% series. The series may have up to two seasonalities. The ARIMA model is
% of the form:
%
% phi(B)*phi_s(B^s)*phi_S(B^S)*(delta*delta_s*delta_S*y_t -mu) =
% th(B)*th_s(B^s)*th_S(B^S)*a_t
%
% In the subdirectory spec, there is a specification file where all the
% options for the ARIMA model are defined and returned in the structue ser.
% The name of this file is given in function arimaestos and passed to this
% function.
% These options include, log transformation criteria, automatic
% identification of ARMA model and differencing operators, automatic
% specification of trading day, Easter effect and leap year (for quarterly
% and monthly series only), outlier search and forecasting, among other
% things.
% No automatic model identification is performed for the second seasonality
% (S). This part must be entered by the user. Automatic model
% identification is performed for the regular and the first seasonal part.
% Output is written in an external file in the subdirectory results.
%
%     INPUTS:
%     dbname : name of the series
%     ser    : a structure, containing the instructions for this function
%     fidr   : an integer, corresponding to the output file
%     ii     : an integer, corresponding to the series currently handled.
%
%  OUTPUTS:
%      outa  : a structure containing model information for the input
%              with fields:
%       title: a string with the name of series
%      nziyip: a 1 x 3 array with number of obs., initial year, initial per.
%        freq: number of seasons
%        orig: original series
%       model: structre with ARIMA model information. In the case of a
%              transfer function, it is the ARIMA model corresponding to
%              the finite linear approximation to the input filters. It has
%              the following fields:
%          lam: = 0, logs are taken, = 1, no logs
%         mean: = 1, a mean is added, = 0, no mean
%            p: degree of regular AR polynomial
%            d: degree of regular differencing
%            q: degree of regular MA polynomial
%           ps: degree of seasonal AR polynomial
%           ds: degree of seasonal differencing
%           qs: degree of seasonal MA polynomial
%         nreg: the number of regression variables
%       result: a structure containing estimation results
%          phi: an array containing the regular AR polynomial coefficients
%         phis: an array containing the seasonal AR polynomial coefficients
%           th: an array containing the regular MA polynomial coefficients
%          ths: an array containing the seasonal MA polynomial coefficients
%        nrout: number of outliers
%            C: critical value for outlier detection
%         nind: observation numbers of the outliers
%          tip: string containing the outlier types
%       matsis: a structure containing the state space form of the model
%       resinf: a structure containing information about the residuals
%           hb: array containing the regression estimates
%           Mb: matrix containing the covariance matrix of the regression
%               estimates
%            Y: matrix containing the total regression effects
%          seb: array containing the regression standard errors
%           tb: array containing the t-values of the regression estimates
%          Yrg: array containing the regression variables that are not
%               outliers
%        Youtg: array containing the outlier variables
%           se: array containing the standard errors of the estimates
%           tt: array containing the t-values of the estimates
%          npr: number of forecasts
%          pry: array containing the forecasts (transformed scale)
%         spry: array containing the standard errors of the forecasts
%         opry: same as pry but in the original scale
%        ospry: same os spry but in the original scale
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

% profile on       %this is to optimize the code. It gives some information
%                  %about the performance of all of the functions.

if nargin < 3
    fidr = 1;
    ii = 0;
elseif nargin < 4
    ii = 0;
end

tic;
clear x

outa.title = dbname;
%
%check for missing values
%
if (~isfield(ser, 'yor'))
    error('field yor must be present in structure ser')
else
    yorm = ser.yor;
end

[yor, Xm, nmiss, idxn] = chmarima(yorm);
outa.yor = yor;
if (nmiss > 0)
    outa.model.nmiss = nmiss;
    outa.model.idxn = idxn;
    rnamesrgm = ['miss', num2str(1)];
    for i = 2:nmiss
        rnamesrgm = char(rnamesrgm, ['miss', num2str(i)]);
    end
    if isfield(ser, 'Y') %Y contains the regression variables
        Y = ser.Y;
        [nY, mY] = size(Y);
        nXm = size(Xm, 1);
        if nY > nXm
            ser.Y = [[Xm; zeros(nY-nXm, nmiss)], Y];
        else
            ser.Y = [Xm, Y];
        end
        if isfield(ser, 'rnamesrg') %names for regression variables
            ser.rnamesrg = char(rnamesrgm, ser.rnamesrg);
        else
            for i = 1:mY
                rnamesrgm = char(rnamesrgm, ['reg', num2str(i)]);
            end
            ser.rnamesrg = rnamesrgm;
        end
    else
        if isfield(ser, 'npr')
            ser.Y = [Xm; zeros(ser.npr, nmiss)];
        else
            ser.Y = Xm;
        end
        ser.rnamesrg = rnamesrgm;
    end
end

if (~isfield(ser, 'bg_year'))
    error('field bg_year must be present in structure ser')
else
    bg_year = ser.bg_year;
end
if (~isfield(ser, 'bg_per'))
    error('field bg_per must be present in structure ser')
else
    bg_per = ser.bg_per;
end
if (~isfield(ser, 'freq'))
    error('field freq must be present in structure ser')
else
    freq = ser.freq;
end

outa.nziyip = [length(yor), bg_year, bg_per];
outa.freq = freq;
outa.orig = yor;


%default values for the program

arimadefval

%start of program

%files for output
if (pr == 1)
    pathc = pwd; %current directory
    if exist([pathc, filesep, 'results'], 'dir') ~= 7
        mkdir(pathc, 'results'); %create directory results
    end
    sname = ['results', filesep, dbname, '.txt'];
    fid = fopen(sname, 'w');
else
    fid = 1;
end


%checking for inconsistencies
if (npr > 0) && (nY < ny + npr) && (nY > 0)
    if (ny ~= nY), disp('regression variables do not have the same length as series'), end
    meanY = mean(Y);
    Y = [Y; zeros(npr, mY)];
    [nY, mY] = size(Y);
    nreg = mY;
    for i = 1:npr, Y(ny+i, :) = meanY;
    end
    if pr == 1
        fprintf(fid, '%s\n', 'Regression matrix is extended with mean values for prediction');
    end
elseif (npr == 0)
    if (ny < nY)
        Y = Y(1:ny, :);
        nY = ny;
    elseif (ny > nY) && ~isempty(Y)
        disp('regression variables have shorter length than series');
        return
    end
end

if ((s ~= 12) && (s ~= 4))
    trad = 0;
    easte = 0;
    leapy = 0;
end


datei = cal(bg_year, bg_per, freq); %initial date and frequency.
if lam == 0, y = log(yor);
else y = yor;
end %logs of the series
flagc0 = 0; %flag for automatic model identification
%default values of C. Also the values for automatic model identification
%if ny <= 50, C0=2.8; elseif (ny> 50 & ny <= 75), C0=3.2; elseif (ny> 75 & ny <= 100),...
%C0=3.5; elseif (ny > 100 & ny <= 200), C0=3.8; elseif (ny >200 & ny<=300), C0=3.9;...
%elseif (ny >300 & ny<=400), C0=4.; elseif (ny >400 & ny<=500), C0=4.1; else, C0=4.2; end
if (C < 0) %critical level for outlier detection
    C = 2.2 + log(log(ny));
end
nreg0 = nreg; %number of initial regression variables
dur = 0; %initial value of dur

if s > 1, maxPQ = 1;
else maxPQ = 0;
end

%create structures
clear parm
parm = mparm(s, S, dr, ds, dS, p, ps, q, qs, qS, ny, nreg, pfix, pvar, lam, flagm, trad, ...
    leapy, easte, dur, ninput, nmiss);

clear iout
nrout = 0;
nind = [];
tip = [];
Yo = [];
ornames = [];
iout = miout(omet, C, delta, schr, nrout, nind, tip, Yo, ornames);
if ((lam < 0)) || (autmid == 1)
    flagc0 = 1; %C0 is used in the log test and automtic model identification
end
if flagc0 == 1, iout.C = C0;
end

clear aenames
[nlag, aenames] = lagaena(parm); %names for Arima estimation printing

clear infm
infm = minfm(f, tr, mvx, tolf, maxit, nu0, jac, prt, chb, inc);

clear inft
inft = minft(fid, fh, wd, nd, scale);
%end of creating structures

if (trad < 0)
    tradt = 1; %flag for TD test
else
    tradt = 0;
end
if (easte < 0)
    eastet = 1; %flag for Easter test
else
    eastet = 0;
end
if (leapy < 0)
    leapyt = 1; %flag for Leap year test
else
    leapyt = 0;
end


%test for logs
nrout0 = 0;
nrout1 = 0;
nroutb = 0;
Y0 = [];
Y1 = [];
Yb = [];
if lam < 0
    miny = min(y);
    %   ctd=constant(ny,1,dr,ds,0,0,s);                  %generate a mean for the series
    ctd = constantx(ny, 1, dr, ds, dS, 0, 0, s, S); %generate a mean for the series
    if nreg > 0, YY = [ctd, Y(1:ny, :)];
    else YY = ctd;
    end
    % outlier detection is performed on both the original and the logged series
    est = 1;
    lpr = 0;
    x0 = cinest(y, YY, parm, est, ols, a, lpr, fid);
    [nrout1, nind1, tip1, Y1] = outlr(y, YY, parm, iout, infm, 0, x0, sp1, sp2, fmarqdt, ols, a);
    if miny <= 0
        lam = 1;
        parm.lam = lam; %series has negative values
        x0 = cinest(y, Y1, parm, est, ols, a, lpr, fid);
        xv = x0(pvar);
        xf = x0(pfix);
        if ~isempty(xv)
            xv = arimaopt(fmarqdt, fid, x0, xv, xf, y, Y1, parm, infm, lpr);
            x0(pvar) = xv;
        end
        [F, e, g1, M] = residual2x(x0, y, Y1, s, S, dr, ds, dS, p, ps, q, qs, qS);
    else
        y0 = log(y);
        x0 = cinest(y0, YY, parm, est, ols, a, lpr, fid);
        [nrout0, nind0, tip0, Y0] = outlr(y0, YY, parm, iout, infm, 0, x0, sp1, sp2, fmarqdt, ols, a);
        % nrout0
        % nrout1
        %both series should have the same number of outliers
        if (nrout0 < nrout1),
            nrout1 = nrout0;
            Y1 = Y0;
        elseif (nrout1 < nrout0)
            nrout0 = nrout1;
            Y0 = Y1;
            %if both series have the same number of outliers, first the one with less LS and
            %then the one with less TC is preferred
        else ...
                cls0 = 0;
            cls1 = 0;
            ctc0 = 0;
            ctc1 = 0;
            for i = 1:nrout0
                if strcmp(tip0(i), 'LS')
                    cls0 = cls0 + 1;
                end,
                if strcmp(tip1(i), 'LS')
                    cls1 = cls1 + 1;
                end,
                if strcmp(tip0(i), 'TC')
                    ctc0 = ctc0 + 1;
                end,
                if strcmp(tip1(i), 'TC')
                    ctc1 = ctc1 + 1;
                end,
            end
            if (cls0 < cls1)
                nrout1 = nrout0;
                Y1 = Y0;
                ctc1 = ctc0;
            elseif (cls1 < cls0)
                nrout0 = nrout1;
                Y0 = Y1;
                ctc0 = ctc1;
            end,
            if (ctc0 < ctc1)
                nrout1 = nrout0;
                Y1 = Y0;
            elseif ...
                    (ctc1 < ctc0), nrout0 = nrout1;
                Y0 = Y1;
            end,
        end
        x0 = cinest(y, Y1, parm, est, ols, a, lpr, fid);
        xv = x0(pvar);
        xf = xi(pfix);
        if ~isempty(xv)
            xv = arimaopt(fmarqdt, fid, x0, xv, xf, y, Y1, parm, infm, lpr);
            x0(pvar) = xv;
        end
        [F, e, g1, M] = residual2x(x0, y, Y1, s, S, dr, ds, dS, p, ps, q, qs, qS);
        Ff1 = F' * F;
        
        x0 = cinest(y0, Y0, parm, est, ols, a, lpr, fid);
        xv = x0(pvar);
        xf = xi(pfix);
        if ~isempty(xv)
            xv = arimaopt(fmarqdt, fid, x0, xv, xf, y0, Y0, parm, infm, lpr);
            x0(pvar) = xv;
        end
        [F, e, g0, M] = residual2x(x0, y0, Y0, s, S, dr, ds, dS, p, ps, q, qs, qS);
        Ff0 = F' * F;
        gmean = exp(sum(y0)/double(ny));
        Ff0 = Ff0 * gmean^2;
        
        % Ff0 and Ff1 are the criteria. There should be a relative difference between the two
        % greater than 2 per cent ?. If not, no transformation is selected
        if (Ff0 < .9875 * Ff1)
            lam = 0;
        else
            lam = 1;
        end
        %    if lam == 0
        %     format long g
        %     aa=abs((cFf1-cFf0)/cFf0)
        %     nrout0,nrout1
        %     cFf0
        %     cFf1
        %     lam
        %     pause
        %     format short
        %    end
        parm.lam = lam;
    end
    if lam == 0, y = log(yor);
    else y = yor;
    end %logs of the series
elseif (autmid == 1) || (tradt == 1) || (leapyt == 1) || (eastet == 1)
    %   ctd=constant(ny,1,dr,ds,0,0,s);                  %generate a mean for the series
    ctd = constantx(ny, 1, dr, ds, dS, 0, 0, s, S); %generate a mean for the series
    if nreg > 0, YY = [ctd, Y(1:ny, :)];
    else YY = ctd;
    end
    % outlier detection is performed before automatic model identification
    est = 1;
    lpr = 0;
    x0 = cinest(y, YY, parm, est, ols, a, lpr, fid);
    [nroutb, nindb, tipb, Yb] = outlr(y, YY, parm, iout, infm, 0, x0, sp1, sp2, fmarqdt, ols, a);
    x0 = cinest(y, Yb, parm, est, ols, a, lpr, fid);
    xv = x0(pvar);
    xf = x0(pfix);
    if ~isempty(xv)
        xv = arimaopt(fmarqdt, fid, x0, xv, xf, y, Yb, parm, infm, lpr);
        x0(pvar) = xv;
    end
    %   [F,e,gb,M]=residual2(x0,y,Yb,s,dr,ds,p,ps,q,qs,1);
    [F, e, gb, M] = residual2x(x0, y, Yb, s, S, dr, ds, dS, p, ps, q, qs, qS);
    %   nroutb
    %   gb
    %   return
end

outa.model.lam = lam;

if isfield(ser, 'usm')
    outa.y = y;
    return
end


%print series title and series
if pr == 1
    if freq == 1
        dateib = cal(1900, 1, 6); %print data in six columns
        inftb = inft;
        inftb.fh = 0; %without years and header
        fnameb = [dbname, ' (starting year is ', num2str(bg_year), ')'];
    else
        dateib = datei;
        inftb = inft;
        fnameb = dbname;
    end
    if nmiss > 0
        lamn = 1;
        prtser(fid, fnameb, yorm, y, ny, dateib, inftb, lamn);
        fnameb = [fnameb, ' with missing values filled'];
    end
    prtser(fid, fnameb, yor, y, ny, dateib, inftb, lam);
end


%correction for outlier and regression effects and test for trading day,
%leap year and easter
if (autmid == 1) || (trad ~= 0) || (leapy ~= 0) || (easte ~= 0)
    %correct series for regression effects
    %correct first for outlier effects. The previous estimates are used for this
    if nrout1 > 0
        if lam == 0
            ng = length(g0);
            ycoro = genycor(y, Y0(:, ...
                nreg+2:nreg+1+nrout0), ny, g0(ng-nrout0+1:ng)); ...
        else
        ng = length(g1);
        ycoro = genycor(y, Y1(:, nreg+2:nreg+1+nrout1), ...
            ny, g1(ng-nrout1+1:ng));
        end
    elseif nroutb > 0
        ng = length(gb);
        ycoro = genycor(y, Yb(:, nreg+2:nreg+1+nroutb), ny, gb(ng-nroutb+1:ng));
    else
        ycoro = y;
    end
    %correct for other regression effects. The previous estimates are used for this
    if nreg > 0
        [nY1, mY1] = size(Y1);
        [nYb, mYb] = size(Yb);
        if mY1 > 0,
            if lam == 0
                ycor = genycor(ycoro, Y0(:, 2:nreg+1), ny, g0(2:nreg+1));
            else
                ycor = genycor(ycoro, Y1(:, 2:nreg+1), ny, g1(2:nreg+1));
            end
        elseif mYb > 0
            ycor = genycor(ycoro, Yb(:, 2:nreg+1), ny, gb(2:nreg+1));
        else
            ycor = ycoro;
        end
    else
        ycor = ycoro;
    end
    %test for trading day, easter and leap year effects
    if (trad ~= 0) || (leapy ~= 0) || (easte ~= 0)
        YY = [];
        nnY = max(nY, ny+npr);
        ct = constantx(ny, 1, dr, ds, dS, 0, 0, s, S); %generate a mean for the series
        YY = [ct, YY];
        if (trad ~= 0)
            if length(tradval) > 1 %select TD variable by BIC
                parm = tradid(ycor, YY, infm, parm, ser, ols, a, tradval, fid, fmarqdt);
                trad = parm.trad;
            else
                trad = tradval(1);
            end
            if trad > 0
                Ytr = trade(bg_year, bg_per, nnY, trad, 0, [], freq);
                YYtr = Ytr(1:ny, :);
                Y = [Y, Ytr];
                [nY, mY] = size(Y);
                YY = [YY, YYtr];
                rnamesr = 1; ...
                    if isempty(rnamesrg)
                    rnamesrg = ['trad', num2str(1)];
                    else
                        rnamesrg = char(rnamesrg, ['trad', num2str(1)]);
                    end
                    for i = 2:trad, rnamesrg = char(rnamesrg, ['trad', num2str(i)]);
                    end
                    nreg = nreg + trad;
                    parm.nreg = nreg;
            else
                parm.trad = 0;
                trad = 0;
                tradt = 0;
            end
        end
        if (leapy ~= 0)
            if leapyt == 1 %test whether leap year is significant by BIC
                parm = leapid(ycor, YY, infm, parm, ser, ols, a, fid, fmarqdt);
                leap = parm.leap;
                leapy = leap;
            else
                leap = leapy;
            end
            if leap > 0
                Yle = genleap(bg_year, bg_per, nnY, freq);
                YYle = Yle(1:ny);
                Y = [Y, Yle];
                [nY, mY] = size(Y);
                YY = [YY, YYle];
                rnamesr = 1; ...
                    if isempty(rnamesrg)
                    rnamesrg = 'leap';
                    else
                        rnamesrg = char(rnamesrg, 'leap');
                    end
                    nreg = nreg + leap;
                    parm.nreg = nreg;
            else
                parm.leap = 0;
                leapy = 0;
                leapyt = 0;
            end
        end
        if (easte ~= 0)
            if length(durval) > 1 %select Easter duration by BIC
                parm = durid(ycor, YY, infm, parm, ser, ols, a, durval, fid, fmarqdt);
                dur = parm.dur;
            else %use Easter duration entered by the user
                dur = durval(1);
            end
            if dur > 0
                Yd = eastdate;
                Ye = east(bg_year, bg_per, nnY, dur, freq, Yd);
                Y = [Y, Ye];
                [nY, mY] = size(Y);
                rnamesr = 1;
                if isempty(rnamesrg)
                    rnamesrg = ['east(', num2str(dur), ')'];
                else
                    rnamesrg = char(rnamesrg, ['east(', num2str(dur), ')']);
                end
                nreg = nreg + 1;
                parm.east = 1;
                parm.nreg = nreg;
                easte = 1;
            else
                parm.east = 0;
                easte = 0;
                eastet = 0;
            end
        end
    end
    
    
    if nreg > 0
        %    ct=constant(ny,1,0,0,0,0,s);
        ct = constantx(ny, 1, 0, 0, 0, 0, 0, s, S);
        YY = [ct, Y(1:ny, :)];
        %estimate regression effects other than outliers
        est = 1;
        [yd, beta] = diffest(ycoro, YY, s, S, 0, 0, dS, est);
        %generate series corrected for all regression effects
        ycor = genycor(ycoro, YY(:, 2:end), ny, beta(2:end));
        %   ss=1:ny; plot(ss,y,ss,ycor,'r');
        %   return
    else
        ycor = ycoro;
    end
    
end


%incorporate a mean if flagm=1 and there is no automatic modeling
if (flagm == 1) && (autmid == 0)
    mY = size(Y, 2);
    %   ctd=constant(ny,1,dr,ds,0,0,s);
    ctd = constantx(ny+npr, 1, dr, ds, dS, 0, 0, s, S);
    if mY > 0
        if nmiss > 0
            Y = [Y(:, 1:nmiss), ctd, Y(:, nmiss+1:mY)];
            rnamesrg0 = rnamesrg;
            rnamesrg = char(rnamesrg0(1:nmiss, :), 'mean');
            if nmiss < mY
                rnamesrg = char(rnamesrg, rnamesrg0(nmiss+1:end, :));
            end
        else
            Y = [ctd, Y];
            if isempty(rnamesrg)
                rnamesrg = 'mean';
                rnamesr = 1;
            else
                rnamesrg = char('mean', rnamesrg);
            end
        end
    else
        Y = ctd;
        if isempty(rnamesrg)
            rnamesrg = 'mean';
            rnamesr = 1;
        else
            rnamesrg = char('mean', rnamesrg);
        end
    end
    [nY, mY] = size(Y);
    nreg = mY;
end
%automatic model identification
if (autmid == 1)
    %identify first dr and ds (degrees for regular and seasonal differencing)
    if (fixdif == 0)
        [nr1, ns1, dr, ds] = crcreg(ycor, freq, maxr);
        %    [freq,lam,dr,ds]
        %    [dr,ds]=curbic(ycor,freq);
        %    [dr ds],pause
    end
    parm.dr = dr;
    parm.ds = ds;
    
    %identify arma model for the differenced series
    %first correct the series for regression effects
    %   ct=constant(ny,1,dr,ds,0,0,s);                 %generate a mean for the series
    ct = constantx(ny, 1, dr, ds, dS, 0, 0, s, S); %generate a mean for the series
    if nreg > 0, YY = [ct, Y(1:ny, :)];
    else YY = ct;
    end
    %estimate regression effects other than outliers. Use the outlier corrected series
    est = 1;
    [yd, beta] = diffest(ycoro, YY, s, S, dr, ds, dS, est);
    %generate series corrected for all regression effects
    ycor = genycor(ycoro, YY, ny, beta);
    %then identify the model
    parm = armaid(ycor, parm, ols, a, maxpq, maxPQ);
    [s, S, dr, ds, dS, p, ps, q, qs, qS, ny, nreg, pfix, pvar] = imparm(parm);
    nr = p + ps + q + qs + qS; %number of arma parameters
    clear aenames
    [nlag, aenames] = lagaena(parm); %names for Arima estimation printing
    flagm = 0;
    %   ct=constant(ny,1,dr,ds,0,0,s);   %generate a mean for the differenced series
    ct = constantx(ny, 1, dr, ds, dS, 0, 0, s, S); %generate a mean for the series
    if nreg > 0, YY = [ct, Y(1:ny, :)];
    else YY = ct;
    end
    est = 1;
    x0 = cinest(y, YY, parm, est, ols, a, 0, fid);
    xv = x0(pvar);
    xf = x0(pfix);
    xv = arimaopt(fmarqdt, fid, x0, xv, xf, y, YY, parm, infm, 0);
    x0(pvar) = xv;
    % check seasonal underdifference
    if (ds == 0) && (ps == 1) && (qs == 1) && (fixdif == 0)
        Rs = -x0(p+ps);
        pps = p + ps;
        ppsq = pps + q;
        %   ns,nr
        %   x,Rs,hm1s,hm1sr,abs(x(pps)-x(ppsq+qs))
        alphamax = .499;
        can = .13;
        alpha2 = min(alphamax, .5-1/(ny^.90)); %.28  .51 .75
        hms = max(ny^(-alphamax), ny^(-alpha2));
        hm1s = 1 - hms;
        if ((Rs > hm1s) && (abs(x0(pps)-x0(ppsq+qs)) > can))
            parm.ds = 1;
            parm.ps = 0;
            ps = 0;
            nr = p + ps + q + qs + qS; %number of arma parameters
            parm.pvar = 1:nr;
            [s, S, dr, ds, dS, p, ps, q, qs, qS, ny, nreg, pfix, pvar] = imparm(parm);
            if pr == 1, prmod11x(fid, s, p, dr, q, ps, ds, qs, S, dS, qS, lam, flagm);
            end
            clear aenames
            [nlag, aenames] = lagaena(parm); %names for Arima estimation printing
            %   ct=constant(ny,1,dr,ds,0,0,s);   %generate a mean for the differenced series
            ct = constantx(ny, 1, dr, ds, dS, 0, 0, s, S); %generate a mean for the series
            if nreg > 0, YY = [ct, Y(1:ny, :)];
            else YY = ct;
            end
            est = 1;
            x0 = cinest(y, YY, parm, est, ols, a, 0, fid);
            xv = x0;
            xv = arimaopt(fmarqdt, fid, x0, xv, xf, y, YY, parm, infm, 0);
            x0 = xv;
        end
    end
    if (ds == 0) && (ps == 1) && (S == 0) && (fixdif == 0)
        [F, e, g, M, A, P, matsis] = residual2x(x0, y, YY, s, S, dr, ds, dS, p, ps, q, qs, qS);
        Ss = e' * e;
        Ff = F' * F;
        ndrs = ny - nmiss - dr - ds * s - dS * S; %residual sum of squares
        nparm = length(xv) + nreg;
        conp = Ss / (ndrs - nparm); %estimated sigma square
        sconp = sqrt(conp);
        ne = length(e);
        nv = length(xv);
        infr = rescomp(e, lag, nv, Ss, conp, sconp, Ff, ndrs, nreg);
        r = infr.r;
        orders = 1:floor(lag/s);
        nrs = ps + qs;
        np = min(2, length(orders));
        Qs = zeros(np, 1);
        pvalQs = zeros(np, 1);
        dfQs = zeros(np, 1);
        for i = 1:np
            Qs(i) = sum((r(s:s:i*s).^2)./(ne - s:-s:ne - i * s)') * ne * (ne + 2);
            dfQs(i) = max(0, i-nrs); %degrees of freedom
            if dfQs(i) > 0
                pvalQs(i) = 1 - gammp(dfQs(i)*.5, Qs(i)*.5);
            else
                pvalQs(i) = 1;
            end
        end;
        if np == 2
            Rs = -x0(p+ps);
            alphamax = .38;
            alpha2 = min(alphamax, .5-1/(ny^.90)); %.28  .51 .75
            hms = max(ny^(-alphamax), ny^(-alpha2));
            hm1s = 1 - hms;
            if ((Rs > hm1s) && (pvalQs(2) < .05))
                parm.ds = 1;
                parm.ps = 0;
                ps = 0;
                parm.qs = 1;
                qs = 1;
                nr = p + ps + q + qs + qS; %number of arma parameters
                parm.pvar = 1:nr;
                [s, S, dr, ds, dS, p, ps, q, qs, qS, ny, nreg, pfix, pvar] = imparm(parm);
                clear aenames
                [nlag, aenames] = lagaena(parm); %names for Arima estimation printing
                ct = constantx(ny, 1, dr, ds, dS, 0, 0, s, S); %generate a mean for the series
                if nreg > 0, YY = [ct, Y(1:ny, :)];
                else YY = ct;
                end
                est = 1;
                x0 = cinest(y, YY, parm, est, ols, a, 0, fid);
                xv = x0;
                xv = arimaopt(fmarqdt, fid, x0, xv, xf, y, YY, parm, infm, 0);
                x0 = xv;
            end
        end
    end
    % end of check seasonal underdifference
    
    % check regular underdifference
    if (dr == 0) && (p == 1) && (q == 1) && (fixdif == 0)
        Rr = -x0(p);
        pps = p + ps;
        ppsq = pps + q;
        %   ns,nr
        %   x,Rs,hm1s,hm1sr,abs(x(pps)-x(ppsq+qs))
        alphamax = .499;
        can = .13;
        if (s == 1)
            alpha2 = min(alphamax, .5-1/(ny^.55));
        else
            alpha2 = min(alphamax, .5-1/(ny^.90)); %.28  .51 .75
            alpha2r = min(alphamax, .5-1/(ny^.47)); %.33  .47
        end
        hms = max(ny^(-alphamax), ny^(-alpha2));
        hm1s = 1 - hms;
        if (s > 1)
            hmsr = max(ny^(-alphamax), ny^(-alpha2r));
            hm1sr = 1 - hmsr;
        end
        if (((Rr > hm1s) && (abs(x0(p)-x0(ppsq))) > can)) || ...
                (ds == 1) && ((Rr > hm1sr) && (abs(x0(p)-x0(ppsq)) > can))
            parm.dr = 1;
            parm.p = 0;
            p = 0;
            nr = p + ps + q + qs + qS; %number of arma parameters
            parm.pvar = 1:nr; %no fixed parameters when autmid = 1
            [s, S, dr, ds, dS, p, ps, q, qs, qS, ny, nreg, pfix, pvar] = imparm(parm);
            clear aenames
            [nlag, aenames] = lagaena(parm); %names for Arima estimation printing
            %   ct=constant(ny,1,dr,ds,0,0,s);   %generate a mean for the differenced series
            ct = constantx(ny, 1, dr, ds, dS, 0, 0, s, S); %generate a mean for the series
            if nreg > 0, YY = [ct, Y(1:ny, :)];
            else YY = ct;
            end
            est = 1;
            x0 = cinest(y, YY, parm, est, ols, a, 0, fid);
            xv = x0;
            xv = arimaopt(fmarqdt, fid, x0, xv, xf, y, YY, parm, infm, 0);
            x0 = xv;
        end
    end
    % end of check regular underdifference
    
    
    %   [F,e,beta,M]=residual2(x0,y,YY,s,dr,ds,p,ps,q,qs,1); nyd=ny-dr-ds*s; myd=length(beta);
    [F, e, beta, M] = residual2x(x0, y, YY, s, S, dr, ds, dS, p, ps, q, qs, qS);
    nyd = ny - dr - ds * s - dS * S;
    myd = length(beta);
    sg = e' * e / (nyd - myd);
    se = sqrt(diag(M)*sg);
    t = beta ./ se;
    %test whether mean is significant
    if (abs(t(1)) > 2.0) && (fixdif == 0)
        if nY > 0
            %     ct=constant(nY,1,dr,ds,0,0,s);
            ct = constantx(nY, 1, dr, ds, dS, 0, 0, s, S); %generate a mean for the series
        else
            %     ct=constant(ny+npr,1,dr,ds,0,0,s);
            ct = constantx(ny+npr, 1, dr, ds, dS, 0, 0, s, S); %generate a mean for the series
        end
        mY = size(Y, 2);
        if mY > 0
            if nmiss > 0
                Y = [Y(:, 1:nmiss), ct, Y(:, nmiss+1:mY)];
                rnamesrg0 = rnamesrg;
                if isempty(rnamesrg0(1:nmiss, :))
                    rnamesrg = 'mean';
                else
                    rnamesrg = char(rnamesrg0(1:nmiss, :), 'mean');
                end
                if nmiss < mY
                    rnamesrg = char(rnamesrg, rnamesrg0(nmiss+1:end, :));
                end
            else
                Y = [ct, Y];
                if isempty(rnamesrg)
                    rnamesrg = 'mean';
                else
                    rnamesrg = char('mean', rnamesrg);
                end
            end
        else
            Y = ct;
            if isempty(rnamesrg)
                rnamesrg = 'mean';
            else
                rnamesrg = char('mean', rnamesrg);
            end
        end
        [nY, mY] = size(Y);
        rnamesr = 1;
        nreg = nreg + 1;
        parm.nreg = nreg;
        flagm = 1;
        parm.flagm = flagm;
    end
    %test whether Easter effect is significant
    feast = 0;
    if (easte == 1)
        if (abs(t(end)) < 1.5 && (eastet == 1))
            [nY, mY] = size(Y);
            if (mY > 1), Y(:, end) = [];
            else Y = [];
            end
            [nY, mY] = size(Y);
            rnamesrg(end, :) = [];
            easte = 0;
            nreg = nreg - 1;
            parm.nreg = nreg;
            parm.east = 0;
            dur = 0;
            parm.dur = 0;
            feast = 1;
        end
    end
    %test whether leap year effect is significant
    fleap = 0;
    if (leapy == 1)
        if ((abs(t(end-feast)) < 1.5) && (leapyt == 1))
            [nY, mY] = size(Y);
            if (mY > 1), Y(:, end) = [];
            else Y = [];
            end
            [nY, mY] = size(Y);
            rnamesrg(end, :) = [];
            leapy = 0;
            nreg = nreg - 1;
            parm.nreg = nreg;
            parm.leap = 0;
            fleap = 1;
        end
    end
    %test whether trading day effect is significant
    if (trad > 0)
        if ((max(abs(t(end-feast-fleap-trad+1:end-feast-fleap))) < 1.5) ...
                && (tradt == 1))
            [nY, mY] = size(Y);
            if (mY > 1), Y(:, end-trad+1:end) = [];
            else Y = [];
            end
            [nY, mY] = size(Y);
            rnamesrg(end-trad+1:end, :) = [];
            nreg = nreg - trad;
            parm.nreg = nreg;
            trad = 0;
            parm.trad = 0;
        end
    end
    %print automatically identified model
    %   if pr == 1, prmod1(fid,s,p,dr,q,ps,ds,qs,lam,flagm); end
    if pr == 1, prmod1x(fid, s, p, dr, q, ps, ds, qs, S, dS, qS, lam, flagm);
    end
else
    %   if pr == 1, prmod2(fid,s,p,dr,q,ps,ds,qs,lam,flagm); end  %print default model
    if pr == 1, prmod2x(fid, s, p, dr, q, ps, ds, qs, S, dS, qS, lam, flagm);
    end %print default model
end

if (flagc0 == 1), iout.C = C;
end

%end of automatic model identification

outa.model.mean = flagm;
outa.model.p = p;
outa.model.d = dr;
outa.model.q = q;
outa.model.ps = ps;
outa.model.ds = ds;
outa.model.qs = qs;
if (S > 0)
    outa.model.S = S;
    outa.model.dS = dS;
    outa.model.qS = qS;
end
if (nreg > 0)
    outa.model.nreg = nreg;
end

%initial conditions
est = 1;
x0 = cinest(y, Y, parm, est, ols, a, pr, fid);


nrout = 0; %number of outliers
if out == 1 %outlier detection
    [nrout, nind, tip, Yo] = outlr(y, Y, parm, iout, infm, npr, x0, sp1, sp2, fmarqdt, ols, a);
    %Yo=[Y outliers], where outliers are the outlier variables detected
    iout.nrout = nrout;
    iout.nind = nind;
    iout.tip = tip;
    iout.Yo = Yo; %update iout structure
    Y = Yo;
    [nY, mY] = size(Y);
    fout = 0;
    %check whether there is an outlier at the end
    for i = 1:nrout, if nind(i) == ny, fout = i;
            break, end, end
    if fout > 0, nind(fout) = [];
        tip(fout, :) = [];
        Yo(:, fout) = [];
        nrout = nrout - 1;
        Y = Yo; ...
            iout.nrout = nrout;
        iout.nind = nind;
        iout.tip = tip;
        iout.Yo = Yo; %update iout structure
        [nY, mY] = size(Y);
    end
end
%print outliers and obtain new initial conditions
if nrout > 0
    ornames = ['out', num2str(1)];
    for i = 2:nrout
        ornames = char(ornames, ['out', num2str(i)]); %generate names for the outliers
    end
    iout.ornames = ornames;
    if pr == 1
        trtout(fid, iout, ser); %print outliers
        if fout > 0
            fprintf(fid, '%s\n', 'There is an outlier at the end of the series');
        end
        fprintf(fid, '\n');
    end
    est = 1;
    x0 = cinest(y, Y, parm, est, ols, a, pr, fid); %new initial conditions
elseif (out == 1) && (pr == 1)
    if omet == 1, met = 'Exact max. likelihood';
    else met = 'Hannan-Rissanen';
    end
    fprintf(fid, '\n%s %3.1f %s\n\n', 'No outliers detected (C = ', C, [', Method is ', ...
        met, '):']);
end
%end of outlier detection

%Arima estimation
xv = x0(pvar);
xf = xi(pfix);
if (nlestim == 1)
    xv = arimaopt(fmarqdt, fid, x0, xv, xf, y, Y, parm, infm, pr);
    x(pvar) = xv;
    x(pfix) = xf;
else
    x(pvar) = xv;
    x(pfix) = xf;
end
%end of Arima estimation


%check seasonal underdifference
if (ds == 0) && (S == 0) && (autmid == 1) && (fixdif == 0) && (freq > 1)
    [F, e, g, M, A, P, matsis] = residual2x(x, y, Y, s, S, dr, ds, dS, p, ps, q, qs, qS);
    Ss = e' * e;
    Ff = F' * F;
    ndrs = ny - nmiss - dr - ds * s - dS * S; %residual sum of squares
    nparm = length(xv) + nreg;
    conp = Ss / (ndrs - nparm); %estimated sigma square
    sconp = sqrt(conp);
    ne = length(e);
    nv = length(xv);
    infr = rescomp(e, lag, nv, Ss, conp, sconp, Ff, ndrs, nreg);
    r = infr.r;
    orders = 1:floor(lag/s);
    nrs = ps + qs;
    np = min(2, length(orders));
    Qs = zeros(np, 1);
    pvalQs = zeros(np, 1);
    dfQs = zeros(np, 1);
    for i = 1:np
        Qs(i) = sum((r(s:s:i*s).^2)./(ne - s:-s:ne - i * s)') * ne * (ne + 2);
        dfQs(i) = max(0, i-nrs); %degrees of freedom
        if dfQs(i) > 0
            pvalQs(i) = 1 - gammp(dfQs(i)*.5, Qs(i)*.5);
        else
            pvalQs(i) = 1;
        end
    end;
    if np == 2
        if ps == 1
            Rs = -x0(p+ps);
            alphamax = .38;
            alpha2 = min(alphamax, .5-1/(ny^.90)); %.28  .51 .75
            hms = max(ny^(-alphamax), ny^(-alpha2));
            hm1s = 1 - hms;
        else
            Rs = 1.;
            hm1s = 0.;
        end
        if ((Rs > hm1s) && (pvalQs(2) < .05))
            parm.ds = 1;
            parm.ps = 0;
            ps = 0;
            parm.qs = 1;
            qs = 1;
            nr = p + ps + q + qs + qS; %number of arma parameters
            parm.pvar = 1:nr; %no fixed parameters when autmid = 1
            [s, S, dr, ds, dS, p, ps, q, qs, qS, ny, nreg, pfix, pvar] = imparm(parm);
            if flagm == 1
                Y = [Y(:, 1:nmiss), Y(:, nmiss+2:end)];
                [~, mY] = size(Y);
                flagm = 0;
                rnamesrg = [rnamesrg(1:nmiss, :); rnamesrg(nmiss+2:end, :)];
                parm.flagm = 0;
            end
            if pr == 1, prmod11x(fid, s, p, dr, q, ps, ds, qs, S, dS, qS, lam, flagm);
            end
            ct = constantx(ny, 1, dr, ds, dS, 0, 0, s, S); %generate a mean for the series
            if nreg > 0, YY = [ct, Y(1:ny, :)];
            else YY = ct;
            end
            est = 1;
            x0 = cinest(y, YY, parm, est, ols, a, 0, fid);
            parm.pvar = pvar;
            clear aenames
            [nlag, aenames] = lagaena(parm); %names for Arima estimation printing
            xv = x0;
            xv = arimaopt(fmarqdt, fid, x0, xv, xf, y, Y, parm, infm, pr);
            x = xv;
        end
    end
end
% end of check seasonal underdifference
% check regular underdifference
if (dr == 0) && (p == 1) && (q == 1) && (fixdif == 0) && (autmid == 1)
    Rr = -x(p);
    pps = p + ps;
    ppsq = pps + q;
    %   ns,nr
    %   x,Rs,hm1s,hm1sr,abs(x(pps)-x(ppsq+qs))
    alphamax = .499;
    can = .13;
    if (s == 1)
        alpha2 = min(alphamax, .5-1/(ny^.55));
    else
        alpha2 = min(alphamax, .5-1/(ny^.90)); %.28  .51 .75
        alpha2r = min(alphamax, .5-1/(ny^.47)); %.33  .47
    end
    hms = max(ny^(-alphamax), ny^(-alpha2));
    hm1s = 1 - hms;
    if (s > 1)
        hmsr = max(ny^(-alphamax), ny^(-alpha2r));
        hm1sr = 1 - hmsr;
    end
    if (((Rr > hm1s) && (abs(x(p)-x(ppsq))) > can)) || ...
            (ds == 1) && ((Rr > hm1sr) && (abs(x(p)-x(ppsq)) > can))
        parm.dr = 1;
        parm.p = 0;
        p = 0;
        nr = p + ps + q + qs + qS; %number of arma parameters
        parm.pvar = 1:nr; %no fixed parameters when autmid = 1
        [s, S, dr, ds, dS, p, ps, q, qs, qS, ny, nreg, pfix, pvar] = imparm(parm);
        if flagm == 1
            Y = [Y(:, 1:nmiss), Y(:, nmiss+2:end)];
            [~, mY] = size(Y);
            flagm = 0;
            rnamesrg = [rnamesrg(1:nmiss, :); rnamesrg(nmiss+2:end, :)];
            parm.flagm = 0;
        end
        if pr == 1, prmod11x(fid, s, p, dr, q, ps, ds, qs, S, dS, qS, lam, flagm);
        end
        ct = constantx(size(y, 1), 1, dr, ds, dS, 0, 0, s, S); %generate a mean for the series
        if nreg > 0, YY = [ct, Y(1:ny, :)];
        else YY = ct;
        end
        est = 1;
        x0 = cinest(y, YY, parm, est, ols, a, 0, fid); %regression without inputs
        parm.pvar = pvar;
        clear aenames
        [nlag, aenames] = lagaena(parm); %names for Arima estimation printing
        xv = x0;
        xv = arimaopt(fmarqdt, fid, x0, xv, xf, y, Y, parm, infm, pr);
        x = xv;
    end
end
% end of check regular underdifference

% check regular overdifference
if (dr >= 1) && (q >= 1) && (fixdif == 0) && (autmid == 1)
    pps = p + ps;
    ppsq = pps + q;
    thr = [1, x(pps+1:ppsq)];
    aa = max(abs(roots(thr)));
    if aa > .98
        parm.dr = parm.dr - 1;
        if (p > 0) || (q > 1)
            parm.q = parm.q - 1;
            q = q - 1;
        end
        nr = p + ps + q + qs + qS; %number of arma parameters
        parm.pvar = 1:nr; %no fixed parameters when autmid = 1
        [s, S, dr, ds, dS, p, ps, q, qs, qS, ny, nreg, pfix, pvar] = imparm(parm);
        if flagm == 0
            ct = constantx(size(y, 1), 1, dr, ds, dS, 0, 0, s, S); %generate a mean for the series
            if nreg > 0, YY = [ct, Y(1:ny, :)];
            else YY = ct;
            end
        else
            YY = Y;
        end
        est = 1;
        x0 = cinest(y, YY, parm, est, ols, a, 0, fid); %regression without inputs
        parm.pvar = pvar;
        clear aenames
        [nlag, aenames] = lagaena(parm); %names for Arima estimation printing
        xv = x0;
        xv = arimaopt(fmarqdt, fid, x0, xv, xf, y, Y, parm, infm, pr);
        x = xv;
        if (flagm == 0)
            [F, e, beta, M] = residual2x(x, y, YY, s, S, dr, ds, dS, p, ps, q, qs, qS);
            nyd = ny - nmiss - dr - ds * s - dS * S;
            myd = length(beta);
            sg = e' * e / (nyd - myd);
            se = sqrt(diag(M)*sg);
            t = beta ./ se;
            %test whether mean is significant
            if (abs(t(1)) > 2.0) && (flagm == 0)
                if nY > 0
                    ct = constantx(nY, 1, dr, ds, dS, 0, 0, s, S); %generate a mean for the series
                else
                    ct = constantx(ny+npr, 1, dr, ds, dS, 0, 0, s, S); %generate a mean for the series
                end
                mY = size(Y, 2);
                if mY > 0
                    if nmiss > 0
                        Y = [Y(:, 1:nmiss), ct, Y(:, nmiss+1:mY)];
                        if isempty(rnamesrg(1:nmiss, :))
                            rnamesrg = 'mean';
                        else
                            rnamesrg = char(rnamesrg(1:nmiss, :), 'mean');
                        end
                        if nmiss < mY
                            rnamesrg = char(rnamesrg, rnamesrg(nmiss+1:end, :));
                        end
                    else
                        Y = [ct, Y];
                        if isempty(rnamesrg)
                            rnamesrg = 'mean';
                        else
                            rnamesrg = char('mean', rnamesrg);
                        end
                    end
                else
                    Y = ct;
                    if isempty(rnamesrg)
                        rnamesrg = 'mean';
                    else
                        rnamesrg = char('mean', rnamesrg);
                    end
                end
                [nY, mY] = size(Y);
                rnamesr = 1;
                nreg = nreg + 1;
                parm.nreg = nreg;
                flagm = 1;
                parm.flagm = flagm;
            end
        end
        if pr == 1, prmod11x(fid, s, p, dr, q, ps, ds, qs, S, dS, qS, lam, flagm);
        end
    end
end
% end of check regular overdifference

%
% get residuals
%
%  [F,e,g,M,A,P,matsis]=residual2(x,y,Y,s,dr,ds,p,ps,q,qs,1);
[F, e, g, M, A, P, matsis] = residual2x(x, y, Y, s, S, dr, ds, dS, p, ps, q, qs, qS);
Ss = e' * e;
Ff = F' * F;
ndrs = ny - nmiss - dr - ds * s - dS * S; %residual sum of squares
nparm = length(xv) + nreg;
conp = Ss / (ndrs - nparm); %estimated sigma square
sconp = sqrt(conp);

%check whether mean is significant
if (flagm == 1) && (fixdif == 0)
    seb = sqrt(diag(M)*conp);
    tb = g ./ seb; %standard errors and t-values
    if (abs(tb(nmiss+1)) < 2.0)
        mY = size(Y, 2);
        if nmiss > 0
            Y = [Y(:, 1:nmiss), Y(:, nmiss+1:mY)];
            rnamesrg = [rnamesrg(1:nmiss, :); rnamesrg(nmiss+1:mY, :)];
        else
            Y = Y(:, 2:mY);
            rnamesrg = rnamesrg(2:mY, :);
        end
        [nY, mY] = size(Y);
        rnamesr = 1;
        nreg = nreg - 1;
        parm.nreg = nreg;
        flagm = 0;
        parm.flagm = flagm;
        if pr == 1, prmod11x(fid, s, p, dr, q, ps, ds, qs, S, dS, qS, lam, flagm);
        end
        if ~isempty(xv)
            xv = arimaopt(fmarqdt, fid, x0, xv, xf, y, Y, parm, infm, pr);
            x(pvar) = xv;
        end
        %   [F,e,gg,M,A,P,matsis]=residual2(x,ycii,Y,s,dr,ds,p,ps,q,qs,1);
        [F, e, g, M, A, P, matsis] = residual2x(x, y, Y, s, S, dr, ds, dS, p, ps, q, qs, qS);
        nparm = length(xv) + nreg;
        Ss = e' * e;
        Ff = F' * F;
        ndrs = ny - nmiss - dr - ds * s - dS * S; %residual sum of squares
        conp = Ss / (ndrs - nparm); %estimated sigma square
        sconp = sqrt(conp);
    end
end


outa.model.result.pvar = pvar;
outa.model.result.pfix = pfix;
outa.model.result.xv = xv;
outa.model.result.xf = xf;
outa.model.result.x = x;
outa.model.result.nlestim = nlestim;
outa.model.result.mvx = mvx;
outa.model.result.fmarqdt = fmarqdt;

if (p > 0)
    outa.model.phi = x(1:p);
end
if (ps > 0)
    outa.model.phis = x(p+1:p+ps);
end
if (q > 0)
    pps = p + ps;
    outa.model.th = x(pps+1:pps+q);
end
if (qs > 0)
    pqps = p + q + ps;
    outa.model.ths = x(pqps+1:pqps+qs);
end
if (qS > 0)
    pqpsqs = p + q + ps + qs;
    outa.model.thS = x(pqpsqs+1:pqpsqs+qS);
end

if (out == 1)
    outa.model.nrout = nrout;
    outa.model.C = C;
    if (nrout > 0)
        outa.model.nind = nind;
        outa.model.tip = tip;
    end
end


outa.model.result.sigma2c = sconp;
outa.model.matsis = matsis;
outa.model.matsis.A = A;
outa.model.matsis.P = P;


%residual computations
if ~isempty(g) && olsres == 1
    e = matsis.olsres;
end
nv = length(xv);
infr = rescomp(e, lag, nv, Ss, conp, sconp, Ff, ndrs, nreg);

%computation of Pierce Qs value when there is one seasonality only
if (s > 1) && (S == 0)
    ne = length(e);
    r = infr.r;
    orders = 1:floor(lag/s);
    nrs = ps + qs;
    np = length(orders);
    Qs = zeros(np, 1);
    pvalQs = zeros(np, 1);
    dfQs = zeros(np, 1);
    for i = 1:np
        Qs(i) = sum((r(s:s:i*s).^2)./(ne - s:-s:ne - i * s)') * ne * (ne + 2);
        dfQs(i) = max(0, i-nrs); %degrees of freedom
        if dfQs(i) > 0
            pvalQs(i) = 1 - gammp(dfQs(i)*.5, Qs(i)*.5);
        else
            pvalQs(i) = 1;
        end
    end;
    infr.Qs = Qs;
    infr.dfQs = dfQs;
    infr.pvalQs = pvalQs;
    %Bell's seasonal dummies F-test
    nb = length(e);
    Ys = seasdm(size(outa.orig, 1), datei);
    Yb = Ys(end-size(e, 1)+1:end, :);
    nYb = size(Yb, 2);
    [bb, ~, er] = bmols(e, Yb);
    Ybb = Yb * bb;
    F = (Ybb' * Ybb / nYb) / (er' * er / (nb - nYb));
    pf = fdis_cdf(F, nYb, nb-nYb);
    %p-value for the F-test
    if (pf < .5)
        pH = 2 * pf;
    else
        pH = 2 * (1 - pf);
    end
    infr.pFtest = pH;
    infr.F = F;
    infr.dfn = nYb;
    infr.dfd = nb - nYb;
    %end of Bell's F-test
end

outa.model.resinf = infr;

%compute standard errors and t-values for the regression part
if mY > 0
    seb = sqrt(diag(M)*conp);
    tb = g ./ seb; %standard errors and t-values
    if nmiss > 0
        Interp = y(idxn) - g(1:nmiss);
        tb(1:nmiss) = NaN(nmiss, 1);
        yinterp = y;
        yinterp(idxn) = Interp;
        sInterp = seb(1:nmiss);
        outa.model.yinterp = yinterp;
        outa.model.Interp = Interp;
        outa.model.sInterp = sInterp;
        %obtain interpolations in the original scale using the log-normal
        %distribution
        oInterp = Interp;
        osInterp = sInterp;
        if lam == 0
            for i = 1:nmiss
                oInterp(i) = exp(Interp(i)+(sInterp(i)^2)/double(2.));
                osInterp(i) = exp(double(2.)*Interp(i)+sInterp(i)^2) ...
                    * (exp(sInterp(i)^2) - double(1.));
            end
            oyinterp = yor;
            oyinterp(idxn) = oInterp;
            outa.model.oyinterp = oyinterp;
            outa.model.oInterp = oInterp;
            outa.model.osInterp = osInterp;
        end
    end
    Mbeta = [g, seb, tb];
    outa.model.hb = g;
    outa.model.Mb = M;
    outa.model.Y = Y;
    outa.model.seb = seb;
    outa.model.tb = tb;
end

% rnamesrg
%compute regression effects other than outliers
Yrg = [];
if nreg > 0
    nypnpr = ny + npr;
    t = 1:nypnpr;
    mreg = nreg - nmiss;
    Yrg = zeros(nypnpr, mreg);
    for i = nmiss + 1:nreg
        Yrg(:, i-nmiss) = Y(t, i) * g(i); %generate matrix of regression effects
    end
else
    mreg = 0;
end

%compute outlier effects
Youtg = [];
if nrout > 0
    nypnpr = ny + npr;
    t = 1:nypnpr;
    Youtg = zeros(nypnpr, nrout);
    for i = 1:nrout
        %    Youtg=[Youtg Y(t,nreg+i)*g(nreg+i)];          %generate matrix of outlier effects
        Youtg(:, i) = Y(t, nreg+i) * g(nreg+i); %generate matrix of outlier effects
    end
end

outa.model.Yrg = Yrg;
outa.model.Youtg = Youtg;

%
%standard errors
if (nlestim == 1)
    if fjac == 0
        %standard errors via second derivatives of log likelihood
        %
        %the following line changed on 28-3-2017
        %
        %      xt=x';
        xt = x(pvar)';
        %     H=fdhess('logF',xt,'fstarima1',y,Y,parm,infm,xf);
        %     SS=pinv(H/2)/(ny-dr-ds*s-nr-nreg);
        H = fdhess('logF', xt, f, y, Y, parm, infm, xf);
        SS = pinv(H/2) / (ny - dr - ds * s - dS * S - nparm - ninput - nreg);
        se = sqrt(abs(diag(SS)))';
    else
        %standard errors via the jacobian
        smvx = infm.mvx;
        infm.mvx = 0;
        %       [f,junk] = fstarima1(x,y,Y,parm,infm,xf);
        %       [J,junk] = fdjac2(infm,x,f,y,Y,parm,infm,xf);
        %       infm.mvx=smvx;
        %       RR=qr(J);
        %       R=RR(1:nparm,1:nparm); Ri=pinv(R);
        xp = x(pvar);
        f = fasttf(xp, y, Y, parm, infm, xf);
        [J, junk] = fdjac2(infm, xp, f, y, Y, parm, infm, xf);
        infm.mvx = smvx;
        RR = qr(J);
        np = length(xp);
        R = RR(1:np, 1:np);
        Ri = pinv(R);
        SS = conp * (Ri * Ri');
        se = sqrt(abs(diag(SS)))';
    end
    tt = zeros(size(x));
else
    se = NaN(size(xv));
    tt = NaN(size(x));
    SS = [];
end

%t-values
%
%the following lines changed on 28-3-2017
%
see = zeros(size(x));
see(pvar) = se;
see(pfix) = NaN;
tt(pfix) = NaN;
if ~isempty(pvar) && (nlestim == 1)
    tt(pvar) = x(pvar) ./ se;
end
z = [x', see', tt', nlag];
%end of change

if ~isempty(x)
    outa.model.se = se;
    outa.model.tt = tt;
    outa.model.SS = SS;
end


if pr == 1
    %print Arima estimation results
    clear in
    in.cnames = char('  Estimate', 'Std. Error', '   T-ratio', 'Lag');
    in.rnames = aenames;
    in.fmt = char('%12.4f', '%12.4f', '%12.4f', '%2.0f');
    in.fid = fid;
    mprint(z, in);
    fprintf(fid, '\nResidual standard error:%11.4f\n\n', sconp);
    %print roots of ar and ma polynomials
    aenames = lagaenar(parm);
    z = rootsarma(x, parm);
    clear in
    in.cnames = char('  Real p.', 'Imag. p.', ' Modulus', 'Argument', ' Period');
    in.rnames = aenames;
    in.fmt = char('%12.4f', '%12.4f', '%12.4f', '%6.4f', '%6.4f');
    in.fid = fid;
    mprint(z, in);
    fprintf(fid, '\n');
    %print correlations of the estimates
    clear in
    DD = diag(sqrt(abs(diag(SS))));
    CC = (DD \ SS) / DD;
    clear aenames
    [~, aenames] = lagaena(parm); %names for Arima estimation parameters
    in.cnames = aenames(2:end, :);
    in.cnames = in.cnames(pvar, :);
    in.rnames = aenames;
    in.rnames = [in.rnames(1, :); in.cnames];
    in.fmt = '%12.4f';
    in.fid = fid;
    if ~isempty(CC)
        fprintf(fid, '\nCorrelations of the estimates:\n');
        mprint(CC, in);
        fprintf(fid, '\n');
    end
    %print regression variables, including outliers
    if mY > 0
        clear in
        in.cnames = char('  Estimate', 'Std. Error', '   T-ratio');
        rnames = char('Parameter');
        if rnamesr == 1
            if ~isempty(rnamesrg)
                rnames = char(rnames, rnamesrg);
            end
        else
            if isempty(rnamesrg)
                for i = 1:nreg, rnames = char(rnames, ['reg', num2str(i)]);
                end
            else
                rnames = char(rnames, rnamesrg);
            end
        end
        if nrout > 0
            rnames = char(rnames, ornames);
        end
        in.rnames = rnames;
        in.fmt = char('%12.5f', '%12.5f', '%8.2f');
        in.fid = fid;
        mprint(Mbeta, in);
        fprintf(fid, '\n');
        %print correlations of the regression estimates
        fprintf(fid, '\nCorrelations of the regression estimates:\n');
        clear in
        DDr = diag(sqrt(abs(diag(M))));
        CCr = (DDr \ M) / DDr;
        in.cnames = rnames(2:end, :);
        in.rnames = rnames;
        in.fmt = '%12.4f';
        in.fid = fid;
        mprint(CCr, in);
        fprintf(fid, '\n');
        
        
        if nmiss > 0
            clear in
            in.cnames = char('  Estimate', 'Std. Error');
            rnamesrgi = ['interp. ', num2str(idxn(1))];
            for i = 2:nmiss
                rnamesrgi = char(rnamesrgi, ['interp. ', num2str(idxn(i))]);
            end
            rnames = char('Interpolated value', rnamesrgi);
            in.rnames = rnames;
            in.fmt = char('%12.5f');
            in.fid = fid;
            mprint([Interp, sInterp], in);
            fprintf(fid, '\n');
            if (lam == 0)
                rnames = char('Interpolated value (levels)', rnamesrgi);
                in.rnames = rnames;
                mprint([oInterp, osInterp], in);
                fprintf(fid, '\n');
            end
        end
    end
    %print residual diagnostics
    printres(fid, infr);
    % print summary
    prsummry(ii, ny, nreg0, fidr, dbname, parm, iout);
end


%Forecasts
if npr > 0
    Z = matsis.Z;
    T = matsis.T;
    H = matsis.H;
    [pry, spry] = predt(ny, npr, sconp, Y, Z, T, H, A, P, g, M);
    outa.model.npr = npr;
    outa.model.pry = pry;
    outa.model.spry = spry;
    
    %obtain forecasts in the original scale using the log-normal
    %distribution
    opry = pry;
    ospry = spry;
    if lam == 0
        for i = 1:npr
            opry(i) = exp(pry(i)+(spry(i)^2)/double(2.));
            ospry(i) = exp(double(2.)*pry(i)+spry(i)^2) * (exp(spry(i)^2) - double(1.));
        end
        outa.model.opry = opry;
        outa.model.ospry = ospry;
    end
    
    if pr == 1
        %print forecasts
        fprintf(fid, '\n');
        clear in
        in.fid = fid;
        if lam == 1
            in.cnames = char('Obs.      ', 'Forecasts ', 'Std. Error');
            in.fmt = char('%5.0f', '%12.4f', '%12.4f');
            mprint([[ny + 1:ny + npr]', pry, spry], in);
        else
            in.cnames = char('Obs.         ', 'Forecasts    ', 'Std. Error   ', ...
                'For. (levels)', 'Std. Error   ');
            in.fmt = char('%5.0f', '%12.4f', '%12.4f', '%19.4f', '%19.4f');
            mprint([[ny + 1:ny + npr]', pry, spry, opry, ospry], in);
        end
        fprintf(fid, '\n');
    end
end

% print elapsed time
time = toc;
if pr == 1, fprintf(fid, 'Elapsed time: %8.2f\n\n', time);
end

%close text external file
if ((fid ~= 1) && (pr == 1)), fclose(fid);
end

%write numerical results
if pr == 1
    tabla = [];
    if npr > 0
        yorx = [yor; opry];
        tabla = [tabla, yorx];
        if lam == 0
            yx = [y; pry];
            tabla = [tabla, yx];
        end
    else
        yorx = yor;
        tabla = [tabla, yorx];
        if lam == 0
            yx = y;
            tabla = [tabla, yx];
        end
    end
    [pY, qY] = size(Y);
    if qY > 0
        tabla = [tabla, Y * g]; %series of regression effects
    end
    nname = ['results', filesep, dbname, '.num'];
    str = ['save ''', nname, ''' tabla -ascii -double'];
    eval(str);
end

%plots
outa.gft = gft;
if gft >= 1
    pathc = pwd; %current directory
    if (gft > 1)
        if exist([pathc, filesep, 'graphs'], 'dir') ~= 7
            mkdir(pathc, 'graphs'); %create directory graphs
        end
    end
    if ~isempty(ycii)
        infr.ycii = ycii;
    end
    dbname = strrep(dbname, '_', '\_');
    %Youtg: outlier effects, Yrg: regression effects.
    plotres(y, Y, g, yor, datei, cw, dbname, gft, nrout, Youtg, mreg, Yrg, infr, s, lam);
    if (nmiss > 0)
        figure
        vnames = char('Interpolated series');
        tsplot(yinterp, datei, vnames);
        disp('press any key to continue');
        pause;
    end
    if npr > 0
        %plot forecasts
        figure
        tt = ny - npr + 1:ny;
        y1 = [y(tt); pry];
        y2 = [y(tt); pry + cw * spry];
        y3 = [y(tt); pry - cw * spry];
        %    ny1=length(y1);
        %    miny=[min(y1); min(y2); min(y3)]; maxy=[max(y1); max(y2); max(y3)];
        %    t=1:ny1;
        %    npt1=min(miny); npt2=max(maxy); step=1.01*(npt2-npt1)/100;
        %    c2=npt1:step:npt2*1.01; c1=npr*ones(size(c2));
        %    plot(t,y1,'r-',t,y2,'-.',t,y3,'--',c1,c2);
        t = -npr + 1:npr;
        vnames = ['Upper 95% band '; ...
            'Forecast       '; ...
            'Lower 95% band '];
        plot(t, y2, '-.', t, y1, 'r-', t, y3, '--');
        legend(vnames);
        set(gca, 'tickdir', 'in');
        set(gca, 'xcolor', 'k');
        set(gca, 'GridLineStyle', ':');
        set(gca, 'Xgrid', 'on');
        axis tight;
        
        
        if lam == 0
            title(['Forecasts of series ', dbname, ' (in logs)'])
            disp('press any key to continue');
            pause;
            figure
            y1 = [yor(tt); opry];
            %        miny=[min(y1)]; maxy=[max(y1)];
            %        t=1:ny1;
            %        npt1=min(miny); npt2=max(maxy); step=1.01*(npt2-npt1)/100;
            %        c2=npt1:step:npt2*1.01; c1=npr*ones(size(c2));
            %        plot(t,y1,'r-',c1,c2);
            t = -npr + 1:npr;
            plot(t, y1, 'r-');
            set(gca, 'tickdir', 'in');
            set(gca, 'xcolor', 'k');
            set(gca, 'GridLineStyle', ':');
            set(gca, 'Xgrid', 'on');
            axis tight;
            title(['Forecasts of series ', dbname])
        else
            title(['Forecasts of series  ', dbname])
        end
        disp('press any key to continue');
        pause;
        if (s > 1)
            %plot forecasts in rates of grothw
            figure
            yt = tasa([yor; opry], s) * 100;
            nyt = length(yt);
            tt = nyt - 2 * npr + 1:nyt;
            %     miny=min(yt(tt)); maxy=max(yt(tt));
            %     t=1:2*npr;
            %     npt1=min(miny); npt2=max(maxy); step=1.01*(npt2-npt1)/100;
            %     c2=npt1:step:npt2*1.01; c1=npr*ones(size(c2));
            %     plot(t,yt(tt),'r-',c1,c2);
            t = -npr + 1:npr;
            plot(t, yt(tt), 'r-');
            set(gca, 'tickdir', 'in');
            set(gca, 'xcolor', 'k');
            set(gca, 'GridLineStyle', ':');
            set(gca, 'Xgrid', 'on');
            axis tight;
            title(['Forecasts of original series ', dbname, ' (rates of growth in percentage)'])
        end
    end
    disp('press any key to continue');
    pause;
    close all
end

outa.tfmodel = [];

% profile off
