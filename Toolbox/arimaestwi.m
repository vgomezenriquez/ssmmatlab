function outa = arimaestwi(dbname, ser, fidr, ii)
%
% function to identify, estimate and forecast a transfer function model for
% one series. The method used for automatic model identificaiton is
% described in Gomez (2009), "Transfer Function Model Identification",
% Boletin de Estadistica e Investigacion Operativa, 25, pp. 109-115.
%
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
%              the finite linear approximation to the input filters.
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
%     tfmodel: structure with transfer function model information.
%       matsis: a structure containing the state space form of the model
%               correponding to the filtered inputs
%       result: a structure containing estimation results
%         nreg: the number of regression variables
%          Yrg: array containing the regression variables that are not
%               outliers
%        Youtg: array containing the outlier variables
%          yci: output corrected by filtered inputs
%        tford: a three column array in which the i-th row has three
%               numbers corresponding to the delay, the numerator degree
%               and the denominator degree of the i-th input filter
%          phi: an array containing the regular AR polynomial coefficients
%         phis: an array containing the seasonal AR polynomial coefficients
%           th: an array containing the regular MA polynomial coefficients
%          ths: an array containing the seasonal MA polynomial coefficients
%          omg: a cell array containing the numerators of the input filters
%          del: a cell array containing the denominators of the input
%               filters
%       resinf: a structure containing information about the residuals
%           se: array containing the standard errors of the estimates
%           tt: array containing the t-values of the estimates
%          npr: number of forecasts
%         dpry: array containing the forecasts of the output corrected by
%               the filtered inputs
%        dspry: array containing the standard errors of the forecasts of
%               the output corrected by the filtered inputs
%          Yin: array containing the input variables
%      modpred: a multiple structure containing the input forecasts if any.
%               the forecasts for each input are given in field .pred
%     modinput: a multiple structure containing the models for the inputs
%               if any (.mod = 0, no model; .mod =1, there is model). The
%               model for each input has fields .alpha, .phi, .theta,
%               .sigma2
%            y: output series in the transformed scale
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
else
    rnamesrgm = [];
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


if (prelivar == 1)
    return
end
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
    if (ny ~= nY)
        disp('regression variables do not have the same length as series')
    end
    meanY = mean(Y);
    Y = [Y; zeros(npr, mY)];
    [nY, mY] = size(Y);
    nreg = mY;
    for i = 1:npr
        Y(ny+i, :) = meanY;
    end
    if pr == 1
        fprintf(fid, '%s\n', 'Regression matrix is extended with mean values for prediction');
        if (ser.ninput > 0)
            fprintf(fid, '%s\n', 'in the first round of the program');
        end
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
    C = C0;
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
if ((lam < 0) || (autmid == 1))
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
    ctd = constant(ny, 1, dr, ds, 0, 0, s); %generate a mean for the series
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
        %    [F,e,g1,M]=residual2(x0,y,Y1,s,dr,ds,p,ps,q,qs,1);
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
%         xv = x0(pvar);
%         xf = xi(pfix);
%         if ~isempty(xv)
%             xv = arimaopt(fmarqdt, fid, x0, xv, xf, y, Y1, parm, infm, lpr);
%             x0(pvar) = xv;
%         end
        [F, e, g1, M] = residual2x(x0, y, Y1, s, S, dr, ds, dS, p, ps, q, qs, qS);
        Ff1 = F' * F; %residual sum of squares
        x0 = cinest(y0, Y0, parm, est, ols, a, lpr, fid);
%         xv = x0(pvar);
%         xf = xi(pfix);
%         if ~isempty(xv)
%             xv = arimaopt(fmarqdt, fid, x0, xv, xf, y0, Y0, parm, infm, lpr);
%             x0(pvar) = xv;
%         end
        [F, e, g0, M] = residual2x(x0, y0, Y0, s, S, dr, ds, dS, p, ps, q, qs, qS);
        Ff0 = F' * F; %residual sum of squares
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
elseif (autmid == 1) | (tradt == 1) | (leapyt == 1) | (eastet == 1)
    ctd = constant(ny, 1, dr, ds, 0, 0, s); %generate a mean for the series
    if nreg > 0, YY = [ctd, Y(1:ny, :)];
    else, YY = ctd;
    end
    % outlier detection is performed before automatic model identification
    est = 1;
    lpr = 0;
    x0 = cinest(y, YY, parm, est, ols, a, lpr, fid);
    xv = x0(pvar);
    xf = x0(pfix);
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
        nym = length(yorm);
        prtser(fid, fnameb, yorm, y, nym, dateib, inftb, lamn);
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
        ct = constant(ny, 1, dr, ds, 0, 0, s); %generate a mean for the differenced series
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
                rnamesr = 1;
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
        ct = constant(ny, 1, 0, 0, 0, 0, s);
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
            nmiss1 = nmiss - npmiss;
            Y = [Y(:, 1:nmiss1), ctd, Y(:, nmiss1+1:mY)];
            rnamesrg0 = rnamesrg;
            rnamesrg = char(rnamesrg0(1:nmiss1, :), 'mean');
            if nmiss < mY
                rnamesrg = char(rnamesrg, rnamesrg0(nmiss1+1:end, :));
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
    if fixdif == 0
        [nr1, ns1, dr, ds] = crcreg(ycor, freq, maxr);
        %    [freq,lam,dr,ds]
        %    [dr,ds]=curbic(ycor,freq);
        %    [dr ds],pause
    end
    parm.dr = dr;
    parm.ds = ds;
    
    %identify arma model for the differenced series
    %first correct the series for regression effects
    ct = constant(ny, 1, dr, ds, 0, 0, s); %generate a mean for the series
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
    ct = constant(ny, 1, dr, ds, 0, 0, s); %generate a mean for the differenced series
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
            parm.pvar = 1:nr; %no fixed parameters when autmid = 1
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
        ndrs = ny - nmiss + npmiss - dr - ds * s - dS * S; %residual sum of squares
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
                parm.pvar = 1:nr; %no fixed parameters when autmid = 1
                [s, S, dr, ds, dS, p, ps, q, qs, qS, ny, nreg, pfix, pvar] = imparm(parm);
                clear aenames
                [nlag, aenames] = lagaena(parm); %names for Arima estimation printing
                ct = constantx(ny, 1, dr, ds, dS, 0, 0, s, S); %generate a mean for the series
                if nreg > 0 
                    YY = [ct, Y(1:ny, :)];
                else
                    YY = ct;
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
    if (dr <= 1) && (p >= 1) && (fixdif == 0)
        aa = roots([1, x0(1:p)]);
        [Rr,ii] = max(real(aa));  
        alphamax = .499;
        if (s == 1)
            alpha2 = min(alphamax, .5-1/(ny^.55));
        else
            alpha2 = min(alphamax, .5-1/(ny^.90)); %.28  .51 .75 .90
            alpha2r = min(alphamax, .5-1/(ny^.38)); %.33 .38 .47
        end
        hms = max(ny^(-alphamax), ny^(-alpha2));
        hm1s = 1 - hms;
        if (s > 1)
            hmsr = max(ny^(-alphamax), ny^(-alpha2r));
            hm1sr = 1 - hmsr;
        end
        if ((Rr > hm1s) && (abs(imag(aa(ii))) < hms)  || ((ds == 1)...
                && (Rr > hm1sr)) && (abs(imag(aa(ii))) < hmsr) )
            parm.dr = dr + 1;
            parm.p = 0;
            p = 0;
            nr = p + ps + q + qs + qS; %number of arma parameters
            if (nr == 0)
                q = 1;
                parm.q = q;
                nr = 1;
            end
            parm.pvar = 1:nr; %no fixed parameters when autmid = 1
            [s, S, dr, ds, dS, p, ps, q, qs, qS, ny, nreg, pfix, pvar] = imparm(parm);
            clear aenames
            [nlag, aenames] = lagaena(parm); %names for Arima estimation printing
            %   ct=constant(ny,1,dr,ds,0,0,s);   %generate a mean for the differenced series
            ct = constantx(ny, 1, dr, ds, dS, 0, 0, s, S); %generate a mean for the series
            if nreg > 0 
                YY = [ct, Y(1:ny, :)];
            else
                YY = ct;
            end
            est = 1;
            x0 = cinest(y, YY, parm, est, ols, a, 0, fid);
            xv = x0;
            xv = arimaopt(fmarqdt, fid, x0, xv, xf, y, YY, parm, infm, 0);
            x0 = xv;
        end
    end
    % end of check regular underdifference
    %   [F,e,beta,M]=residual2(x0,y,YY,s,dr,ds,p,ps,q,qs,1);
    [F, e, beta, M] = residual2x(x0, y, YY, s, S, dr, ds, dS, p, ps, q, qs, qS);
    nyd = ny - nmiss + npmiss - dr - ds * s - dS * S;
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
                nmiss1 = nmiss - npmiss;
                Y = [Y(:, 1:nmiss1), ct, Y(:, nmiss1+1:mY)];
                rnamesrg0 = rnamesrg;
                if isempty(rnamesrg0(1:nmiss1, :))
                    rnamesrg = 'mean';
                else
                    rnamesrg = char(rnamesrg0(1:nmiss1, :), 'mean');
                end
                if nmiss1 < mY
                    rnamesrg = char(rnamesrg, rnamesrg0(nmiss1+1:end, :));
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
        if ser.ninput > 0
            fprintf(fid, '%s\n', 'The following observation numbers correspond to the first round series');
            sero = ser;
            sero.bg_year = bg_year;
            sero.bg_per = bg_per;
            trtout(fid, iout, sero); %print outliers
        else
            trtout(fid, iout, ser); %print outliers
        end
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
    x0 = x;
else
    x(pvar) = xv;
    x(pfix) = xf;
    x0 = x;
end
%end of Arima estimation

outa.model.mean = flagm;
outa.model.p = p;
outa.model.d = dr;
outa.model.q = q;
outa.model.ps = ps;
outa.model.ds = ds;
outa.model.qs = qs;
if (nreg > 0)
    outa.model.nreg = nreg;
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

if (out == 1)
    outa.model.nrout = nrout;
    outa.model.C = C;
    if (nrout > 0)
        outa.model.nind = nind;
        outa.model.tip = tip;
    end
end

%
% get residuals
%
%  [F,e,g,M,A,P,matsis]=residual2(x,y,Y,s,dr,ds,p,ps,q,qs,1);
[F, e, g, M, A, P, matsis] = residual2x(x, y, Y, s, S, dr, ds, dS, p, ps, q, qs, qS);
Ss = e' * e;
Ff = F' * F;
ndrs = ny - nmiss + npmiss - dr - ds * s - dS * S; %residual sum of squares
nparm = length(xv) + nreg;
conp = Ss / (ndrs - nparm); %estimated sigma square
sconp = sqrt(conp);
gm = g;

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
        if npmiss > 0
            if npmiss1 <= nmiss
                nmnp = nmiss - npmiss;
                tb(1:nmnp) = NaN(nmnp, 1);
            end
        else
            tb(1:nmiss) = NaN(nmiss, 1);
        end
    end
    Mbeta = [g, seb, tb];
    outa.model.hb = g;
    outa.model.Mb = M;
    outa.model.Y = Y;
    outa.model.seb = seb;
    outa.model.tb = tb;
    
    
    %backward elimination for transfer function identification
    if (backwd == 1) && (myl > 0) %myl is the number of lags.
        flag = 1;
        %inreg is the number of initial regression variables, included the mean
        inreg = initreg + flagm;
        [m0, n0] = size(nelinput); %nelinput is the number of weights for each input variable
        firstel0 = [];
        cont = inreg;
        for i = 1:m0
            firstel0 = [firstel0; abs(tb(cont+1))];
            cont = cont + nelinput(i);
        end
        %cnelinput is the number of weights for each remaining input variable
        cnelinput = nelinput;
        while flag == 1
            [mi, ni] = size(cnelinput);
            firstel = [];
            firstelin = [];
            cont = inreg;
            for i = 1:mi
                firstel = [firstel; abs(tb(cont+1))];
                firstelin = [firstelin; cont + 1];
                cont = cont + cnelinput(i);
            end
            %     firstel
            %     firstelin
            cont0 = 0;
            for i = 1:m0
                if firstel0(i) < 1.e11
                    cont0 = cont0 + 1;
                    firstel0(i) = firstel(cont0);
                end
            end
            [mtv0, j0] = min(firstel0);
            [mtv, j] = min(firstel);
            i = firstelin(j);
            [nY, mY] = size(Y);
            %     firstel0
            %     mtv,j,i,myl
            %     nelinput
            %     cnelinput
            %     pause
            if (mtv < Cb) && (i <= myl+inreg)
                nreg = nreg - 1;
                myl = myl - 1;
                nelinput(j0) = nelinput(j0) - 1;
                if nelinput(j0) == 0, firstel0(j0) = 1.e11;
                end
                cnelinput(j) = cnelinput(j) - 1;
                if cnelinput(j) == 0
                    naux = [];
                    for k = 1:mi
                        if k ~= j
                            naux = [naux; cnelinput(k)];
                        end
                    end
                    cnelinput = naux;
                end
                if i > 1
                    Y = Y(:, [1:i - 1, i + 1:mY]);
                    if rnamesr == 1
                        [nrg, mr] = size(rnamesrg);
                        rnamesrg = rnamesrg([1:i - 1, i + 1:nrg], :);
                    end
                    %      [beta,tv]=btval([],[y Y]);  yc=y-Y*beta;
                    xv = x0(pvar);
                    xf = x0(pfix);
                    if ~isempty(xv)
                        xv = arimaopt(fmarqdt, fid, x0, xv, xf, y, Y, parm, infm, 0);
                        x(pvar) = xv;
                        x(pfix) = xf;
                    else
                        x(pvar) = xv;
                        x(pfix) = xf;
                    end
                    if (p > 0)
                        outa.model.phi = x(1:p);
                    end
                    if (q > 0)
                        outa.model.th = x(p+1:p+q);
                    end
                    if (ps > 0)
                        pq = p + q;
                        outa.model.phis = x(pq+1:pq+ps);
                    end
                    if (qs > 0)
                        pqps = p + q + ps;
                        outa.model.ths = x(pqps+1:pqps+qs);
                    end
                    %       [F,e,g,M,A,P,matsis]=residual2(x,y,Y,s,dr,ds,p,ps,q,qs,1);
                    [F, e, g, M, A, P, matsis] = residual2x(x, y, Y, s, S, dr, ds, dS, p, ps, q, qs, qS);
                    Ss = e' * e;
                    Ff = F' * F;
                    ndrs = ny - nmiss + npmiss - dr - ds * s - dS * S; %residual sum of squares
                    nparm = length(xv) + nreg;
                    conp = Ss / (ndrs - nparm); %estimated sigma square
                    sconp = sqrt(conp);
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
                    seb = sqrt(diag(M)*conp);
                    tb = g ./ seb; %standard errors and t-values
                    %       if nmiss > 0
                    %        tb(1:nmiss)=NaN(nmiss,1);
                    %       end
                    if nmiss > 0
                        if npmiss > 0
                            if npmiss1 <= nmiss
                                nmnp = nmiss - npmiss;
                                tb(1:nmnp) = NaN(nmnp, 1);
                            end
                        else
                            tb(1:nmiss) = NaN(nmiss, 1);
                        end
                    end
                    Mbeta = [g, seb, tb];
                    outa.model.hb = g;
                    outa.model.Mb = M;
                    outa.model.Y = Y;
                    outa.model.seb = seb;
                    outa.model.tb = tb;
                else
                    if mY == 1
                        Y = [];
                        flag = 0;
                    else
                        Y = Y(:, i+1:mY);
                        if rnamesr == 1
                            [nrg, mr] = size(rnamesrg);
                            rnamesrg = rnamesrg([1:i - 1, i + 1:nrg], :);
                        end
                        %       [beta,tv]=btval([],[y Y]);  yc=y-Y*beta;
                        xv = x0(pvar);
                        xf = x0(pfix);
                        if ~isempty(xv)
                            xv = arimaopt(fmarqdt, fid, x0, xv, xf, y, Y, parm, infm, 0);
                            x(pvar) = xv;
                            x(pfix) = xf;
                        else
                            x(pvar) = xv;
                            x(pfix) = xf;
                        end
                        if (p > 0)
                            outa.model.phi = x(1:p);
                        end
                        if (q > 0)
                            outa.model.th = x(p+1:p+q);
                        end
                        if (ps > 0)
                            pq = p + q;
                            outa.model.phis = x(pq+1:pq+ps);
                        end
                        if (qs > 0)
                            pqps = p + q + ps;
                            outa.model.ths = x(pqps+1:pqps+qs);
                        end
                        %         [F,e,g,M,A,P,matsis]=residual2(x,y,Y,s,dr,ds,p,ps,q,qs,1);
                        [F, e, g, M, A, P, matsis] = residual2x(x, y, Y, s, S, dr, ds, dS, p, ps, q, qs, qS);
                        Ss = e' * e;
                        Ff = F' * F;
                        ndrs = ny - nmiss + npmiss - dr - ds * s - dS * S; %residual sum of squares
                        nparm = length(xv) + nreg;
                        conp = Ss / (ndrs - nparm); %estimated sigma square
                        sconp = sqrt(conp);
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
                        seb = sqrt(diag(M)*conp);
                        tb = g ./ seb; %standard errors and t-values
                        %         if nmiss > 0
                        %          tb(1:nmiss)=NaN(nmiss,1);
                        %         end
                        if nmiss > 0
                            if npmiss > 0
                                if npmiss1 <= nmiss
                                    nmnp = nmiss - npmiss;
                                    tb(1:nmnp) = NaN(nmnp, 1);
                                end
                            else
                                tb(1:nmiss) = NaN(nmiss, 1);
                            end
                        end
                        Mbeta = [g, seb, tb];
                        outa.model.hb = g;
                        outa.model.Mb = M;
                        outa.model.Y = Y;
                        outa.model.seb = seb;
                        outa.model.tb = tb;
                    end
                end
            else
                flag = 0;
            end
        end
    end
    %end of backward elimination for transfer function
    [nn, mm] = size(nelinput);
    inputidx = zeros(nn, mm);
    inputidx = [];
    for i = 1:nn
        if nelinput(i) ~= 0
            inputidx = [inputidx; i];
        end
    end
    parm.inputidx = inputidx;
    outa.model.inputidx = inputidx;
end

%  rnamesrg,Mbeta,pause
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
    t = 1:ny + npr;
    for i = 1:nrout
        Youtg = [Youtg, Y(t, nreg+i) * g(nreg+i)]; %generate matrix of outlier effects
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
        SS = pinv(H/2) / (ny - dr - ds * s - nr - ninput - nreg);
        se = sqrt(abs(diag(SS)))';
    else
        %standard errors via the jacobian
        smvx = infm.mvx;
        infm.mvx = 0;
        %       [f,junk] = fstarima1(x,y,Y,parm,infm,xf);
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
if ~isempty(pvar)
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
            rnames = char(rnames, rnamesrg);
        else
            if isempty(rnamesrg)
                for i = 1:nreg, rnames = char(rnames, ['reg', num2str(i)]);
                end
            else
                rnames = char(rnames, rnamesrg);
            end
        end
        if nrout > 0, rnames = char(rnames, ornames);
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
    end
    %print residual diagnostics
    printres(fid, infr);
    % print summary
    prsummry(ii, ny, nreg0, fidr, dbname, parm, iout);
    fprintf(fid, '\nSecond Round of the Program:\n\n');
end

%transfer function identification and estimation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%second round for the program
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if isfield(ser, 'ninput')
    ninput = ser.ninput;
    if (ninput > 0) && (tfident == 1)
        [mn, nn] = size(nelinput);
        Yina = []; %firstelina=[]; %Yina will contain the significant input variables
        delaya = [];
        maa = [];
        ara = [];
        for i = 1:mn
            if nelinput(i) == 0 %an input is not significant
                ninput = ninput - 1;
            else
                Yina = [Yina, Yin(:, i)]; %firstelina=[firstelina; firstelin(i)];
            end
        end
        Yin = Yina;
        ser.Yin = Yin;
        [nYin, mYin] = size(Yin); %firstelin=firstelina;
        ser.ninput = ninput;
        ser.delay = delaya;
        ser.ma = maa;
        ser.ar = ara;
        if (npr > 0)
            for i = 1:ninput
                ser.modpred(i).pred = ser.modpred(inputidx(i)).pred;
            end
        end
    end
end

%update regression and input variables, taking the forecast horizon into
%account
if ninput > 0
    xm = x;
    nregm = nreg;
    prm = pr;
    pvarm = pvar;
    pfixm = pfix;
    %  nvinput=(nlagtf+1)*ninput;
    nvinput = sum(cnelinput);
    nreg = nreg - nvinput;
    parm.nreg = nreg;
    parm = tfivparm(Yin, parm);
    mY = mY - nvinput + npmiss;
    initreg = initreg + npmiss;
    if ~isempty(rnamesrgm0(1:npmiss, :))
        rnamesrg = char(rnamesrgm0(1:npmiss, :), rnamesrg);
    end
    nyinpr = nyin + npr;
    if mY > 0
        if (initreg > 0)
            Y = Yreg;
        else
            Y = [];
        end
        if (flagm == 1)
            ct = constant(nyinpr, 1, dr, ds, 0, 0, s);
            if mY > 0
                if nmiss > 0
                    Y = [Y(:, 1:nmiss), ct];
                else
                    Y(:, 1) = ct;
                end
            else
                Y = ct;
            end
            %     Y=[ct Y];
        end
        %set up rest of regression variables
        if (trad > 0)
            Ytr = trade(bg_year, bg_per, nyinpr, trad, 0, [], freq);
            Y = [Y, Ytr];
        end
        if (leapy > 0)
            Yle = genleap(bg_year, bg_per, nyinpr, freq);
            Y = [Y, Yle];
        end
        if (easte > 0)
            Yd = eastdate;
            Ye = east(bg_year, bg_per, nyinpr, dur, freq, Yd);
            Y = [Y, Ye];
        end
        if (nrout > 0)
            bb = 1;
            aa = [1, -delta];
            % vdelta is a vector containing the delta^i weights
            vdelta = poldiv(bb, aa, nyinpr);
            Yo = [];
            tip = iout.tip;
            nind = iout.nind;
            for i = 1:nrout
                tipo = tip(i, :);
                nindo = nind(i, :) + nlagtf;
                if strcmp(tipo, 'AO')
                    v1 = zeros(nyinpr, 1);
                    v1(nindo) = 1;
                    Yo = [Yo, v1];
                elseif strcmp(tipo, 'TC')
                    v2 = zeros(nyinpr, 1);
                    v2(nindo:nyinpr) = vdelta(1:nyinpr-nindo+1);
                    Yo = [Yo, v2];
                elseif strcmp(tipo, 'LS')
                    v3 = zeros(nyinpr, 1);
                    v3(nindo:nyinpr) = ones(nyinpr-nindo+1, 1);
                    Yo = [Yo, v3];
                end
            end
            Y = [Y, Yo];
        end
        if (rnamesr == 1)
            fpi = flagm + initreg;
            rnamesrgax = rnamesrg(1:fpi, :);
            [nrg, mrg] = size(rnamesrg);
            if (nrg > fpi + nvinput)
                if isempty(rnamesrgax)
                    rnamesrgax = char(rnamesrg(fpi+nvinput+1:end, :));
                else
                    rnamesrgax = char(rnamesrgax, rnamesrg(fpi+nvinput+1:end, :));
                end
            end
            rnamesrg = rnamesrgax;
        end
    else
        Y = [];
        rnamesr = 0;
        rnamesrg = [];
    end
    [nY, mY] = size(Y);
    
    %transfer function identification.
    nvinput = (nlagtf + 1) * ninput;
    tford = zeros(ninput, 3);
    if tfident > 0
        nlagtf1 = nlagtf + 1;
        %gm, extended vector of regression estimates if
        %backwd = 1. It contains the zeros.
        gm = zeros(nreg+nvinput, 1);
        gm(1:nreg) = g(1:nreg);
        if backwd == 1
            cont = 0;
            for i = 1:ninput
                %    i
                delay = nlagtf1 - cnelinput(i);
                tford(i, 1) = delay;
                cont = cont + cnelinput(i);
                gm(nreg+nlagtf1*(i - 1)+delay+1:nreg+nlagtf1*i) = ...
                    g(nreg+cont-cnelinput(i)+1:nreg+cont);
            end
        else
            gm = g;
            for i = 1:ninput
                delay = ser.delay(i);
                tford(i, 1) = delay;
            end
        end
        for jj = 1:ninput
            dbm = 0;
            dam = 0;
            mdis = 1.d10;
            %     for kk=jj+1:ninput
            %      tford(kk,2:3)=[0 2]; %as a first approximation, order two in the denominator
            %     end
            %     ord=nr;
            %     for kk=1:jj-1
            %      ord=ord+tford(kk,2)+1+tford(kk,3);
            %     end
            for db = 0:maxndtf %order of numerator
                for da = 0:maxddtf %order of denominator
                    %       tford(jj,2:3)=[db da];
                    %       pvar=pvarm; pfix=pfixm;
                    ng = length(gm(nreg+nlagtf1*(jj - 1)+tford(jj, 1)+1:nreg+nlagtf1*jj));
                    prs = db + da + 1;
                    %       if ng - prs - 1 > 0
                    mhr = max(da, db);
                    if (ng - mhr - 1 >= da)
                        [num, den] = shank(gm(nreg+nlagtf1*(jj - 1)+tford(jj, 1)+1:nreg+nlagtf1*jj), db, da);
                        %        [xx,pvar,pfix,parm]=xmparm(ninput,x0,xm,xf,pvar,pfix,nr,nlagtf,gm,tford,parm,nreg);
                        %        num=xx(ord+1:ord+1+db);
                        %        if da > 0
                        %         den=[1;xx(ord+2+db:ord+1+db+da)'];
                        %        else
                        %         den=1;
                        %        end
                        %      num
                        %      den
                        psi = poldiv(num, den, ng-1)';
                        piw = gm(nreg+(nlagtf + 1)*(jj - 1)+tford(jj, 1)+1:nreg+(nlagtf + 1)*jj);
                        %      dis = sum((psi-piw).^2).^0.5;
                        dng = double(ng);
                        %      dis = log(sum((psi-piw).^2)/dng) + double(prs)*log(dng)/dng;  %BIC
                        dis = log(sum((psi - piw).^2)/dng) + 2 * double(prs) / double(ng-prs-1); %corrected AIC
                    else
                        dis = 1.d10;
                        %         fprintf(1,'\n%s %3.0f\n','Too few weights to estimate transfer function model for input ',jj);
                        %         fprintf(1,'%s %2.0f %s %2.0f\n','with numerator degree',db,'and denominator degree',da);
                        %        return
                    end
                    if (dis < mdis)
                        mdis = dis;
                        dbm = db;
                        dam = da;
                    end
                end
            end
            %     dbm,dam
            tford(jj, 2:3) = [dbm, dam]; %tford contains the estimated orders
        end
        pr = prm;
    else
        for i = 1:ninput
            tford(i, 1) = ser.delay(i);
            tford(i, 2) = ser.ma(i);
            tford(i, 3) = ser.ar(i);
        end
    end
    %end of transfer function identification
    
    if lam == 0, y = log(yin);
    else y = yin;
    end %restore original series
    ny = length(y);
    parm.ny = ny;
    %regression (Y), input (Yin), outliers (Yo)
    parm = tfivparm(Yin, parm);
    pvar = pvarm;
    pfix = pfixm;
    [xx, pvar, pfix, parm] = xmparm(ninput, x0, xm, xf, pvar, pfix, nr, nlagtf, gm, tford, parm, nreg);
    clear aenames
    [nlag, aenames] = lagaena(parm); %names for Arima estimation printing
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %
    %search for outliers in the first nlagtf observations
    nrout0 = 0;
    if (out == 1) && (ninput > 0) %outlier detection
        nvinput = sum(cnelinput);
        nreg = nregm - nvinput;
        parm.nreg = nreg;
        xv = xx(pvar);
        xf = xx(pfix);
        %transfer function estimation before outlier detection
        f = 'fasttf';
        infm.f = f;
        if ~isempty(xv)
            xv = arimaopt(fmarqdt, fid, xx, xv, xf, y, Y, parm, infm, pr);
            xx(pvar) = xv;
        end
        sp1 = 1;
        sp2 = nlagtf; %search for outliers in the initial span
        ppsqqsqS = p + ps + q + qs + qS;
        x0 = xx(1:ppsqqsqS);
        parm0 = parm;
        parm0.ninput = 0;
        parm0.pvar = pvarm;
        parm0.pfix = pfixm;
        [nrout0, nind0, tip0, Yo0] = ...
            outlr(y, Y, parm0, iout, infm, npr, x0, sp1, sp2, fmarqdt, ols, a);
        iout.nrout = nrout + nrout0;
        iout.nind = [iout.nind; nind0];
        iout.tip = [iout.tip; tip0];
        iout.Yo = Yo0; %update iout structure
        Y = Yo0;
        [nY, mY] = size(Y);
    end
    %print outliers
    if nrout0 > 0
        for i = nrout + 1:nrout + nrout0
            ornames = strvcat(ornames, ['out', num2str(i)]);
        end
        iout.ornames = ornames;
        nrout = nrout + nrout0;
        if pr == 1
            fprintf(fid, '\n');
            trtout(fid, iout, ser); %print outliers
            fprintf(fid, '\n');
        end
    elseif (out == 1) && (pr == 1) && (ninput > 0)
        %    if omet == 1, met='Exact max. likelihood'; else met='Hannan-Rissanen'; end
        %    fprintf(fid,'\n%s %3.1f %s\n\n','No outliers detected in initial stretch (C = ',C,[', Method is '...
        %      met '):']);
        if nrout > 0
            iout.nind = iout.nind + nlagtf;
            fprintf(fid, '\n');
            trtout(fid, iout, ser); %print outliers
            fprintf(fid, '\n');
        end
    end
    %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %
    %transfer function estimation
    %
    %this line added on 14-3-2008
    x = zeros(size(xx));
    xv = xx(pvar);
    xf = xx(pfix);
    f = 'fasttf';
    infm.f = f;
    if pr == 1
        fprintf(fid, '\n');
    end
    if (nlestim == 1)
        xv = arimaopt(fmarqdt, fid, xx, xv, xf, y, Y, parm, infm, pr);
        x(pvar) = xv;
        x(pfix) = xf;
    else
        x(pvar) = xv;
        x(pfix) = xf;
    end
    %end of transfer function estimation
    %
    
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
                ps0 = ps;
                qs0 = qs;
                nr0 = p + ps + q + qs + qS;
                parm.ds = 1;
                parm.ps = 0;
                ps = 0;
                parm.qs = 1;
                qs = 1;
                nr = p + ps + q + qs + qS; %number of arma parameters
                pvar0 = parm.pvar;
                parm.pvar = 1:nr; %no fixed parameters when autmid = 1
                [s, S, dr, ds, dS, p, ps, q, qs, qS, ny, nreg, pfix, pvar] = imparm(parm);
                if flagm == 1
                    Y = [Y(:, 1:nmiss), Y(:, nmiss+2:end)];
                    [~, mY] = size(Y);
                    nregm = nregm - 1;
                    flagm = 0;
                    rnamesrg = [rnamesrg(1:nmiss, :); rnamesrg(nmiss+2:end, :)];
                    parm.flagm = 0;
                    nreg = nreg -1;
                    parm.nreg = nreg;
                end
                if pr == 1, prmod11x(fid, s, p, dr, q, ps, ds, qs, S, dS, qS, lam, flagm);
                end
                ct = constantx(ny, 1, dr, ds, dS, 0, 0, s, S); %generate a mean for the series
                if nreg > 0, YY = [ct, Y(1:ny, :)];
                else YY = ct;
                end
                est = 1;
                x0 = cinest(y, YY, parm, est, ols, a, 0, fid); %regression without inputs
                x0 = [x0, x(nr0+1:end)];
                pvar = [pvar, pvar0(nr0+1:end) + ps - ps0 + qs + qs0];
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
    if (dr == 0) && (p >= 1) && (fixdif == 0) && (autmid == 1)
        aa = roots([1, x0(1:p)]);
        [Rr,ii] = max(real(aa)); 
        alphamax = .499;
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
        if ((Rr > hm1s) && (abs(imag(aa(ii))) < hms)  || ((ds == 1)...
                && (Rr > hm1sr)) && (abs(imag(aa(ii))) < hmsr) )
            p0 = p;
            nr0 = p + ps + q + qs + qS;
            parm.dr = 1;
            parm.p = 0;
            p = 0;
            nr = p + ps + q + qs + qS; %number of arma parameters
            pvar0 = parm.pvar;
            parm.pvar = 1:nr; %no fixed parameters when autmid = 1
            [s, S, dr, ds, dS, p, ps, q, qs, qS, ny, nreg, pfix, pvar] = imparm(parm);
            if flagm == 1
                Y = [Y(:, 1:nmiss), Y(:, nmiss+2:end)];
                [~, mY] = size(Y);
                nregm = nregm - 1;
                flagm = 0;
                rnamesrg = [rnamesrg(1:nmiss, :); rnamesrg(nmiss+2:end, :)];
                parm.flagm = 0;
                nreg = nreg -1;
                parm.nreg = nreg;
            end
            if pr == 1, prmod11x(fid, s, p, dr, q, ps, ds, qs, S, dS, qS, lam, flagm);
            end
            ct = constantx(size(y, 1), 1, dr, ds, dS, 0, 0, s, S); %generate a mean for the series
            if nreg > 0, YY = [ct, Y(1:ny, :)];
            else YY = ct;
            end
            est = 1;
            x0 = cinest(y, YY, parm, est, ols, a, 0, fid); %regression without inputs
            x0 = [x0, x(nr0+1:end)];
            pvar = [pvar, pvar0(nr0+1:end) + p - p0];
            parm.pvar = pvar;
            clear aenames
            [nlag, aenames] = lagaena(parm); %names for Arima estimation printing
            xv = x0;
            xv = arimaopt(fmarqdt, fid, x0, xv, xf, y, Y, parm, infm, pr);
            x = xv;
        end
    end
    % end of check regular underdifference
    
    
%     % check regular overdifference
%     if (dr >= 1) && (q >= 1) && (fixdif == 0) && (autmid == 1)
%         pps = p + ps;
%         ppsq = pps + q;
%         thr = [1, x(pps+1:ppsq)];
%         aa = max(real(roots(thr)));
%         if aa > .98
%             parm.dr = parm.dr - 1;
%             if (p > 0) || (q > 1)
%                 parm.q = parm.q - 1;
%                 q = q - 1;
%             end
%             nr = p + ps + q + qs + qS; %number of arma parameters
%             pvar0 = parm.pvar;
%             parm.pvar = 1:nr; %no fixed parameters when autmid = 1
%             [s, S, dr, ds, dS, p, ps, q, qs, qS, ny, nreg, pfix, pvar] = imparm(parm);
%             if flagm == 0
%                 ct = constantx(size(y, 1), 1, dr, ds, dS, 0, 0, s, S); %generate a mean for the series
%                 if nreg > 0, YY = [ct, Y(1:ny, :)];
%                 else YY = ct;
%                 end
%             else
%                 YY = Y;
%             end
%             est = 1;
%             x0 = cinest(y, YY, parm, est, ols, a, 0, fid); %regression without inputs
%             x0 = [x0, x(nr+1:end)];
%             pvar = [pvar, pvar0(nr+1:end)];
%             parm.pvar = pvar;
%             clear aenames
%             [nlag, aenames] = lagaena(parm); %names for Arima estimation printing
%             xv = x0;
%             xv = arimaopt(fmarqdt, fid, x0, xv, xf, y, Y, parm, infm, pr);
%             x = xv;
%             if (flagm == 0)
%                 [F, e, beta, M] = residual2x(x, y, YY, s, S, dr, ds, dS, p, ps, q, qs, qS);
%                 nyd = ny - nmiss - dr - ds * s - dS * S;
%                 myd = length(beta);
%                 sg = e' * e / (nyd - myd);
%                 se = sqrt(diag(M)*sg);
%                 t = beta ./ se;
%                 %test whether mean is significant
%                 if (abs(t(1)) > 2.0) && (flagm == 0)
%                     if nY > 0
%                         ct = constantx(nY, 1, dr, ds, dS, 0, 0, s, S); %generate a mean for the series
%                     else
%                         ct = constantx(ny+npr, 1, dr, ds, dS, 0, 0, s, S); %generate a mean for the series
%                     end
%                     mY = size(Y, 2);
%                     if mY > 0
%                         if nmiss > 0
%                             Y = [Y(:, 1:nmiss), ct, Y(:, nmiss+1:mY)];
%                             if isempty(rnamesrg(1:nmiss, :))
%                                 rnamesrg = 'mean';
%                             else
%                                 rnamesrg = char(rnamesrg(1:nmiss, :), 'mean');
%                             end
%                             if nmiss < mY
%                                 rnamesrg = char(rnamesrg, rnamesrg(nmiss+1:end, :));
%                             end
%                         else
%                             Y = [ct, Y];
%                             if isempty(rnamesrg)
%                                 rnamesrg = 'mean';
%                             else
%                                 rnamesrg = char('mean', rnamesrg);
%                             end
%                         end
%                     else
%                         Y = ct;
%                         if isempty(rnamesrg)
%                             rnamesrg = 'mean';
%                         else
%                             rnamesrg = char('mean', rnamesrg);
%                         end
%                     end
%                     [nY, mY] = size(Y);
%                     rnamesr = 1;
%                     nreg = nreg + 1;
%                     parm.nreg = nreg;
%                     flagm = 1;
%                     parm.flagm = flagm;
%                 end
%             end
%             if pr == 1, prmod11x(fid, s, p, dr, q, ps, ds, qs, S, dS, qS, lam, flagm);
%             end
%         end
%     end
%     % end of check regular overdifference
    
    
    %
    % get residuals
    %
    nvinput = sum(cnelinput);
    nreg = nregm - nvinput;
    parm.nreg = nreg;
    xv = x(pvar);
    xf = x(pfix);
    %yci: the series of filtered inputs
    parmc = parm;
    parmc.npr = 1;
    [F, xv, e, gg, M, yci] = fasttf(xv, y, Y, parmc, infm, xf);
    %ycii: series corrected by filtered inputs in the original scale
    ycii = deltafil([y(1:dr+s*ds+S*dS); yci], dr, ds, 0, 0, s);
    ycii = deltafil(ycii, 0, dS, 0, 0, S);
    %   [F,e,gg,M,A,P,matsis]=residual2(x,ycii,Y,s,dr,ds,p,ps,q,qs,1);
    [F, e, gg, M, A, P, matsis] = residual2x(x, ycii, Y, s, S, dr, ds, dS, p, ps, q, qs, qS);
    
    nparm = length(xv) + nreg;
    Ss = e' * e;
    Ff = F' * F;
    ndrs = ny - nmiss - dr - ds * s - dS * S; %residual sum of squares
    conp = Ss / (ndrs - nparm); %estimated sigma square
    sconp = sqrt(conp);
    
    %check whether mean is significant
    if (flagm == 1) && (fixdif == 0)
        seb = sqrt(diag(M)*conp);
        tb = gg ./ seb; %standard errors and t-values
        if (abs(tb(nmiss+1)) < 2.0)
            if flagm == 1
                Y = [Y(:, 1:nmiss), Y(:, nmiss+2:end)];
                [~, mY] = size(Y);
                nregm = nregm - 1;
                flagm = 0;
                rnamesrg = [rnamesrg(1:nmiss, :); rnamesrg(nmiss+2:end, :)];
                parm.flagm = 0;
            end
            if pr == 1, prmod11x(fid, s, p, dr, q, ps, ds, qs, S, dS, qS, lam, flagm);
            end
            if ~isempty(xv)
                xv = arimaopt(fmarqdt, fid, x0, xv, xf, y, Y, parm, infm, pr);
                x(pvar) = xv;
            end
            nreg = nregm - nvinput;
            parm.nreg = nreg;
            [F, xv, e, gg, M, yci] = fasttf(xv, y, Y, parmc, infm, xf);
            %ycii: series corrected by filtered inputs in the original scale
            ycii = deltafil([y(1:dr+s*ds+S*dS); yci], dr, ds, 0, 0, s);
            ycii = deltafil(ycii, 0, dS, 0, 0, S);
            %   [F,e,gg,M,A,P,matsis]=residual2(x,ycii,Y,s,dr,ds,p,ps,q,qs,1);
            [F, e, gg, M, A, P, matsis] = residual2x(x, ycii, Y, s, S, dr, ds, dS, p, ps, q, qs, qS);
            nparm = length(xv) + nreg;
            Ss = e' * e;
            Ff = F' * F;
            ndrs = ny - nmiss - dr - ds * s - dS * S; %residual sum of squares
            conp = Ss / (ndrs - nparm); %estimated sigma square
            sconp = sqrt(conp);
        end
    end
    
    
    outa.tfmodel.matsis = matsis;
    outa.tfmodel.matsis.A = A;
    outa.tfmodel.matsis.P = P;
    outa.tfmodel.result.pvar = pvar;
    outa.tfmodel.result.pfix = pfix;
    outa.tfmodel.result.xv = xv;
    outa.tfmodel.result.xf = xf;
    outa.tfmodel.result.x = x;
    outa.tfmodel.result.fmarqdt = fmarqdt;
    outa.tfmodel.result.sigma2c = sconp;
    outa.tfmodel.result.Y = Y;
    outa.tfmodel.result.parm = parm;
    outa.tfmodel.result.infm = infm;
    
    %compute regression effects other than outliers
    Yrg = [];
    if nreg > 0
        nypnpr = ny + npr;
        t = 1:nypnpr;
        mreg = nreg - nmiss;
        Yrg = zeros(nypnpr, mreg);
        for i = nmiss + 1:nreg
            Yrg(:, i-nmiss) = Y(t, i) * gg(i); %generate matrix of regression effects
        end
    else
        mreg = 0;
    end
    
    
    %compute outlier effects
    Youtg = [];
    if nrout > 0
        t = 1:ny + npr;
        for i = 1:nrout
            Youtg = [Youtg, Y(t, nreg+i) * gg(nreg+i)]; %generate matrix of outlier effects
        end
    end
    outa.tfmodel.nreg = nreg;
    outa.tfmodel.Yrg = Yrg;
    outa.tfmodel.Youtg = Youtg;
    outa.tfmodel.yci = yci;
    
    
    %residual computations
    if ~isempty(gg) && olsres == 1
        e = matsis.olsres;
    end
    nv = length(xv);
    infr = rescomp(e, lag, nv, Ss, conp, sconp, Ff, ndrs, nreg);
    
    %computation of Pierce Qs value when there is one seasonality only
    if (s > 1) && (S == 0)
        ne = length(e);
        r = infr.r;
        orders = 1:floor(lag/s);
        nr = ps + qs;
        np = length(orders);
        Qs = zeros(np, 1);
        pvalQs = zeros(np, 1);
        dfQs = zeros(np, 1);
        for i = 1:np
            Qs(i) = sum((r(s:s:i*s).^2)./(ne - s:-s:ne - i * s)') * ne * (ne + 2);
            dfQs(i) = max(0, i-nr); %degrees of freedom
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
    
    if mY > 0
        seb = sqrt(diag(M)*conp);
        tb = gg ./ seb; %standard errors and t-values
        if nmiss > 0
            Interp = y(idxn) - gg(1:nmiss);
            tb(1:nmiss) = NaN(nmiss, 1);
            yinterp = y;
            yinterp(idxn) = Interp;
            sInterp = seb(1:nmiss);
            outa.tfmodel.yinterp = yinterp;
            outa.tfmodel.Interp = Interp;
            outa.tfmodel.sInterp = sInterp;
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
                outa.tfmodel.oInterp = oInterp;
                outa.tfmodel.osInterp = osInterp;
            end
        end
        Mbeta = [gg, seb, tb];
        outa.tfmodel.hb = gg;
        outa.tfmodel.Mb = M;
        outa.tfmodel.Y = Y;
        outa.tfmodel.seb = seb;
        outa.tfmodel.tb = tb;
    end
    g = gg;
    %
    
    
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
            SS = pinv(H/2) / (ny - dr - ds * s - nr - ninput - nreg);
            se = sqrt(abs(diag(SS)))';
        else
            %standard errors via the jacobian
            smvx = infm.mvx;
            infm.mvx = 0;
            %        [f,junk] = fstarima1(x,y,Y,parm,infm,xf);
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
    if ~isempty(pvar)
        tt(pvar) = x(pvar) ./ se;
    end
    
    z = [x', see', tt', nlag];
    
    
    if pr == 1
        %print transfer function estimation results
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
                for i = 1:nreg, rnames = char(rnames, ['reg', num2str(i)]);
                end
            end
            if nrout > 0, rnames = char(rnames, ornames);
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
                mprint([Interp, seb(1:nmiss)], in);
                fprintf(fid, '\n');
            end
        end
        %print residual diagnostics
        printres(fid, infr);
        % print summary
        prsummry(ii, ny, nreg0, fidr, dbname, parm, iout);
    end
    outa.tfmodel.tford = tford;
    cont = 0;
    if (p > 0)
        outa.tfmodel.phi = x(1:p);
        cont = cont + p;
    end
    if (q > 0)
        outa.tfmodel.th = x(p+1:p+q);
        cont = cont + q;
    end
    if (ps > 0)
        pq = p + q;
        outa.tfmodel.phis = x(pq+1:pq+ps);
        cont = cont + ps;
    end
    if (qs > 0)
        pqps = p + q + ps;
        outa.tfmodel.ths = x(pqps+1:pqps+qs);
        cont = cont + qs;
    end
    map = tford(:, 2);
    arp = tford(:, 3);
    for i = 1:ninput
        xomgi = [];
        xdeli = [];
        for j = 0:map(i)
            xomgi = [xomgi, x(cont+1)];
            cont = cont + 1;
        end
        xomg{i} = xomgi;
        for j = 1:arp(i)
            xdeli = [xdeli, x(cont+1)];
            cont = cont + 1;
        end
        xdel{i} = xdeli;
    end
    outa.tfmodel.omg = xomg;
    outa.tfmodel.del = xdel;
    outa.tfmodel.resinf = infr;
    outa.tfmodel.se = se;
    outa.tfmodel.tt = tt;
    outa.tfmodel.SS = SS;
end
%end of transfer function identification and estimation

%Forecasts
if npr > 0
    if ninput == 0
        %Arima
        Z = matsis.Z;
        T = matsis.T;
        H = matsis.H;
        [pry, spry] = predt(ny, npr, sconp, Y, Z, T, H, A, P, g, M);
        outa.model.npr = npr;
        outa.model.pry = pry;
        outa.model.spry = spry;
    else
        %Transfer function
        Z = matsis.Z;
        T = matsis.T;
        H = matsis.H;
        %the following forecasts are the forecasts of the disturbance term,
        %that is, the original series minus the filtered inputs
        [pry, spry] = predt(ny, npr, sconp, Y, Z, T, H, A, P, g, M);
        outa.tfmodel.npr = npr;
        outa.tfmodel.dpry = pry;
        outa.tfmodel.dspry = spry;
        parmx = parm; %copy parm structure into parmx
        %extend input variables with forecasts and compute sum of mses of
        %forecasts of filtered inputs
        sprx = zeros(npr, 1);
        delay = tford(:, 1);
        ma = tford(:, 2);
        ar = tford(:, 3);
        outa.tfmodel.Yin = Yin;
        for i = 1:ninput
            outa.tfmodel.modpred(i).pred = ser.modpred(i).pred;
            parmx.inputv(1:ny+npr, i) = [Yin(:, i); ser.modpred(i).pred];
            outa.tfmodel.modinput(i).mod = ser.modinput(i).mod;
            if ser.modinput(i).mod == 1
                phix = ser.modinput(i).phi;
                thetax = ser.modinput(i).theta;
                sigmai = ser.modinput(i).sigma2;
                outa.tfmodel.modinput(i).phi = ser.modinput(i).phi;
                outa.tfmodel.modinput(i).theta = ser.modinput(i).theta;
                outa.tfmodel.modinput(i).sigma2 = ser.modinput(i).sigma2;
                %convolution of input model with filter
                b = delay(i);
                org = p + ps + q + qs + qS;
                omega = x(org+1:org+ma(i)+1);
                if b > 0 %omegab = z^b*omega(z)
                    omegab = [zeros(1, b), omega];
                else
                    omegab = omega;
                end
                delta = 1;
                if ar(i) > 0
                    xtf = x(org+ma(i)+2:org+ma(i)+ar(i)+1);
                    delta = [delta, xtf];
                end
                numerator = fliplr(conv(fliplr(omegab), thetax));
                denominator = fliplr(conv(fliplr(delta), phix));
                psiw = poldiv(numerator, denominator, npr-1); %weights for the input mses
                for j = 1:npr
                    sprx(j) = sprx(j) + (psiw(1:j) * psiw(1:j)') * sigmai;
                end
            end
        end
        outa.tfmodel.y = y;
        yx = [y; zeros(npr, 1)]; %extend output with zeros
        parmx.npr = npr; %pass npr to fasttf
        %obtain series of filtered inputs extended with forecasts
        [fx, xvx, ex, gx, Mx, ydi] = fasttf(xv, yx, Y, parmx, infm, xf);
        ydii = deltafil([y(1:dr+s*ds+S*dS); ydi], dr, ds, 0, 0, s);
        ydii = deltafil(ydii, 0, dS, 0, 0, S);
        prx = -ydii(end-npr+1:end); %sum of forecasts of filtered inputs
        for i = 1:npr %output forecasts
            pry(i) = pry(i) + prx(i); %forecasts
            spry(i) = sqrt(spry(i)^2+sprx(i)); %mse of forecasts
        end
        outa.tfmodel.pry = pry;
        outa.tfmodel.spry = spry;
    end
    
    %obtain forecasts in the original scale using the log-normal
    %distribution
    opry = pry;
    ospry = spry;
    if lam == 0
        for i = 1:npr
            opry(i) = exp(pry(i)+(spry(i)^2)/double(2.));
            ospry(i) = exp(double(2.)*pry(i)+spry(i)^2) * (exp(spry(i)^2) - double(1.));
        end
        if ninput == 0
            outa.model.opry = opry;
            outa.model.ospry = ospry;
        else
            outa.tfmodel.opry = opry;
            outa.tfmodel.ospry = ospry;
        end
    end
    
    if pr == 1
        %print forecasts
        fprintf(fid, '\n');
        clear in
        in.fid = fid;
        if lam == 1
            in.cnames = char('Obs.      ', 'Forecasts ', 'Std. Error');
            in.fmt = char('%5.0f', '%12.4f', '%12.4f');
            mprint([(ny + 1:ny + npr)', pry, spry], in);
        else
            in.cnames = char('Obs.         ', 'Forecasts    ', 'Std. Error   ', ...
                'For. (levels)', 'Std. Error   ');
            in.fmt = char('%5.0f', '%12.4f', '%12.4f', '%19.4f', '%19.4f');
            mprint([(ny + 1:ny + npr)', pry, spry, opry, ospry], in);
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

yor = ser.yor;

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
            nyor = length(yor);
            yx = y(end-nyor+1:end);
            tabla = [tabla, yx];
        end
    end
    [pY, qY] = size(Y);
    if qY > 0
        [ntabla, mtabla] = size(tabla);
        tabla = [tabla, Y(1:ntabla, :) * g]; %series of regression effects
    end
    nname = ['results', filesep, dbname, '.num'];
    str = ['save ''', nname, ''' tabla -ascii -double'];
    eval(str);
end

if ~isempty(ycii)
    infr.ycii = ycii;
    outa.tfmodel.resinf = infr;
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
    dbname = strrep(dbname, '_', '\_');
    %Youtg: outlier effects, Yrg: regression effects.
    plotres(y, Y, g, yor, datei, cw, dbname, gft, nrout, Youtg, mreg, Yrg, infr, s, lam);
    if (nmiss > 0)
        figure
        vnames = char('Interpolated series');
        tsplot(yinterp, datei, vnames);
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
        vnames = char('Upper 95% band','Forecast','Lower 95% band');
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

% profile off
