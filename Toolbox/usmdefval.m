%script file containing structural model default values

if isfield(ser, 'Y') %Y contains the regression variables
    Y = ser.Y;
    if isfield(ser, 'rnamesrg') %names for regression variables
        rnamesrg = ser.rnamesrg;
    else
        rnamesrg = [];
    end
    if ~isempty(rnamesrg)
        rnamesr = 1;
    else
        rnamesr = 0;
    end
else
    Y = [];
    rnamesr = 0;
    rnamesrg = [];
end
ny = length(yor); %series length
[nY, mY] = size(Y);
nreg = mY; %nreg = number of initial regression parameters.
%The mean is not included
if isfield(ser, 'Ycomp') %Y contains the regression variables
    Ycompi = ser.Ycomp;
    mYcomp = size(Ycompi, 2);
    if mYcomp ~= mY
        disp('the number of elements in Ycomp should be equal to ')
        disp('the number of columns in Y')
        return
    end
    Ycomp = zeros(1, mYcomp);
    for i = 1:mYcomp
        if strcmp(Ycompi{i}, 'level')
            Ycomp(i) = 1;
        elseif strcmp(Ycompi{i}, 'slope')
            Ycomp(i) = 2;
        elseif strcmp(Ycompi{i}, 'seas')
            Ycomp(i) = 3;
        elseif strcmp(Ycompi{i}, 'cycle')
            Ycomp(i) = 4;
        elseif strcmp(Ycompi{i}, 'ar')
            Ycomp(i) = 5;
        elseif strcmp(Ycompi{i}, 'irreg')
            Ycomp(i) = 6;
        elseif isempty(Ycompi{i})
            Ycomp(i) = 0;
        else
            disp('the names in Ycomp are incorrect')
            return
        end
    end
else
    Ycomp = [];
end

if isfield(ser, 'W')
    W = ser.W;
else
    W = [];
end

if isfield(ser, 'npr')
    npr = ser.npr; %number of forecasts
else
    npr = 0;
end

ninput = 0;
inc = 0;


%generation of names for regression variables
if (rnamesr == 0) && (nreg > 0)
    rnamesrg = ['reg', num2str(1)];
    for i = 2:nreg
        rnamesrg = char(rnamesrg, ['reg', num2str(i)]);
    end
    rnamesr = 1;
end


yin = []; % we will store here the original series when there are lagged variables
nyin = 0;
Yreg = Y; % and here the original regression variables
initreg = nreg;
chb = 0; %chb=1, in fstlkhev OLS estimator and MSE are computed


if isfield(ser, 'comp')
    comp = ser.comp; %number of forecasts
else
    disp('structural models require a comp field in structure ser')
    return
end

if ~isfield(ser.comp, 'sqrtfil')
    ser.comp.sqrtfil = 0;
end

comp.freq = freq;
datei = cal(bg_year, bg_per, freq);
comp.datei = datei;


%series default parameters
s = freq; %Arima orders

if isfield(ser, 'lam')
    lam = ser.lam;
else
    lam = -1;
end

cw = 1.96; %coefficient for the confidence bands
ols = 0;
a = 2.0; %1.5;         %parameters for the Hannan-Rissanen method
% ols = 0 use the Levinson-Durbin algorithm
%     = 1 use OLS
% a = the exponent in log(n)^a for the length of the long AR
% tramo uses a=2.0


lag = min(36, max([floor(.15*ny), ...
    3 * s, 10]));%number of lags for autocorrelations

if isfield(ser, 'nlestim')
    nlestim = ser.nlestim; %nonlinear estimation
else
    nlestim = 1;
end

tolf = 1e-4; %tolf:  a parameter used for stopping

if isfield(ser, 'olsres') %olsres :  =1 use OLS residuals
    olsres = ser.olsres; %          =0 do not use OLS residuals
else
    olsres = 0;
end


fh = 1;
wd = 9;
nd = 2;
scale = 1; %table parameters:
%fh: flag for header and years
%wd: format width
%nd: number of decimal points
%scale: =1 scale data if necessary
%       =0 do not scale data
if isfield(ser, 'gft')
    gft = ser.gft;
else
    gft = 0;
end


if isfield(ser, 'pr') %printing results in external file
    pr = ser.pr; % pr = 1 print in external file
else %    = 0 do not print in external file
    pr = 1;
end
