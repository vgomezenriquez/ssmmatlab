%Script file for a univariate model with a complex seasonal pattern.
%The series is weekly US gasoline data in thousands of barrels per day,
%from February 1991 to July 2005. The data was used in De Livera, Hyndman
%and Snyder (2001), "Forecasting Time Series With Complex Seasonal Patterns
%Using Exponential Smoothing", Journal of the American Statistical
%Association, 1513-1527.

clear

data = load(fullfile('data', 'gasoline.dat'));
yor = data;
npr = 261; %number of forecasts
lam = 1; %no logs are taken
nreg = 0;
tname = 'gasoline';
fname = fullfile('results', 'gasoline.txt');

y = yor(1:end-npr, :); %transform series
Y = [];
comp.level = [1, 0.1, NaN];
comp.slope = [1, 0., 0];
comp.seasp{1} = [365.25 / 7, 20, 0., 0];
comp.ar = [1, .1, NaN];
comp.arp = [-.1; NaN];
% comp.irreg=[1 .1  NaN];
comp.sqrtfil = 1;


if npr > 0
    mpr = npr;
    npr = 0;
else
    mpr = 0;
end


[strm, ferror] = suusmm(comp, y, Y, npr);
if ferror > 0
    return
end

%estimate model
%
[resultm, strm] = usmestimm(y, strm);

disp(' ');
disp('******************** Results from estimation ********************');
disp(' ');
resultm.tv = resultm.tv';
mprintr(resultm)
disp(' ')
fprintf(1, '%s %9.4f\n', 'Concentrated parameter:', sqrt(resultm.sigma2c));
disp('press any key to continue')
pause


%create estimated model
xvf = resultm.xvf;
xf = resultm.xf;
[X, Z, G, W, T, H, ins, ii, ferror] = pr2usmm(xvf, xf, strm);

disp(' ');
disp('***** Estimated Model  *****');
disp(' ');
disp('Matrix T:')
disp(T)
disp('Matrix Z:')
disp(Z)
disp('Matrix G:')
disp(G)
disp('Matrix H:')
disp(H)
disp('More estimation and diagnostic details are in file "gasoline.txt"')
disp('in the subdirectory "results"')
disp('press any key to continue')
pause

%residual diagnostics
e = resultm.e;
F = resultm.F;
Ss = e' * e;
Ff = F' * F;
ne = length(e); %residual sum of squares
Pevf = resultm.Pevf; %prediction error variance
% disp('standard error (finite sample)')
SPevf = resultm.SPevf;
ny = length(y);
pvar = strm.pvar;
nr = length(pvar);
X = strm.X;
[junk, nbeta] = size(X);
ndrs = ne + nbeta;
freq = 1;
lagl = min(36, max([floor(.2*ny), 3 * freq, 10]));
infr = rescomp(e, lagl, nr, Ss, Pevf, SPevf, Ff, ndrs, nbeta);

%plot residual diagnostics
plotres([], [], [], [], [], 1.96, 'residuals', 1, 0, [], 0, [], infr, 1, 1);
close all

%file for output
fid = fopen(fname, 'w');
%print estimation results
%the date is immaterial here, only the number of columns (ncol)
bg_year = 1900;
bg_per = 1;
ncols = 10;
datei = cal(bg_year, bg_per, ncols);
printusmerm(fid, datei, tname, yor, y, ny, lam, strm, resultm, nreg, nbeta);

%print residual diagnostics
printres(fid, infr);
%close external file
if fid ~= 1
    fclose(fid);
end


%compute forecasts
if mpr > 0
    %hb, Mb, A and P are in structure result. Here, hb is the vector of
    %regression estimates and Mb is the matrix of mse. A is the
    %estimated state vector, x_{t|t-1}, obtained with the Kalman filter at the
    %end of the sample and P is the matrix of standard errors.
    hb = resultm.h;
    Mb = resultm.M;
    A = resultm.A;
    P = resultm.P;
    
    npr = mpr;
    [strm, ferror] = suusmm(comp, y, Y, npr);
    Xp = strm.X;
    Wp = strm.W;
    if ~isempty(Xp)
        Xp = Xp(end-npr+1:end, :);
    end
    if ~isempty(Wp)
        Wp = Wp(end-npr+1:end, :);
    end
    cw = 1.96;
    m = 1; %number of series
    [pry, mypr, alpr, malpr] = ssmpred(npr, m, A, P, Xp, Z, G, Wp, T, H, hb, Mb);
    spry = zeros(m, npr);
    sconp = sqrt(resultm.sigma2c);
    for i = 1:npr
        spry(:, i) = sqrt(diag(mypr(:, :, i))) * sconp;
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
    end
    %plot forecasts
    out.pry = pry;
    out.spry = spry;
    out.opry = opry;
    out.ospry = ospry;
    out.y = y;
    out.yor = yor;
    out.ny = ny;
    out.npr = npr;
    out.cw = cw;
    out.tname = tname;
    out.lam = lam;
    out.s = freq;
    pfctsusm(out);
end
