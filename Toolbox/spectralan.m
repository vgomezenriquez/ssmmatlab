function spr = spectralan(y, per, win, corlag, graph, vnames, width, wina)
%
% Spectral analysis
%
% This programa computes the spectrum, coherence, phase delay, gain and
% cross correlations between a reference cycle and other cycles. If only
% one series is input, the periodogram and the autocorrelations are
% computed.
%
%       INPUTS :
%------------------
%            y : (ly x ny) matrix with the series;
%                if ny = 1, univariate spectral analysis and computation
%                of autocorrelations of y are performed,
%                if ny > 1, multivariate spectral analysis
%                and computation of cross-correlations are performed;
%                the program assumes that the first column contains
%                the reference series
%          per : frequency of the data (number of seasons)
%                if per < 0, it is set to 1 
%          win : window function used for (cross-)periodogram smoothing
%                0, no window is applied (nonsmoothed periodogram).
%                1, the Blackman-Tukey window
%                2, the Parzen window (default)
%                3, the Tukey-Hanning window
%                if win < 0, it is set to 2
%       corlag : number of leads and lags at which the
%                auto-/cross-correlations are computed; 
%                if corlag <= 0 or >= length(y), it is set to length(y)-1
%        graph : 0, do not produce graphs
%                1, produce graphs in the original scale
%                2, produce graphs in logarithms (default)
%                if graph < 0, it is set to 2
%       vnames : string cell array with names for the series; the program
%                assumes that their order coincides with the order in y;
%                default: refseries, series1, series2,...
%                if vnames <= 0, the program generates the names according
%                to the default
%       width : window width factor (1/3 by default)
%               if width <= 0, it is set to 1/3
%        wina : "a" parameter for Blackman-Tukey window (0.23 by default)
%               if wina <= 0, it is set to 0.23
%       OUTPUT :
%------------------
%       spr    : structure containing the following fields
%                always
%                .y      : input matrix y
%                .per    : input per
%                .names  : input names or names created by the program
%                .frq    : frequencies
%                depending on the size of y
%                fields for ny = 1:
%                .f      : (smoothed) periodogram of the series
%                .cr     : autocorrelations
%                fields for ny > 2,
%                .f      : matrix with columns containing the (smoothed)
%                          periodogram of the reference series and all the
%                          other series
%                .cr     : matrix containing the autocorrelations  of the
%                          reference series and the cross correlations
%                          between the reference series and all the other
%                          series
%                .co     : matrix with columns containing coherence
%                          between the ref. series and all the other series
%                .ga     : matrix with columns containing the gain
%                          between the ref. series and all the other series
%                .ph     : matrix with columns containing the phase delay
%                          between the ref. series and all the other series
%                .mc     : array containing the maximum coherence values
%                          between the reference series and all the other
%                          series
%                .phmc   : array containing the phase delays corresponding
%                          to the maximum coherence
%                .mcor   : array containing the maximum correlations
%                          between the reference series and all the other
%                          series
%                .imcor  : time indices corresponding to the maximum cross
%                          correlations
%                .mpa    : array containing the mean phase angle in
%                          radians between the reference series and all the
%                          other series
%                .mph    : array containing the mean phase delay between
%                          the reference series and all the other series
%
% Copyright (c) 22 September 2017 by Victor Gomez
% Ministerio de Hacienda, Direccion Gral. de Presupuestos,
% Subdireccion Gral. de Analisis y P.E.,
% Alberto Alcocer 2, 1-P, D-34, 28046, Madrid, SPAIN.
% Phone : +34-915835439
% E-mail: vgomez@sepg.minhap.es
%
% The author assumes no responsibility for errors or damage resulting from the use
% of this code. Usage of this code in applications and/or alterations of it should
% be referenced. This code may be redistributed if nothing has been added or
% removed and no money is charged. Positive or negative feedback would be appreciated.
%

[ly, ny] = size(y);
spr.y = y;
spr.per = per;
if nargin < 8
    wina = 0.23;
elseif wina <= 0
    wina = 0.23;
end
if nargin < 7
    width = 1./3.;
elseif width <= 0 
    width = 1./3.;
end
if win < 0
    win = 2;
end
if nargin < 6 || vnames <= 0
    dvnames = cell(1, ny);
    dvnames{1} = 'refseries';
    if ny > 1
        for i = 2:ny
            dvnames{i} = ['series', num2str(i-1)];
        end
    end
else
    dvnames = vnames;
end
spr.names = dvnames;
if nargin < 5
    graph = 2;
elseif graph < 0
    graph = 2; 
end
if nargin < 4
   corlag = ly - 1;
elseif corlag < 0 || corlag >= ly
   corlag = ly - 1;
end
if nargin < 3
  win = 2;
elseif win < 0 
  win = 2;
end
if nargin < 2
   per = 1;
elseif per <0
   per = 1;   
end

refc = y(:, 1); % reference series
rfname = dvnames{1}; % name of the reference series

% Compute frequencies
n = ly;
np = floor(n/2);
nfrq = np + 1;
frq = zeros(nfrq, 1);
for i = 0:np
    frq(i+1) = 2 * pi * i / n;
end
spr.frq = frq;
spr.f = zeros(nfrq, ny);
spr.cr = zeros(2*corlag, ny);
%---------------------------------------------
% Periodogram of reference series
[f, frq] = periodg(refc, win, width, wina);
spr.f(:, 1) = f;
% Autocorrelations of reference series
cr = zeros(2*corlag, 1);
ind = -corlag + 1:corlag;
for j = ind
    [cro, stdx, stdx] = croscor(refc, refc, j);
    cr(corlag+j) = cro;
end
spr.cr(:, 1) = cr;
% frequency band for which the results are displayed;
% it corresponds to business cycle periodicities (periods between 1.5 and 8
% years)
fb2 = (2 * pi) / (per * 1.5);
fb1 = (2 * pi) / (per * 8);
frqb = (frq >= fb1 & frq <= fb2); %frequencies in the cyclical band
%lines for the business cycle frequency band
cc = ones(100, 2);
cc(:, 1) = cc(:, 1) * (fb1);
cc(:, 2) = cc(:, 2) * (fb2);
%------------------------------------------------------------------
% Plot spectrum of reference series
if graph >= 1
    figure
    if graph == 1
        ll = (max(f) - min(f)) / 99;
        dd = min(f):ll:max(f);
        plot(frq, f, cc(:, 1), dd, cc(:, 2), dd)
        legend(['spectrum ', rfname])
    elseif graph == 2
        lf = log(max(f,1e-10));
        ll = (max(lf) - min(lf)) / 99;
        dd = min(lf):ll:max(lf);
        plot(frq, lf, cc(:, 1), dd, cc(:, 2), dd) 
        legend(['spectrum (in logs) ', rfname])
    end
    disp('press any key to continue')
    pause
    close all
end

if ny > 1
    spr.co = zeros(nfrq, ny-1);
    spr.ph = zeros(nfrq, ny-1);
    spr.ga = zeros(nfrq, ny-1);
    spr.mc = zeros(1, ny-1);
    spr.phmc = zeros(1, ny-1);
    spr.mcor = zeros(1, ny-1);
    spr.imcor = zeros(1, ny-1);
    spr.mpa = zeros(1, ny-1);
    spr.mph = zeros(1, ny-1);
    %     if out == 1
    %        fid=fopen(outf,'w');
    %        fid=1;
    %        fprintf(fid,'=====================================================\n');
    %        fprintf(fid,'\n               TIME-DOMAIN ANALYSIS\n');
    %     end
    for i = 2:ny % Multivariate analysis
        comc = y(:, i);
        rname = dvnames{i};
        %-----------------------------------------------------
        % Compute cross correlations
        cr = zeros(2*corlag, 1);
        ind = -corlag + 1:corlag;
        for j = ind
            [cro, stdx, stdx] = croscor(refc, comc, j);
            cr(corlag+j) = cro;
        end
        [mcor, j] = max(abs(cr)); %maximum positive correlation
        imcor = ind(j);
        mac = cr(corlag+imcor);
        %------------------------------------------------------
        % spectral analysis
        [co, ph, ga, fx, f, frq] = crosspan(refc, comc, win, width, wina);
        phcb = ph(frqb); %phase delays in the cyclical band
        mph = mean(phcb); %mean phase delay in the cyclical band
        [mc, jc] = max(co(frqb)); %maximum coherence in the cyclical band
        phmc = phcb(jc); %phase delay corresponding to maximum coherenc
        frc = frq(frqb);
        alpha = phcb .* frc; %phase angles in the cyclical band
        malpha = mean(exp(1i*alpha)); %mean of unit vectors
        mpa = angle(malpha); %mean phase angle in radians in the cycl. band
        % Save results
        spr.mcor(i-1) = mac;
        spr.imcor(i-1) = imcor;
        spr.cr(:, i) = cr;
        spr.f(:, i) = f;
        spr.co(:, i-1) = co;
        spr.ga(:, i-1) = ga;
        spr.ph(:, i-1) = ph;
        spr.mc(i-1) = mc;
        spr.phmc(i-1) = phmc;
        spr.mpa(i-1) = mpa;
        spr.mph(i-1) = mph;
        fid3 = 1;
        fprintf(fid3, '\n   %s', rname);
        fprintf(fid3, '\n   maximum correlation at t =%3i', imcor);
        fprintf(fid3, '\n   correlation =%7.3f', mac);
        fprintf(fid3, '\n   phase delay corresponding to maximum coherence =%7.3f', phmc);
        fprintf(fid3, '\n   maximum coherence =%7.3f', mc);
        fprintf(fid3, '\n   mean phase angle in radians (circular statistics) =%7.3f', mpa);
        fprintf(fid3, '\n   mean phase delay in the cyclical band =%7.3f\n', mph);
        disp('press any key to continue')
        pause
        %------------------------------------------------------------------
        % Plot spectrum
        if graph >= 1
            figure
            if graph == 1
                ll = (max(f) - min(f)) / 99;
                dd = min(f):ll:max(f);
                plot(frq, f, cc(:, 1), dd, cc(:, 2), dd)
                legend(['spectrum ', rfname])
            elseif graph == 2
                lf = log(max(f,1e-10));
                ll = (max(lf) - min(lf)) / 99;
                dd = min(lf):ll:max(lf);
                plot(frq, lf, cc(:, 1), dd, cc(:, 2), dd) 
                legend(['spectrum (in logs) ', rfname])
            end
            disp('press any key to continue')
            pause
            figure
            plot(frq, co, cc(:, 1), dd, cc(:, 2), dd)
            legend(['coherence ', rfname, '-', rname])
            disp('press any key to continue')
            pause
            ll = (max(ph) - min(ph)) / 99;
            dd = min(ph):ll:max(ph);
            figure
            plot(frq, ph, cc(:, 1), dd, cc(:, 2), dd)
            legend(['phase delay ', rfname, '-', rname])
            disp('press any key to continue')
            pause
            close all
        end
    end
else
    return
end
