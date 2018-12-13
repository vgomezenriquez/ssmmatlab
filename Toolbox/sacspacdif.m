function [c0s, cvs, rs, fis, pcs] = sacspacdif(y, tname, dr, ds, freq, lag, cw, fplot)
%*************************************************************************
%  This function outputs cv0s, cvs, rs, fis and pcs described below.
%  In addition, it optionally plots the series before and after
%  differencing, as well as the sample autocorrelations and the sample
%  partial correlations of the differenced series.
%
%    INPUTS :
%         y : series
%     tname : name of the series
%        dr : number of regular differences
%        ds : number of seasonal differences
%      freq : frequency of the data
%       lag : number of lags up to which cvs,rs,fis,pcs of the differenced
%             y are computed
%        cw : critical value of the standard normal distribution
%     fplot : =1 plot series (default), =0 do not plot series
%
%   OUTPUTS :
%       c0s : sample variance of differenced y
%       cvs : sample autocovariances of differenced y
%        rs : sample autocorrelations of differenced y
%       fis : an (1 x lag) vector containing the AR(lag) polynomial
%       pcs : sample partial correlations of differenced y
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

if (nargin == 7)
    fplot = 1;
end
s = freq;
yor = y;
if dr > 0
    for i = 1:dr
        y = diferm(y, 1);
    end
end
if ds == 1
    y = diferm(y, s);
end
yc = y - mean(y) * ones(size(y)); %series must be mean-corrected
ic = 1;
[c0s, cvs, rs] = autcov(yc, lag, ic); %autocorrelations
[fis, pcs] = durlev(c0s, cvs); %partial autocorrelations
if (fplot == 1)
    % set(0, 'DefaultAxesFontSize', 10);
    subplot(2, 2, 1)
    plot(yor);
    title(tname);
    axis('tight');
    subplot(2, 2, 2)
    plot(y);
    dfname = ['Differenced series: (', num2str(dr), ',', num2str(ds), ')'];
    title(dfname);
    axis('tight');
    subplot(2, 2, 3)
    bar(rs);
    title('Sample autocorrelations');
    n = length(y);
    t = 1:lag;
    hold on;
    sea = zeros(lag, 1);
    sea(1) = 1 / sqrt(n);
    k = 0; %standard errors
    for i = 1:lag - 1
        k = k + rs(i) * rs(i);
        sea(i+1) = sqrt((1 + 2 * k)/n);
    end
    plot(t, sea*cw);
    hold on;
    plot(t, -sea*cw);
    hold off;
    subplot(2, 2, 4)
    bar(pcs);
    title('Sample partial autocorrelations');
    hold on
    sep = ones(lag, 1) * 1 / sqrt(n); %standard errors
    plot(t, sep*cw);
    hold on;
    plot(t, -sep*cw);
    hold off;
end