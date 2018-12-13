function arimasplot(Y, dr, dru, s, dsu, lag, rp, pcp, cw)
%**************************************************************************
% Auxiliary function called in arimasimul_d.m to plot the different series
%
%  INPUTS:
%        y : array containing the simulated ARIMA series
%       dr : number of regular differences
%      dru : number of regular differences entered by the user
%        s : seasonal frequency
%      dsu : number of seasonal differences entered by the user
%      lag : number of lags for the autocorrelations and partial autocor.
%       rp : theoretical autocorrelations
%      pcp : theoretical partial autocorrelations
%       cw : confidence interval parameter
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

y = Y;
if dr > 0
    for i = 1:dru
        y = diferm(y, 1);
    end
end
if dsu == 1
    y = diferm(y, s);
end
% set(0, 'DefaultAxesFontSize', 10);
subplot(3, 2, 1)
plot(Y);
title('Simulated series');
axis('tight');
subplot(3, 2, 2)
plot(y);
dfname = ['Differenced series: (', num2str(dru), ',', num2str(dsu), ')'];
title(dfname);
axis('tight');
yc = y - mean(y) * ones(size(y)); %series must be mean-corrected
ic = 1;
[c0s, cvs, rs] = autcov(yc, lag, ic); %autocorrelations
[fis, pcs] = durlev(c0s, cvs); %partial autocorrelations
subplot(3, 2, 3)
bar(rp);
title('Popul. autocorrelations');
subplot(3, 2, 4)
bar(pcp);
title('Popul. partial autocorrelations');
subplot(3, 2, 5)
bar(rs);
title('Sample autocorrelations');
n = length(y);
t = 1:lag;
hold on; sea=zeros(lag,1); sea(1)=1/sqrt(n); sm=0; %standard errors
for i = 1:lag - 1
    sm = sm + rs(i) * rs(i);
    sea(i+1) = sqrt((1 + 2 * sm)/n);
end
plot(t, sea*cw);
hold on;
plot(t, -sea*cw);
hold off;
subplot(3, 2, 6)
bar(pcs);
title('Sample partial autocorrelations');
hold on
sep = ones(lag, 1) * 1 / sqrt(n); %standard errors
plot(t, sep*cw);
hold on;
plot(t, -sep*cw);
hold off;