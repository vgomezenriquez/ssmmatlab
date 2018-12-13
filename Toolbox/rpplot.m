function rpplot(r, p, sea, sep, c, fname)
%**************************************************************************
% This function creates plots of sample autocorrelations and sample partial
% correlations
%
%    INPUTS:
%        r : autocorrelations
%        p : partial autocorrelations
%      sea : standard errors of autocorrelations
%      sep : standard error of partial correlations
%        c : critical value of the standard normal distribution
%    fname : label used in the legend
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
%**************************************************************************

n = length(r);
t = 1:n;
subplot(2, 1, 1)
plot(t, zeros(1, n));
legend(fname);
hold on
for i = 1:n
    c1 = i * ones(1, n);
    c2 = 0:r(i) / (n - 1):r(i);
    plot(c1, c2);
    hold on
end
plot(t, c*sea, 'r');
hold on;
plot(t, -c*sea, 'r');
title('Sample autocor.')
hold off


subplot(2, 1, 2)
plot(t, zeros(1, n));
legend(fname);
hold on
for i = 1:n
    c1 = i * ones(1, n);
    c2 = 0:p(i) / (n - 1):p(i);
    plot(c1, c2);
    hold on
end
plot(t, c*sep, 'r');
hold on;
plot(t, -c*sep, 'r');
title('Sample partial autocor.')
hold off
