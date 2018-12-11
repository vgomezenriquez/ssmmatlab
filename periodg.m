function [f, frq] = periodg(x, win)
%        This function computes the (smoothed) periodogram
%
%     INPUTS:
%------------
%         x : series
%       win : window used for smoothing the periodogram;
%             = 0 : no smoothing is performed
%             = 1 : Blackman-Tukey window
%             = 2 : Parzen window
%             = 3 : Tukey-Hanning window
%    OUTPUTS:
%      f    : (smoothed) periodogram
%      frq  : array containing the frequencies
%
% Copyright (c) 21 July 2003 by Victor Gomez
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
n = length(x);
cxx0 = croscov(x, x, 0);
np = floor(n/2);
frq = zeros(np+1, 1);
for i = 0:np
    frq(i+1) = 2 * pi * i / n; %frequencies
end
f = zeros(1, np+1);
if win >= 1
    if win == 1
        [w, m] = blacktu(n);
    elseif win == 2
        [w, m] = parzen(n);
    elseif win == 3
        [w, m] = tukhan(n);
    end
    cxx = zeros(1, m);
    for i = 1:m
        cxx(i) = croscov(x, x, i);
    end
    for i = 1:np
        aux = 0;
        for j = 1:m
            aux = aux + w(j) * cxx(j) * cos(2*pi*j*i/n);
        end
        f(i+1) = cxx0 / (2 * pi) + aux / pi;
    end
    aux = cxx * w;
    f(1) = cxx0 / (2 * pi) + aux / pi;
else
    m = n - 1;
    cxx = zeros(1, m);
    for i = 1:m
        cxx(i) = croscov(x, x, i);
    end
    for i = 1:np
        aux = 0;
        for j = 1:m
            aux = aux + cxx(j) * cos(2*pi*j*i/n);
        end
        f(i+1) = cxx0 / (2 * pi) + aux / pi;
    end
    f(1) = (n / (2 * pi)) * (mean(x))^2;
end
