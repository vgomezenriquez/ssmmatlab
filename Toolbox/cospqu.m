function [c, q] = cospqu(x, y, win, width, wina)
%
%        This function computes the (smoothed) cross-periodogram
%        and the quadrature spectrum
%
%     INPUTS:
%------------
%        x,y : series
%        win : window used for smoothing the periodogram;
%             = 0 : no smoothing is performed
%             = 1 : Blackman-Tukey window
%             = 2 : Parzen window
%             = 3 : Tukey-Hanning window
%               if win < 0, it is set to 2
%    width : window width factor (1/3 by default)
%             if width <= 0, it is set to 1/3
%     wina : "a" parameter for Blackman-Tukey window (0.23 by default)
%             if wina <= 0, it is set to 0.23
%
%
%    OUTPUTS:
%------------
%       c    : cospectrum
%       q    : quadrature spectrum
%
% Copyright (c) 21 July 2003 by Victor Gomez
% Ministerio de Hacienda, Direccion Gral. de Presupuestos,
% Subdireccion Gral. de Analisis y P.E.,
% Alberto Alcocer 2, 1-P, D-34, 28046, Madrid, SPAIN.
% Phone : +34-915835439
% E-mail: vgomez@sepg.minhac.es
%
% The author assumes no responsibility for errors or damage resulting from the use
% of this code. Usage of this code in applications and/or alterations of it should
% be referenced. This code may be redistributed if nothing has been added or
% removed and no money is charged. Positive or negative feedback would be appreciated.
%
if nargin < 5
    wina = 0.23;
elseif wina <= 0
    wina = 0.23;
end
if nargin < 4
    width = 1./3.;
elseif width <= 0 
    width = 1./3.;
end
if nargin < 3
    win = 2;
elseif win < 0
    win = 2;
end
n = length(x);
cxy0 = croscov(x, y, 0);
np = floor(n/2);
c = zeros(1, np+1);
q = zeros(1, np+1);
if win >= 1
    if win == 1
        [w, m] = blacktu(n, width, wina);
    elseif win == 2
        [w, m] = parzen(n, width);
    elseif win == 3
        [w, m] = tukhan(n, width);
    end
    cxy = zeros(1, m);
    cyx = zeros(1, m);
    for i = 1:m
        cxy(i) = croscov(y, x, i);
        cyx(i) = croscov(x, y, i);
    end
    for i = 1:np
        auxc = 0;
        auxq = 0;
        for j = 1:m
            auxc = auxc + w(j) * (cxy(j) + cyx(j)) * cos(2*pi*j*i/n);
            auxq = auxq + w(j) * (cxy(j) - cyx(j)) * sin(2*pi*j*i/n);
        end
        c(i+1) = (cxy0 + auxc) / (2 * pi);
        q(i+1) = auxq / (2 * pi);
    end
    auxc = (cxy + cyx) * w;
    c(1) = (cxy0 + auxc) / (2 * pi);
    q(1) = 0;
else
    cxy = zeros(1, np+1);
    cyx = zeros(1, np+1);
    for i = 1:np
        cxy(i) = croscov(y, x, i);
        cyx(i) = croscov(x, y, i);
    end
    for i = 1:np
        auxc = 0;
        auxq = 0;
        for j = 1:np
            auxc = auxc + (cxy(j) + cyx(j)) * cos(2*pi*j*i/n);
            auxq = auxq + (cxy(j) - cyx(j)) * sin(2*pi*j*i/n);
        end
        c(i+1) = (cxy0 + auxc) / (2 * pi);
        q(i+1) = auxq / (2 * pi);
    end
    auxc = sum(cxy+cyx);
    c(1) = (cxy0 + auxc) / (2 * pi);
    q(1) = 0;
end
q = -q; %change sign of quadrature spectrum
