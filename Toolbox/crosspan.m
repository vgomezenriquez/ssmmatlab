function [co, ph, ga, fx, fy, frq] = crosspan(x, y, win)
%
%        This function performs a cross spectral analysis.
%        See Granger's book for definitions.
%
%     INPUTS:
%------------
%         x : reference series
%         y : other series
%       win : window used for smoothing the periodogram;
%             = 0 : no smoothing is performed
%             = 1 : Blackman-Tukey window
%             = 2 : Parzen window
%             = 3 : Tukey-Hanning window
%    OUTPUTS:
%------------
%        co = coherence
%        ph = phase delay
%        ga = gain
%        fx = periodogram for x
%        fy = periodogram for y
%        frq= frequencies
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
fx = periodg(x, win);
fy = periodg(y, win);
[c, q] = cospqu(x, y, win);
[co, ph, ga] = cohepha(c, q, fx, fy);
n = length(x);
np = floor(n/2);
frq = zeros(np+1, 1);
for i = 0:np
    frq(i+1) = 2 * pi * i / n;
end
%
% phase delay function
%
for i = 1:np
    ph(i+1) = ph(i+1) / frq(i+1);
end
