function ct = deltafil(x, dr, ds, dc, xc, s)
%**************************************************************************
%       This function generates a filtered variable
%       for an ARIMA model that can have as filter
%
%       D(B) = (1-B)^dr (1-B^s)^ds (1-2cos(xc)B+B^2)^dc,
%
%       where 0<=dr<=2, 0<=ds<=1, and 0<=dc<=5.
%
%       The filtered variable, ct, satisfies
%
%       D(B)ct_t = x_t
%
%  INPUTS:
%      x : input series
%     dr : regular differences
%     ds : seasonal differences
%     dc : integer such that dc*2 is order of (1-2cos(xc)B+B^2)^dc
%     xc : frequency in the polynomial (1-2cos(xc)B+B^2);
%          0 < xc < pi => |cos(xc)| < 1, implying that the polynomial
%          has two complex conjugate roots with modulus one
%      s : frequency of the data
%
% OUTPUTS:
%     ct : filtered series
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


drc = dr + 2 * dc;
dss = ds * s;
d = drc + dss;
[N, junk] = size(x);
ct = x;
pol = 1;
if dr > 0
    for j = 1:dr
        pol = conv(pol, [-1, 1]);
    end
end
if ds > 0
    pol = conv(pol, [-1, zeros(1, s-1), 1]);
end
if dc > 0
    for j = 1:dc
        pol = conv(pol, [1, -2 * cos(xc), 1]);
    end
end
pol = fliplr(pol);
m = length(pol);
for i = d + 1:N
    sum = ct(i);
    for j = 1:m - 1
        sum = sum - pol(j+1) * ct(i-j);
    end
    ct(i) = sum;
end