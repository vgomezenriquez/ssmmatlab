function ct = constantx(N, cons, dr, ds, dS, dc, xc, s, S)
%*************************************************************************
%       This function generates a constant variable
%       for an ARIMA model that can have as nonstationary
%       autoregressive part
%
%       D(B) = (1-B)^dr (1-B^s)^ds (1-B^S)^dS (1-2cos(xc)B+B^2)^dc,     (1)
%
%       where 0<=dr<=2, 0<=ds<=1, 0<=dS<=1, and 0<=dc<=5.
%
%       The model is
%
%       D(B)z_t = cons + ARMA                                           (2)
%
%       and the generated variable is
%
%       ct = cons/D(B)                                                  (3)
%
%  INPUTS:
%      N : length of the data vector
%   cons : constant in eq.(2)
%     dr : regular differences
%     ds : first seasonal differences
%     dS : second seasonal differences
%     dc : integer such that dc*2 is order of (1-2cos(xc)B+B^2)^dc
%     xc : frequency in the polynomial (1-2cos(xc)B+B^2);
%          0 < xc < pi => |cos(xc)| < 1, implying that the polynomial
%          has two complex conjugate roots with modulus one
%      s : first frequency of the data
%      S : second frequency of the data
%
% OUTPUTS:
%     ct : constant given by eq.(3)
%
% Copyright (c) 21 July 2015 by Victor Gomez
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

d = dr + ds * s + dS * S + dc * 2;
ct = ones(N, 1) * cons;
if d > 0
    ct(1:d) = zeros(d, 1);
end
if dr > 0
    for j = 1:dr
        for i = d + 1:N
            ct(i) = ct(i) + ct(i-1);
        end
    end
end
if ds > 0
    for i = d + 1:N
        ct(i) = ct(i) + ct(i-s);
    end
end
if dS > 0
    for i = d + 1:N
        ct(i) = ct(i) + ct(i-S);
    end
end
if dc > 0
    for j = 1:dc
        for i = d + 1:N
            ct(i) = ct(i) + 2 * cos(xc) * ct(i-1) - ct(i-2);
        end
    end
end
