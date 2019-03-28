function [w, m] = blacktu(n, m, a)
%
%        This function computes the weights for the
%        Blackman-Tukey window
%
%     INPUTS:
%         n : lentgh of the series; required input to compute window lag
%             size if m is not input to blacktu
%         m : window lag size
%         a : parameter in the window function
%
%    OUTPUTS:
%         w : weights of the Blackman-Tukey window
%         m : window lag size;
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
if nargin == 1
    m=floor(n^(.756));
    a = .25;
elseif nargin == 2
    a = .25;
end
w = zeros(m, 1);
for i = 1:m
    w(i) = 1 - 2 * a + 2 * a * cos(pi*(i / m));
end
