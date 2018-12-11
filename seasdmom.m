function Y = seasdmom(N, datei)
%
%
%       This function generates a matrix of seasonal dummies orothogonal to
%       the mean
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
%

ip = datei.beg_per;
s = datei.freq;
Y = repmat(-1/(s - 1), N, s-1);
for i = 1:N
    ind = mod(ip, s);
    if ind > 0
        Y(i, ind) = 1;
    end
    ip = ip + 1;
end