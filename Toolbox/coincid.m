function ncoin = coincid(nind, eind)
%
% this function detects the number of coincidences between the
% present and the past search for outliers
%
% Input arguments:
% nind: array containing the time index of the present outliers
% eind: array containing the time index of the previous outliers
% Output arguments:
% ncoin: number of coincidences
%
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
%
ncoin = 0;
n = length(nind);
m = length(eind);
if n > 0
    for i = 1:m
        for j = 1:n
            if eind(i) == nind(j)
                ncoin = ncoin + 1;
            end
        end
    end
end
