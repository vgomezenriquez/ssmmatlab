function [meind, metip] = sinfelo(eind, etip, nind, ntip, imin)
%
% this function stores information on the eliminated outlier
%
% Input arguments:
% eind: array containing the time index of the previously eliminated outliers
% etip: array containing the type of the previously eliminated outliers
% Output arguments:
% mind: array incorporating the new time index
% mtip: array incorporating the new type
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

meind = eind;
metip = etip;
neind = length(meind);
rep = 0;
if neind > 0
    for i = 1:neind
        if nind(imin) == meind(i)
            rep = 1;
        end
    end
else
    rep = 0;
end
if rep == 0
    meind = [meind; nind(imin)];
    metip = [metip; ntip(imin)];
end
