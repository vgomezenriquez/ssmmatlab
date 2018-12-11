function [y, Xm, nmiss, idnx] = chmarima(y)
%
%    this function checks whether there are missing values in the series
%
%        INPUTS:
%        y: an array containing the input series
%
%        OUTPUTS:
%        y: an array containing the input series with the missing values
%           replaced with tentative values
%       Xm: a regression matrix whose columns have zeros except for the
%           observation numbers of the missing values in which it has ones
%    nmiss: number of missing observations
%     idxn: index for missing values
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
idn = isnan(y);
nmiss = sum(idn);
if nmiss == 0
    Xm = [];
    idnx = [];
else
    n = length(y);
    nn = find(~idn);
    idnx = find(idn);
    yy = y(nn);
    nyy = length(yy);
    nl = min(10, nyy);
    Xm = zeros(n, nmiss);
    cont = 0;
    for i = 1:n
        if idn(i) == 1
            cont = cont + 1;
            Xm(i, cont) = 1;
            nm = find(nn > i);
            if (~isempty(nm))
                nmp = nm(1);
            else
                nmp = nn(nyy);
            end
            nr = max(1, min(nyy, nmp-nl)):min(nyy, nmp+nl-1);
            y(i) = mean(yy(nr));
        end
    end
end