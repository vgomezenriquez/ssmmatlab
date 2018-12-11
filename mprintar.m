function mprintar(ar, info, tit, strt)
% PURPOSE: print an (np,mp,kp) array in formatted form
%---------------------------------------------------
% USAGE:     mprint(ar,info)
% where: ar   = (np,mp,kp) array to be printed
%        tit         = a character string (tilte) for each (np,mp) matrix
%        strt       = an integer to start the counting of the (np,mp)
%                        matrices
%        info      = a structure containing printing options
%        info.begr = beginning row to print,    (default = 1)
%        info.endr = ending row to print,       (default = np)
%        info.begc = beginning column to print, (default = 1
%        info.endc = ending column to print,    (default = mp*kp)
%        info.cnames = an (mp*kpr x 1) string vector of names for columns (optional)
%                      e.g. info.cnames = strvcat('col1','col2');
%                      (default = no column headings)
%        info.rnames = an (np+1 x 1) string vector of names for rows (optional)
%                      e.g. info.rnames = strvcat('Rows','row1','row2');
%                      (default = no row labels)
%        info.fmt    = a format string, e.g., '%12.6f' or '%12d' (default = %10.4f)
%                      or an (mp*kp x 1) string containing formats
%        info.fid    = file-id for printing results to a file
%                      (defaults to the MATLAB command window)
%                      e.g. fid = fopen('file.out','w');
%        info.rflag  = 1 for row #'s printed, 0 for no row #'s (default = 0)
%        info.width  = # of columns before wrapping occurs (default = 80)
%---------------------------------------------------
% Copyright (c) 2 June 2016 by Victor Gomez
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
if (nargin < 4)
    parf = 0;
else
    parf = 1;
    contb = strt;
end
[np, mp, kp] = size(ar);
names = cell(1, mp*kp);
cont = 0;
for i = 1:kp
    cont = cont + 1;
    if parf == 1
        names{cont} = [tit, '(', num2str(contb), '):'];
        contb = contb + 1;
    else
        names{cont} = tit;
    end
    for j = 2:mp
        cont = cont + 1;
        names{cont} = ' ';
    end
end
info.cnames = char(names);
mprint(ar, info);
