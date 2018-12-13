function prmod11x(fid, s, p, dr, q, ps, ds, qs, S, dS, qS, lam, flagm)
%
% this function prints in the file fid the ARIMA model specification of the
% form (p,dr,q) (ps,ds,qs)_s (0,dS,qS)_S, adding information as to whether
% logs have been taken and whether a mean has been included.
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

fprintf(fid, '\n%22s %1i %1s %1i %1s %1i %1s', 'Model changed to: (', p, ',', dr, ',', q, ')');
if s > 1
    fprintf(fid, '%1s %1i %1s %1i %1s %1i %2s %i', '(', ps, ',', ds, ',', qs, ')_', s);
end
if S > 1
    fprintf(fid, '%1s %1i %1s %1i %1s %1i %2s %i', '(', 0, ',', dS, ',', qS, ')_', S);
end
if lam == 0
    fprintf(fid, '%9s', ' in logs');
    if flagm == 1
        fprintf(fid, '%15s', 'and with mean ');
    end
else
    if flagm == 1
        fprintf(fid, '%15s', 'with mean ');
    end
end
fprintf(fid, '\n\n');
