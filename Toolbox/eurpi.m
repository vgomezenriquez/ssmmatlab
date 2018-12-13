function [nr, ferror] = eurpi(Pi, hm1)
%
%
%   This function estimates the number of unit roots in Pi
%
% Input arguments:
%         Pi : an m x m matrix
%         hm1: a number with which the absolute value of the roots is
%              compared
% Output arguments:
%         nr: an integer, the number of unit roots
%
% Copyright (c) 21 July 2003 by Victor Gomez
% Ministerio de Hacienda, Direccion Gral. de Presupuestos,
% Subdireccion Gral. de Analisis y P.E.,
% Alberto Alcocer 2, 1-P, D-34, 28046, Madrid, SPAIN.
% Phone : +34-915835439
% E-mail: VGomez@sepg.minhaP.es
%
% The author assumes no responsibility for errors or damage resulting from the use
% of this code. Usage of this code in applications and/or alterations of it should
% be referenced. This code may be redistributed if nothing has been added or
% removed and no money is charged. Positive or negative feedback would be appreciated.
%


nr = 0;
[s, ss] = size(Pi);
if (s ~= ss)
    ferror = 1;
    disp('Pi should be square in eurpi');
    return
end


[Q, T] = schur(Pi, 'complex'); % 'real' is the default
aa = max(abs(eig(T)));
if (aa > 1)
    %  disp('matrix Pir unstable')
    for kk = 1:s
        if abs(T(kk, kk)) > 1, T(kk, kk) = 1 / T(kk, kk);
        end
    end
    Pi = Q * T * Q';
    [Q, T] = schur(Pi, 'complex');
end
%number of unit roots determination
for i = 1:s
    %  abstii=abs(T(i,i));
    %   if (abstii > hm1) && (abs(abstii-1.) < 1-hm1)
    if (-real(T(i, i)) > hm1) && (abs(imag(T(i, i))) < 1 - hm1)
        nr = nr + 1;
    end
end
%eigenvalues are sorted in descending absolute value
% if s > 1
%  [Qt,Tt,ap] = SortSchur(Q,T,0+0i); %Qt*Tt*Qt' ,Pi
% else
%  Qt=Q; Tt=T;
% end