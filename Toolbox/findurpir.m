function [nr, ferror] = findurpir(Pi, Th, hm1, can)
%
%
%   This function checks whether there are unit roots in Pi
%
% Input arguments:
%                  Pi: an m x m matrix
%                  Th: an m x m matrix
%                 hm1: a number with which the absolute value of the roots
%                      is compared
%                 can: a small number to handle cancellation
%  Output: nr= an integer, number of unit roots
%         ferror= a flag for errors
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


nr = 0;
ferror = 0;
[s, ss] = size(Pi);
if (s ~= ss)
    ferror = 1;
    disp('Pi should be square in findurpir');
    return
end

if ~isempty(Th)
    [nt, mt] = size(Th);
    if (nt ~= mt)
        ferror = 1;
        disp('Th shoulb be square in findurpir');
        return
    end
end

[Q, T] = schur(Pi); % 'real' is the default
aa = max(abs(eig(T)));
if (aa > 1)
    %  disp('matrix Pi unstable')
    %  [Q,T]=rsf2csf(Q,T);
    for kk = 1:s
        if abs(T(kk, kk)) > 1, T(kk, kk) = 1 / T(kk, kk);
        end
    end
    %  Pi=real((Q*T)/Q);
    Pi = Q * T * Q';
    [Q, T] = schur(Pi);
end

if (s > 1)
    %eigenvalues are sorted in descending absolute value
    [Qt, Tt, ap] = SortSchur(Q, T, 0+0i);
else
    Tt = T;
    Qt = Q;
end
Qs = Qt;
%extract eigenvalues
if isOctave
    E = sort(eig(Tt));
else
    E = ordeig(Tt);
end


if ~isempty(Th)
    [Qt, Tt] = schur(Th);
    aa = max(abs(eig(Tt)));
    if (aa > 1)
        %  disp('matrix Th unstable')
        %  [Qt,Tt]=rsf2csf(Qt,Tt);
        for kk = 1:nt
            if abs(Tt(kk, kk)) > 1, Tt(kk, kk) = 1 / Tt(kk, kk);
            end
        end
        %   Th=real((Qt*Tt)/Qt);
        Th = Qt * Tt * Qt';
        [Qt, Tt] = schur(Th);
    end
    if (nt > 1)
        %eigenvalues are sorted in descending absolute value
        [Qtt, Ttt, apt] = SortSchur(Qt, Tt, 0+0i);
    else
        Ttt = Tt; %Qtt=Qt;
    end
    %extract eigenvalues
    if isOctave
        Et = sort(eig(Ttt));
    else
        Et = ordeig(Ttt);
    end
    if (-real(E(1)) > hm1) && (abs(E(1)-Et(1)) > can) ...
            && (abs(imag(E(1))) < 1 - hm1)
        nr = nr + 1;
    end
else
    if (-real(E(1)) > hm1) && (abs(imag(E(1))) < 1 - hm1)
        nr = nr + 1;
    end
end
