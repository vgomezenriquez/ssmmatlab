function [K, ierror] = nullref(R, Indx)
%
%---------------------------------------------------
% USAGE: [K,ierror] = nullref(R,Indx)
% where:    R = an n x m matrix in column echelon form obtained after applying
%               housref on A'. Thus, if [Q,r,Indx,ierror] = housref(A'),
%               then R=r'.
%         Ind = an index containing the l.i. rows (0) and the l.d.
%               rows (1) of R.
%---------------------------------------------------
% RETURNS:
%           K = a matrix in reversed row echelon form containing a basis of the
%               left null-space of R. Therefore, K*R=0;
%        ierror =1, dimension mismatch in R and Indx
%               =0, there are no errors on input
%---------------------------------------------------
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
[nr, mr] = size(R);
[ni, mi] = size(Indx);
K = [];
ierror = 0;
if (nr ~= mi)
    disp('R and Indx have wrong dimensions in nullref')
    ierror = 1;
    return
end
A = [];
p = sum(Indx);
li = 0;
ld = 0;
K = zeros(p, nr);
nli = min(mi-p, mr);
for i = 1:mi
    if Indx(i) == 0
        A = [A; R(i, 1:nli)];
        li = li + 1;
    else
        ld = ld + 1;
        k = -R(i, 1:li) / A(:, 1:li);
        K(ld, i) = 1;
        sm = 0;
        for j = 1:i - 1
            if Indx(j) == 0
                sm = sm + 1;
                K(ld, j) = k(sm);
            end
        end
    end
end
