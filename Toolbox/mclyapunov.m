function [X, ferror] = mclyapunov(A, B, C)
%
%
%        This function computes the solution to the continuous time
%        LYAPUNOV equation
%        AX - XB = C
%
% Input parameters:
% A       = a nxn matrix
% B       = a mxm matrix
% C       = a nxm matrix
% Output parameters:
% X       = the solution of the continuous time Lyapunov equation
%
% The method is as follows. First, compute the Schur decomposition of
% A and B (complex version). Then, solve recursively a system of linear
% equations.
% QaTaQa'X-XQbTbQb'=C;  TaQa'XQb - Qa'XQbTb = Qa'CQb;
% Z=Qa'XQb; TaZ - ZTb = Qa'CQb;  X=real(QaZQb');
%
%
% Copyright (c) 21 July 2003 by Victor Gomez
% Ministerio de Hacienda, Direccion Gral. de Presupuestos,
% Subdireccion Gral. de Analisis y P.E.,
% Alberto Alcocer 2, 1-P, D-34, 28046, Madrid, SPAIN.
% Phone : +34-915835439
% E-mail: VGomez@sepg.minhap.es
%
%
% The author assumes no responsibility for errors or damage resulting from the use
% of this code. Usage of this code in applications and/or alterations of it should
% be referenced. This code may be redistributed if nothing has been added or
% removed and no money is charged. Positive or negative feedback would be appreciated.
%

ferror = 0;

[Qa, Ta] = schur(A, 'complex');
[Qb, Tb] = schur(B, 'complex');
[n, m] = size(C);
X = zeros(n, m);
Cc = Qa' * C * Qb;
%  Qa,Ta,Qb,Tb,Cc
for i = n:-1:1
    for j = 1:m
        sum = Cc(i, j);
        for k = i + 1:n
            sum = sum - Ta(i, k) * X(k, j);
        end
        for p = 1:j - 1
            sum = sum + Tb(p, j) * X(i, p);
        end
        d = Ta(i, i) - Tb(j, j);
        if (abs(sum) < eps)
            X(i, j) = 0;
        elseif (abs(d) < eps)
            disp('Some eigenvalues of Ta and Tb equal in mclyapunov')
            ferror = 3;
            return
        else
            X(i, j) = sum / d;
        end
    end
end
X = real(Qa*X*Qb');
