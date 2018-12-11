function [P, ferror] = mlyapunov(F, Q, eigmod)
%
%
%        This function computes the solution to the LYAPUNOV equation
%        P=FPF' + Q
%
% Input parameters:
% F       = a nxn matrix
% Q       = a nxn symmetric matrix
% eigmod  = a number with which the modulus of the eigenvalues of F will be compared
% Output parameters:
% P       = the solution of the Lyapunov equation
%
% Copyright (c) 21 July 2003 by Victor Gomez
% Ministerio de Hacienda, Direccion Gral. de Presupuestos,
% Subdireccion Gral. de Analisis y P.E.,
% Alberto Alcocer 2, 1-P, D-34, 28046, Madrid, SPAIN.
% Phone : +34-915835439
% E-mail: VGomez@sepg.minhap.es
%
% If there is an eigenvalue with modulus greater than or equal to
% eigmod, the system is not solved. In this case, nonstat is equal to the
% number of nonstationary eigenvalues.

%
% The author assumes no responsibility for errors or damage resulting from the use
% of this code. Usage of this code in applications and/or alterations of it should
% be referenced. This code may be redistributed if nothing has been added or
% removed and no money is charged. Positive or negative feedback would be appreciated.
%

P = [];
ferror = 0;

if nargin < 3, eigmod = 1;
end
[U, F] = schur(F', 'complex');
F = F';
nonstat = size(find(abs(diag(F)) >= eigmod), 1);
if nonstat > 0
    disp('nonstable F in mlyapunov');
    ferror = 1;
    return
end
n = size(F, 1);
P = zeros(n);
Q = U' * Q * U;
for j = 1:n
    for i = j:n
        d = 1 - F(i, i) * F(j, j)';
        suma = F(i, 1:i) * P(1:i, 1:j) * F(j, 1:j)' + Q(i, j);
        if (abs(suma) < eps)
            P(i, j) = 0;
        elseif (abs(d) < eps)
            disp('Ill-conditioned P in mlyapunov');
            ferror = 2;
            return
        else
            P(i, j) = suma / d;
            P(j, i) = P(i, j)';
        end
    end
end
P = real(U*P*U');
