function [z, rx1] = armafil(y, omega, delta, b, inc)
%
% This function filters the series y_t  using the filter nu(z) =
% z^b*omega(z)/delta(z), where omega(z)=omega_0 + omega_1*z + omega_2*z^2 +
% ....+omega_q*z^q and delta(z) = 1 + delta_1*z + ... + delta_p*z^p.
%
% Input arguments:
% y: the input series
% omega: a (q+1) array containing the coefficients of omega(z) in ascending
%        order
% delta: a (p+1) array containing the coefficients of delta(z) in ascending
%        order
% b: the delay of the filter
% inc: = 1 initial conditions for the filter are estimated
%      = 0 initial conditions equal to zero
%
% Output arguments:
%   z: the output series
% rx1: the design matrix for the initial state x_1
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

p = length(delta) - 1;
q = length(omega) - 1;
qb = q + b;
r = max(p, qb);
% p,q,r
n = length(y);
z = zeros(n, 1);
rx1 = [];
if r == 0
    for i = 1:n
        z(i) = omega(1) * y(i);
    end
    return
end
%
% Setup matrices for the state space form
%
F = zeros(r);
F(1:r-1, 2:r) = eye(r-1);
if p > 0
    F(r, r-p+1:r) = -fliplr(delta(2:end));
end
if b > 0
    omega = [zeros(1, b), omega];
end
psi = poldiv(omega, delta, r);
G = psi(2:end)';
% omega, delta, b, F, G , inc
FF = sparse(F(r, :));
x = zeros(r, 1);
if inc == 1
    %this is to estimate x_1=beta by regression later. rx1 contains the
    %design matrix for beta.
    A = eye(r);
    rx1 = A(1, :);
    for i = 1:n - 1
        wk1 = FF * A;
        A(1:r-1, :) = A(2:r, :);
        A(r, :) = wk1;
        rx1 = [rx1; A(1, :)];
    end
end
% initial conditions by backward filtering starting with zeros
% ry=flipud(y);
% for i=1:n
%  wk1=FF*x;
%  x(1:r-1)=x(2:r);
%  x(r)=wk1;
%  x = x+G*ry(i);
% end
% x

%z is the filtered series with x_1=0;
for i = 1:n
    z(i) = x(1) + psi(1) * y(i);
    wk1 = FF * x;
    x(1:r-1) = x(2:r);
    x(r) = wk1;
    x = x + G * y(i);
end
