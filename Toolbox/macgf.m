function [c, ierror] = macgf(phi, th, Sigma, nc)
%
% This function computes the autocovariance function of a VARMA process
% phi(B)z_t=th(B)a_t
% The parameter nc is the number of desired autocovariances plus one, because
% the variance is included: c(0),c(1),...,c(nc-1)
%
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
%
[np, mp, kp] = size(phi);
[nt, mt, kt] = size(th);
[ns, ms] = size(Sigma);
c = [];
ierror = 0;
if (np ~= nt | mp ~= mt | np ~= mp | ns ~= np)
    disp('wrong dimensions of phi, th or Sigma in macgf');
    ierror = 1;
    return
end
s = np;
c = zeros(s, s, nc);
p = kp - 1;
q = kt - 1;
for i = kp:-1:2
    if (phi(:, :, i) == 0)
        p = p - 1;
        phi = phi(:, :, 1:i-1);
    else
        break
    end
end
for i = kt:-1:2
    if (th(:, :, i) == 0)
        q = q - 1;
        th = th(:, :, 1:i-1);
    else
        break
    end
end
[np, mp, kp] = size(phi);
[nt, mt, kt] = size(th);
% phi
% th
% p
% q
r = max(p, q);


ieyep = 0;
if ~any(any(phi(:, :, 1)-eye(np)))
    ieyep = 1;
    A1 = [];
else
    A1 = phi(:, :, 1);
end

ths = th;
if any(any(Sigma-eye(ns)))
    for i = 1:kt
        ths(:, :, i) = th(:, :, i) * Sigma;
    end
end
% ths
% th
G = pmmulbf(ths, th); %generating function of MA part
[ng, mg, kg] = size(G);
Gm = G(:, :, q+1:end);
if (p == 0)
    nm = min(nc, q+1);
    c = zeros(np, np, nc);
    c(:, :, 1:nm) = Gm(:, :, 1:nm);
    for i = nm + 1:nc
        c(:, :, i) = zeros(np, np);
    end
    if (ieyep == 0)
        for i = 1:nm
            c(:, :, i) = A1 \ c(:, :, i) / A1';
        end
    end
    return
end

X = zeros(np, np, r+1);
Gme = zeros(np, np, p+1);
if (q < p)
    Gme(:, :, 1:q+1) = Gm;
    for i = q + 2:p + 1
        Gme(:, :, i) = zeros(np, np);
    end
elseif (q > p)
    %compute B_i, i=p+1,..,q and reduce the problem to degree p on both
    %sides of the equation
    %compute the necessary Psi weights
    Psi0 = th(:, :, 1);
    if (ieyep == 0), Psi0 = A1 \ Psi0;
    end %first weight
    Psi = th(:, :, 2) - phi(:, :, 2) * Psi0;
    if (ieyep == 0), Psi = A1 \ Psi;
    end %second weight
    phix = zeros(np, np, q+1);
    phix(:, :, 1:p+1) = phi;
    for i = p + 2:q + 1, phix(:, :, i) = zeros(np, np);
    end
    %  phix
    for i = 2:q - (p + 1) %rest of the weights
        A = [];
        for j = i:-1:2
            A = [A, phix(:, :, j)];
        end
        AA = -phix(:, :, i+1);
        if (ieyep == 0), AA = AA * Psi0;
        end
        Psix = th(:, :, i+1) + AA - A * Psi;
        if (ieyep == 0), Psix = A1 \ Psix;
        end
        Psi = [Psi; Psix];
    end
    %  Psi0
    %  Psi
    for i = p + 1:q %compute B_i, i=p+1,..,q
        X(:, :, i+1) = X(:, :, i+1) + ths(:, :, i+1) * Psi0';
        for j = i + 1:q
            X(:, :, i+1) = X(:, :, i+1) + ths(:, :, j+1) * Psi((j - i - 1)*np+1:(j - i)*np, :)';
        end
    end
    %  X
    Gme(:, :, 1) = Gm(:, :, 1); %change Gamma_i, i=1,...,p of the MA part
    for i = 2:p + 1
        summ = Gm(:, :, i);
        for j = p + 2:p + 2
            summ = summ - X(:, :, j) * phi(:, :, j-i+1)';
        end
        Gme(:, :, i) = summ;
    end
else
    Gme = Gm;
end
% Gme

%solve the symmetric polynomial matrix equation of degree p
[Y, ierror] = sympmeq(phi, Gme);
%insert the solution in X
for i = 1:p + 1
    X(:, :, i) = Y(:, :, i);
end
% X
% Gme

%compute covariances
c(:, :, 1) = 2 * X(:, :, 1);
if (ieyep == 0), c(:, :, 1) = A1 \ c(:, :, 1);
end
AA = c(:, :, 1);
c(:, :, 1) = .5 * (AA + AA'); %enforce symmetry
for i = 2:min(p+1, nc)
    summ = X(:, :, i) - .5 * phi(:, :, i) * AA;
    for j = 2:i - 1
        summ = summ - phi(:, :, i+1-j) * c(:, :, j);
    end
    c(:, :, i) = summ;
    if (ieyep == 0), c(:, :, i) = A1 \ c(:, :, i);
    end
end
if nc > p + 1
    for i = p + 2:min(r+1, nc)
        summ = X(:, :, i) - phi(:, :, p+1) * c(:, :, i-p);
        for j = 2:p
            summ = summ - phi(:, :, p+2-j) * c(:, :, i-p-1+j);
        end
        c(:, :, i) = summ;
        if (ieyep == 0), c(:, :, i) = A1 \ c(:, :, i);
        end
    end
    for i = r + 2:nc
        summ = -phi(:, :, p+1) * c(:, :, i-p);
        for j = 2:p
            summ = summ - phi(:, :, p+2-j) * c(:, :, i-p-1+j);
        end
        c(:, :, i) = summ;
        if (ieyep == 0), c(:, :, i) = A1 \ c(:, :, i);
        end
    end
end