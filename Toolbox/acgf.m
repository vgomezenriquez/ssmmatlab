function c = acgf(phi, th, nc)
%
% This function computes the autocovariance function of an ARMA process
% phi(B)z_t=th(B)a_t
% The parameter nc is the number of desired autocovariances plus one,
% because the variance is included: g(0),g(1),...,g(nc-1)
%
%        Input parameters:
%      phi: a (p+1 x 1) array, where p is the degree of phi(z), containing
%           the ocefficients of phi(z) in degree descending order
%      th : a (q+1 x 1) array, where q is the degree of th(z), containing
%      nc : an integer, the number of desired covariances plus one
%           the ocefficients of th(z) in degree descending order
%        Output parameters:
%      c  : a (nc x 1) array containing the autocovariances in the order
%           g(0),g(1),...,g(nc-1)
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
p = length(phi) - 1;
q = length(th) - 1;
r = max(p, q);
if r == 0
    c = [(th / phi)^2, zeros(1, nc-1)];
    return
end
thi = fliplr(th);
if q < r
    thi = [thi, zeros(1, r-q)];
end
g = mpbf(thi, thi);
if p == 0
    nm = min(q+1+nc, length(g));
    c = g(q+1:nm);
    c = [c, zeros(1, nc-length(c)+1)];
else
    phii = fliplr(phi);
    if p < r
        phii = [phii, zeros(1, r-p)];
    end
    A = zeros(r+1);
    for i = 1:r + 1
        A(i, 1) = phii(i);
        A(i, r+1) = phii(r-i+2);
        if i > 1
            for k = 2:i
                A(i, k) = A(i, k) + phii(i-k+1);
                A(i, r+2-k) = A(i, r+2-k) + phii(r+1-i+k);
            end
        end
    end
    %  [Q,R]=qr(A)
    %  qg=Q'*g(1:r+1)';
    %  g=R\qg;
    g = A \ g(1:r+1)';
    g = g';
    c = poldiv(fliplr(g), phii(1:p+1), nc-1);
    c(1) = 2 * c(1);
end
