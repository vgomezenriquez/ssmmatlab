function bic = cbic(x, yd, nd, s, p, ps, q, qs)
%
% this function computes the BIC criterion of an ARMA model
%
% Input arguments:
% x : array containing model parameters
% yd: vector containing the data
% nd: length of yd
% s:  number of seasons
% p:  AR order
% ps: order of the AR of order s
% q:  order of the regular MA
% qs: order of the MA of order s
%
% Output arguments:
% bic: the bic criterion
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

prs = p + s * ps;
qrs = q + s * qs;
pmps = p + ps;
qmqs = q + qs;
if pmps + qmqs == 0
    bic = log(yd'*yd/double(nd));
    return
end
ndat = max(prs, qrs);
%    [F,e,g,M]=residual2(x,yd(1:ndat),[],s,0,0,p,ps,q,qs,1);
[F, e] = residual3(x, yd(1:ndat), [], s, p, ps, q, qs);
res = zeros(nd, 1);
res(1:ndat) = e;
for i = 1:nd - ndat
    sum = yd(i+ndat);
    for k = 1:p
        sum = sum + x(k) * yd(i+ndat-k);
    end
    for k = 1:ps
        sum = sum + x(p+k) * yd(i+ndat-k*s);
        for ik = 1:p
            sum = sum + x(p+k) * x(ik) * yd(i+ndat-k*s-ik);
        end
    end
    for k = 1:q
        sum = sum - x(pmps+k) * res(i+ndat-k);
    end
    for k = 1:qs
        sum = sum - x(pmps+q+k) * res(i+ndat-s*k);
        for ik = 1:q
            sum = sum - x(pmps+q+k) * x(pmps+ik) * res(i+ndat-k*s-ik);
        end
    end
    res(i+ndat) = sum;
end
dn = double(nd);
dln = log(dn);
var = res' * res / dn;
bic = log(var) + double(pmps+qmqs) * dln / dn;
