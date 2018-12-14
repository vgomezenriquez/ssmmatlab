function [x, laphi, x2] = hanris(yc, s, S, p, ps, q, qs, qS, ols, a)
%
% this function applies the Hannan-Rissanen method to obtain
% ARMA estimates
%
% Input arguments:
% yd   : a vector containing the series
% s    : seasonality
% S    : second seasonality
% p    : degree of AR polynomial
% ps   : degree of AR seasonal polynomial
% q    : degree of MA polynomial
% qs   : degree of MA seasonal polynomial
% qS   : degree of MA second seasonal polynomial
% ols  : = 1, perform OLS, = 0, use the Durbin Levinson algorithm
% a    : an integer, the degree of the AR approximation in the first step
%        of the Hanna-Rissanen method.
%
% Output arguments:
% x       : second or third step estimates (if process has an MA part and
%           the second step model is invertible)
% laphi   : a vector containing the estimates of the AR approximation in
%           first step of the Hanna-Rissanen method
% x2      : second step estimates
%
%
% Copyright (c) 21 July 2014 by Victor Gomez
% Ministerio de Hacienda y A.P., Direccion Gral. de Presupuestos,
% Subdireccion Gral. de Analisis y P.E.,
% Alberto Alcocer 2, 1-P, D-34, 28046, Madrid, SPAIN.
% Phone : +34-915835439
% E-mail: vgomez@sepg.minhap.es
%
% The author assumes no responsibility for errors or damage resulting from the use
% of this code. Usage of this code in applications and/or alterations of it should
% be referenced. This code may be redistributed if nothing has been added or
% removed and no money is charged. Positive or negative feedback would be appreciated.
%
laphi = [];
[n, m] = size(yc);
qt = q + s * qs + S * qS;
pt = p + s * ps;
if qt > 0
    % compute lag length
    N = floor(max(log(n)^a, max(pt, 2*qt)));
    N = round(min(N, n-n/4));
    % create lag matrix
    Y = [];
    for i = 1:N
        Y = [Y, yc(N-i+1:n-i)];
    end
    if ols == 1
        % perform OLS to estimate residuals
        phi = bols(yc(N+1:n), Y);
    else
        % use the Durbin Levinson algorithm
        [c0, cv, r] = autcov(yc, N, 0);
        [phi, pc] = durlev(c0, cv);
    end
    laphi = phi';
    a = zeros(n, 1);
    a(1) = yc(1);
    for i = 2:N
        a(i) = yc(i) - yc(i-1:-1:1)' * phi(1:i-1);
    end
    a(N+1:n) = yc(N+1:n) - Y * phi;
end
% perform OLS to estimate ARMA parameters
m = round(max(pt+1, qt+1));
% first regression
Y = [];
for i = 1:p
    Y = [Y, -yc(m-i:n-i)];
end
if ps > 0
    for i = s:s + p
        Y = [Y, -yc(m-i:n-i)];
    end
end
if qt > 0
    for i = 1:q
        Y = [Y, a(m-i:n-i)];
    end
    if qs > 0
        for i = s:s + q
            Y = [Y, a(m-i:n-i)];
        end
    end
    if qS > 0
        for i = S:S + q
            Y = [Y, a(m-i:n-i)];
        end
    end
    if qs > 0 && qS > 0
        for i = s + S:s + S + q
            Y = [Y, a(m-i:n-i)];
        end
    end
end
phi = bols(yc(m:n), Y);

psp = p + ps + ps * p;
qsq = q + qs + qs * q;
qSq = qS + qS * q;

%the following lines added 13-06-2008
%if the model is nonstationary or noninvertible, we do
%not go into the second step.
% if the model is noninvertible, the recursions to obtain the derivatives
% can explode.

x = zeros(1, p+ps+q+qs+qS);
if p > 0
    x(1:p) = phi(1:p);
end
if ps > 0
    x(p+1:p+ps) = phi(p+1:p+ps);
end
if q > 0
    x(p+ps+1:p+ps+q) = phi(psp+1:psp+q);
end
if qs > 0
    x(p+ps+q+1) = phi(psp+q+1);
end
if qS > 0
    x(p+ps+q+qs+1) = phi(psp+qsq+1);
end
x2 = x;

chkma = chkroots(x(p+ps+1:end), 0, 0, q, qs, qS); %check invertibility
if chkma == 1
    %   xx=invroots(x(p+ps+1:end),0,0,q,qs,qS);     %invert roots
    %   x=[x(1:p+ps) xx];
    return
end
%end of lines added 13-06-2008


% second regression (to correct the biases)
if qt > 0
    %form new residuals and auxiliary variables
    eta = zeros(n, 1);
    xi = zeros(n, 1);
    a(1) = yc(1);
    eta(1) = a(1);
    xi(1) = a(1);
    for i = 2:n
        sum = yc(i);
        sume = 0;
        sumx = 0;
        if p > 0
            pi = max(1, i-p);
            t = i - 1:-1:pi;
            nt = length(t);
            tp = 1:nt;
            sum = sum + yc(t)' * phi(tp);
            sume = sume - eta(t)' * phi(tp);
        end
        if ps > 0 && i - s > 0
            pi = max(1, i-s-p);
            t = i - s:-1:pi;
            nt = length(t);
            tp = p + 1:p + nt;
            sum = sum + yc(t)' * phi(tp);
            sume = sume - eta(t)' * phi(tp);
        end
        if q > 0
            pi = max(1, i-q);
            t = i - 1:-1:pi;
            nt = length(t);
            tp = psp + 1:psp + nt;
            sum = sum - a(t)' * phi(tp);
            sumx = sumx - xi(t)' * phi(tp);
        end
        if qs > 0 && i - s > 0
            pi = max(1, i-s-q);
            t = i - s:-1:pi;
            nt = length(t);
            tp = psp + q + 1:psp + q + nt;
            sum = sum - a(t)' * phi(tp);
            sumx = sumx - xi(t)' * phi(tp);
        end
        if qS > 0 && i - S > 0
            pi = max(1, i-S-q);
            t = i - S:-1:pi;
            nt = length(t);
            tp = psp + qsq + 1:psp + qsq + nt;
            sum = sum - a(t)' * phi(tp);
            sumx = sumx - xi(t)' * phi(tp);
        end
        if qs > 0 && qS > 0 && i - s - S > 0
            pi = max(1, i-s-S-q);
            t = i - s - S:-1:pi;
            nt = length(t);
            tp = psp + qsq + qSq + 1:psp + qsq + qSq + nt;
            sum = sum - a(t)' * phi(tp);
            sumx = sumx - xi(t)' * phi(tp);
        end
        a(i) = sum;
        eta(i) = sume + sum;
        xi(i) = sumx + sum;
    end
    Y = [];
    for i = 1:p
        Y = [Y, [zeros(i, 1); -eta(1:n-i)]];
    end
    if ps > 0
        for i = s:s + p
            Y = [Y, [zeros(i, 1); -eta(1:n-i)]];
        end
    end
    if qt > 0
        for i = 1:q
            Y = [Y, [zeros(i, 1); xi(1:n-i)]];
        end
        if qs > 0
            for i = s:s + q
                Y = [Y, [zeros(i, 1); xi(1:n-i)]];
            end
        end
        if qS > 0
            for i = S:S + q
                Y = [Y, [zeros(i, 1); xi(1:n-i)]];
            end
        end
        if qs > 0 && qS > 0
            for i = s + S:s + S + q
                Y = [Y, [zeros(i, 1); xi(1:n-i)]];
            end
        end
    end
    phis = bols(a, Y);
    phi = phi + phis;
end

x = zeros(1, p+ps+q+qs+qS);
if p > 0
    x(1:p) = phi(1:p);
end
if ps > 0
    x(p+1:p+ps) = phi(p+1:p+ps);
end
if q > 0
    x(p+ps+1:p+ps+q) = phi(psp+1:psp+q);
end
if qs > 0
    x(p+ps+q+1) = phi(psp+q+1);
end
if qS > 0
    x(p+ps+q+qs+1) = phi(psp+qsq+1);
end
