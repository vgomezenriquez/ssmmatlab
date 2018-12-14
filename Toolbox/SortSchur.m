function [Q, R, ap] = SortSchur(Q, R, z);

ap = [];
r = find(abs(diag(R, -1)) > 100*eps);
n = size(R, 1);
s = 1:n + 1;
s(r+1) = [];
N = length(s) - 1;
for k = 1:N;
    w = s(k):s(k+1) - 1;
    if length(w) == 2
        co = 1;
        si = 0;
        X = R(w, w);
        if X(1, 1) ~= X(2, 2);
            tau = (X(1, 2) + X(2, 1)) / (X(1, 1) - X(2, 2));
            off = sqrt(tau^2+1);
            ve = [tau - off, tau + off];
            [d, we] = min(abs(ve));
            co = 1 / sqrt(1+ve(we)^2);
            si = ve(we) * co;
        end
        Qs = [co, -si; si, co];
        R(:, w) = R(:, w) * Qs;
        R(w, :) = Qs' * R(w, :);
        Q(:, w) = Q(:, w) * Qs;
    end
end
for k = 1:N;
    p(k) = min(abs(z-eig(R(s(k):s(k+1)-1, s(k):s(k+1)-1))));
end
swl = [];
for k = 1:N - 1;
    [dummy, j] = max(p(k:N));
    for l = j + k - 2:-1:k;
        p(l+[0, 1]) = p(l+[1, 0]);
        swl = [swl, l];
    end
end
for k = swl;
    v = s(k):s(k+1) - 1;
    w = s(k+1):s(k+2) - 1;
    nrA = norm(R([v, w], [v, w]), inf);
    [p, q] = size(R(v, w));
    Ip = eye(p);
    Iq = eye(q);
    r = [];
    for j = 1:q
        r = [r; R(v, w(j))];
    end
    K = kron(Iq, R(v, v)) - kron(R(w, w)', Ip);
    Pe = [];
    Qe = [];
    for ka = 1:p * q - 1;
        [a, er] = max(abs(K(ka:p*q, ka:p*q)));
        [dummy, ce] = max(abs(a));
        cl = ce + ka - 1;
        rw = er(ce) + ka - 1;
        K([ka, rw], :) = K([rw, ka], :);
        K(:, [ka, cl]) = K(:, [cl, ka]);
        Pe(ka) = rw;
        Qe(ka) = cl;
        if K(ka, ka) ~= 0;
            rs = ka + 1:p * q;
            K(rs, ka) = K(rs, ka) / K(ka, ka);
            K(rs, rs) = K(rs, rs) - K(rs, ka) * K(ka, rs);
        end
    end
    H = tril(K')';
    L = tril(K, -1) + eye(p*q);
    e = min(abs(diag(H)));
    sigp = 1:p * q;
    for ka = 1:p * q - 1;
        sigp([ka, Pe(ka)]) = sigp([Pe(ka), ka]);
    end
    r = e * r(sigp);
    x = (H \ (L \ r));
    sigq = 1:p * q;
    for ka = 1:p * q - 1;
        sigq([ka, Qe(ka)]) = sigq([Qe(ka), ka]);
    end
    x(sigq) = x;
    X = [];
    for ka = 1:q
        X = [X, x((ka - 1)*p+1:ka*p)];
    end
    [Qr, Rr] = qr([-X; e * Iq]);
    R(:, [v, w]) = R(:, [v, w]) * Qr;
    R([v, w], :) = Qr' * R([v, w], :);
    Q(:, [v, w]) = Q(:, [v, w]) * Qr;
    s(k+1) = s(k) + s(k+2) - s(k+1);
    v = s(k):s(k+1) - 1;
    w = s(k+1):s(k+2) - 1;
    if length(v) == 2
        co = 1;
        si = 0;
        X = R(v, v);
        if X(1, 1) ~= X(2, 2);
            tau = (X(1, 2) + X(2, 1)) / (X(1, 1) - X(2, 2));
            off = sqrt(tau^2+1);
            ve = [tau - off, tau + off];
            [d, we] = min(abs(ve));
            co = 1 / sqrt(1+ve(we)^2);
            si = ve(we) * co;
        end
        Qs = [co, -si; si, co];
        R(:, v) = R(:, v) * Qs;
        R(v, :) = Qs' * R(v, :);
        Q(:, v) = Q(:, v) * Qs;
    end
    if length(w) == 2
        co = 1;
        si = 0;
        X = R(w, w);
        if X(1, 1) ~= X(2, 2);
            tau = (X(1, 2) + X(2, 1)) / (X(1, 1) - X(2, 2));
            off = sqrt(tau^2+1);
            ve = [tau - off, tau + off];
            [d, we] = min(abs(ve));
            co = 1 / sqrt(1+ve(we)^2);
            si = ve(we) * co;
        end
        Qs = [co, -si; si, co];
        R(:, w) = R(:, w) * Qs;
        R(w, :) = Qs' * R(w, :);
        Q(:, w) = Q(:, w) * Qs;
    end
    ap(k) = norm(R(w, v), inf) / (10 * eps * nrA);
end
R = R - tril(R, -2);
for j = 2:N;
    R(s(j), s(j)-1) = 0;
end
