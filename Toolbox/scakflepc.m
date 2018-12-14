function [e, f, hb, Mb, A, P, qyy, R] = scakflepc(y, X, Z, G, W, T, H, ins, i, chb)
%
%
%       This function applies the two stage Kalman filter to the series y
%       for prediction and profile likelihood evaluation corresponding to
%       the model
%
%       y_t = X_t*beta + Z*alpha_t + G*epsilon_t
%       alpha_{t+1}= W_t*beta + T*alpha_t + H*epsilon_t,
%
%       where epsilon_t is (0,I),
%
%       with initial state
%
%       alpha_1 = delta
%
%       that is considered fixed and unknown. It is assumed that the state
%       space model is in innovations form and P_t = 0 for all t. That is,
%       G=Sigma^{1/2} and H=K*Sigma^{1/2}, where Sigma =
%       Sigma^{1/2}*Sigma^{1/2 '} is the Cholesky decomposition of the
%       covariance matrix of the innovations. Therefore, the recursions are
%       simplified and only the innovations and the state estimators need
%       to be updated. No collapsing takes place. No missing values are
%       allowed.
%
%       Input parameters:
%       y:     an (n x p) matrix of observations;
%       X    : a  (p x nbeta) if it is time invariant;
%              it can be []
%       Z    : a  (p x nalpha) matrix if it is time invariant
%       G    : a  (p x nepsilon) matrix if it is time invariant
%       W    : an (n*nalpha x nbeta) matrix containing the W_t matrices;
%              an (nalpha x nbeta) matrix if it is time invariant;
%              it can be []
%       T    : an (nalpha x nalpha) matrix if it time invariant
%       H    : an (nalpha x nepsilon) if it is time invariant
%       ins: an nalpha x (cc+cw0+ca1+cca1) matrix containing the initial
%            state information, according to array i below
%       i    : a  1 x 4 array containing 4 integers, i=[cc cw0 ca1 cca1],
%            where
%            cc   = 0
%            cw0  = 0
%            ca1  = 0
%            cca1 = nalpha
%       chb= 1 compute hb and Mb
%            0 do not compute hb and Mb
%
%       Output parameters:
%       e   : residual vector (Q'_2*y)
%       f   : factor by which the residuals are to be multiplied
%             for minimization of the nonlinear sum of squares
%       hb  : the beta estimator
%       Mb  : the Mse of the beta estimator
%        A  : the estimated augmented state vector at the end of filtering
%        P  : the Mse of A at the end of filtering
%      qyy  : Q'_1*y in the QR decomposition to estimate beta
%        R  : R in the QR decomposition to estimate beta
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
% system dimensions
%
[n, p] = size(y);
[mx, nx] = size(X);
[mz, nalpha] = size(Z);
[mg, neps] = size(G);
[mw, nw] = size(W);
[mt, nt] = size(T);
[mh, nh] = size(H);
[mi, ni] = size(ins);
[mc, nc] = size(i);
nbeta = max(nx, nw);
cc = i(1);
cw0 = i(2);
ca1 = i(3);
cca1 = i(4);
%
% check for inconsistencies
%
np = n * p;
nnalpha = n * nalpha;
if nbeta > np
    disp('too many regressors in scakfle2')
    return
end
if mi ~= nalpha
    disp('the number of rows in ins is incorrect')
    return
end
if mx > p || mz > p || mg > p || mw > nalpha ...
        || mt > nalpha || mh > nalpha
    %
    % system matrices are time varying
    %
    if mx > p && mx ~= np
        disp('the number of rows in X is incorrect')
        return
    end
    if mz > p && mz ~= np
        disp('the number of rows in Z is incorrect')
        return
    end
    if mg > p && mg ~= np
        disp('the number of rows in G is incorrect')
        return
    end
    if mw > nalpha && mw ~= nnalpha
        disp('the number of rows in W is incorrect')
        return
    end
    if mt > nalpha && mt ~= nnalpha
        disp('the number of rows in T is incorrect')
        return
    end
    if mh > nalpha && mh ~= nnalpha
        disp('the number of rows in H is incorrect')
        return
    end
end
if nw ~= nx && nw ~= 0 && nx ~= 0
    disp('the number of columns in W and X should agree')
    return
end
if nt ~= nalpha
    disp('the number of columns in T and Z should agree')
    return
end
if nh ~= neps
    disp('the number of columns in G and H should agree')
    return
end
if mi ~= nalpha
    disp('the number of rows in ins and the number of columns in Z should agree')
    return
end
if nc ~= 4 || mc ~= 1
    disp('c should be a 1 x 4 matrix')
    return
end
if nbeta == 0 && cw0 ~= 0
    disp('the number of rows in W_0 should be zero')
    return
end
if ni ~= cc + cw0 + ca1 + cca1
    disp('the number of columns in ins should be ')
    disp(cc+cw0+ca1+cca1)
    return
end


if (mg ~= p) || (mg ~= neps)
    disp('matrix G should be lower triangular in scakflepc')
    return
else
    if ((G(1, 1) - 1.) ~= 0)
        disp('matrix G should have a one in the (1,1) position in scakflepc')
        return
    end
end


%
% initial covariance matrix
%
P = zeros(nalpha, nalpha);
%
% initial state
%
if cw0 == 0
    if nbeta == 0
        W0 = [];
    else
        W0 = zeros(nalpha, nbeta);
    end
end
if ca1 == 0
    a1 = zeros(nalpha, 1);
end
if cca1 == 0
    disp('cca1 should not be zero in scakflepc')
    return
else
    ndelta = cca1;
    aa1 = ins(:, cc+cw0+ca1+1:cc+cw0+ca1+cca1);
end
A = [-aa1, -full(W0), a1];

%
% augmented Kalman filter
%
SQT = zeros(np, ndelta+nbeta+1);
f = 1;
if p == 1
    %series is univariate
    for i = 1:n
        %   i
        ip = (i - 1) * p + 1:i * p;
        if mx > p
            XX = X(ip, :);
        else
            XX = X;
            if nbeta > 0 && isempty(X) == 1
                XX = zeros(p, nbeta);
            end
        end
        ZZ = Z;
        if mw > nalpha
            ia = (i - 1) * nalpha + 1:i * nalpha;
            WW = W(ia, :);
        else
            WW = W;
            if nbeta > 0 && isempty(W) == 1
                WW = zeros(nalpha, nbeta);
            end
        end
        TT = T;
        K = H;
        YY = y(i, :)';
        V = [zeros(p, ndelta), full(XX), YY] - ZZ * A;
        DD = V;
        %
        % SQT updating
        %
        SQT(ip, :) = DD;
        A = [zeros(nalpha, ndelta), -full(WW), zeros(nalpha, 1)] + TT * A + K * V;
    end
else
    %series is multivariate
    for i = 1:n
        %   i
        ip = (i - 1) * p + 1:i * p;
        if mx > p
            XX = X(ip, :);
        else
            XX = X;
            if nbeta > 0 && isempty(X) == 1
                XX = zeros(p, nbeta);
            end
        end
        ZZ = Z;
        if mw > nalpha
            ia = (i - 1) * nalpha + 1:i * nalpha;
            WW = W(ia, :);
        else
            WW = W;
            if nbeta > 0 && isempty(W) == 1
                WW = zeros(nalpha, nbeta);
            end
        end
        TT = T;
        K = H;
        YY = y(i, :)';
        V = [zeros(p, ndelta), full(XX), YY] - ZZ * A;
        DD = G \ V;
        %
        % SQT updating
        %
        SQT(ip, :) = DD;
        A = [zeros(nalpha, ndelta), -full(WW), zeros(nalpha, 1)] + TT * A + K * DD;
    end
end
%
% computation of residuals, hat beta and its Mse
%
nbeta = nbeta + ndelta;
[ns, ms] = size(SQT);
if ms < nbeta
    error('Singular model in SCAKFLEPC')
end
[Q, R] = qr(SQT(:, 1:nbeta));
qy = Q' * SQT(:, nbeta+1);
qyy = qy(1:nbeta);
e = qy(nbeta+1:ns);
if chb == 1
    %design matrix may be singular sometimes
    M = R(1:nbeta, :);
    Mb = pinv(M);
    hb = Mb * qy(1:nbeta);
    Mb = Mb * Mb';
    if p == 1
        f = 1.;
    else
        f = (prod(diag(G)))^(1 / p);
    end
end
