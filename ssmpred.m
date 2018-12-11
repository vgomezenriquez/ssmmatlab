function [ypr, mypr, alpr, malpr] = ssmpred(n, p, A, P, X, Z, G, W, T, H, g, M)
%
%
%        This function computes n forecasts of the state vector and the
%        observations corresponding to the model
%
%        y_t = X_t*beta + Z_t*alpha_t + G_t*epsilon_t
%        alpha_{t+1}= W_t*beta + T_t*alpha_t + H_t*epsilon_t,
%
%        where epsilon_t is (0,sigma^2I),
%
%        with initial state
%
%        alpha_1= c + W_0*beta + a_1 + A_1*delta
%
%        where c is (0,Omega) and delta is (0,kI) (diffuse).
%
%        Input parameters:
%        n  : number of forecasts
%        p  : number of variables
%        A  : the estimated augmented state vector at the end of filtering
%        P  : the Mse of A at the end of filtering
%        X  : an (n*p x nbeta) matrix containing the X_t matrices if it
%             is time-varying. A (p x nbeta) if it is time invariant.
%             It can be []
%        Z  : an (n*p x nalpha) matrix containing the Z_t matrices if it
%             is time-varying. A (p x nalpha) matrix if it is time
%             invariant
%        G  : an (n*p x nepsilon) matrix containing the G_t matrices if
%             it is time-varying. A (p x nepsilon) matrix if it is time
%             invariant
%        W  : an (n*nalpha x nbeta) matrix containing the W_t matrices if
%             is time-varying. An (nalpha x nbeta) matrix if it is time
%             invariant. It can be []
%        T  : an (n*nalpha x nalpha) matrix containing the T_t matrices
%             if it is time-varying. An (nalpha x nalpha) matrix if it time
%             invariant.
%        H  : an (n*nalpha x nepsilon) matrix containing the H_t matrices
%             if it is time-varying. An (nalpha x nepsilon) if it is time
%             invariant.
%        g   : the beta estimator
%        M   : the Mse of the beta estimator
%
%        Output parameters:
%        ypr : a (p x n) matrix containing the forecasts of the
%              observations
%        mypr: a (p x p x n) array containing the Mse of forecasts
%        alpr: an (nalpha x n) matrix containing the forecasts of the state
%              vector
%       malpr: an (nalpha x nalpha x n) array containing the Mse of alpr
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
[mx, nx] = size(X);
[mz, nalpha] = size(Z);
[mg, neps] = size(G);
[mw, nw] = size(W);
[mt, nt] = size(T);
[mh, nh] = size(H);
nbeta = max(nx, nw);
alpr = zeros(nalpha, n);
malpr = zeros(nalpha, nalpha, n);
ypr = zeros(p, n);
mypr = zeros(p, p, n);
%
% check for inconsistencies
%
%this part added on 24-3-2017
np = n * p;
nnalpha = n * nalpha;
%end of addition

if mx > p | mz > p | mg > p | mw > nalpha ...
        | mt > nalpha | mh > nalpha
    %
    % system matrices are time varying
    %
    if mx > p & mx ~= np
        disp('the number of rows in X is incorrect')
        return
    end
    if mz > p & mz ~= np
        disp('the number of rows in Z is incorrect')
        return
    end
    if mg > p & mg ~= np
        disp('the number of rows in G is incorrect')
        return
    end
    if mw > nalpha & mw ~= nnalpha
        disp('the number of rows in W is incorrect')
        return
    end
    if mt > nalpha & mt ~= nnalpha
        disp('the number of rows in T is incorrect')
        return
    end
    if mh > nalpha & mh ~= nnalpha
        disp('the number of rows in H is incorrect')
        return
    end
end
if nw ~= nx & nw ~= 0 & nx ~= 0
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

%
% augmented Kalman filter
%
%A
%P
[junk, mA] = size(A);
ndelta = mA - nbeta - 1;
%
for i = 1:n
    %   i
    ip = (i - 1) * p + 1:i * p;
    ia = (i - 1) * nalpha + 1:i * nalpha;
    if mx > p
        XX = X(ip, :);
    else
        XX = X;
        if nbeta > 0 & isempty(X) == 1
            XX = zeros(p, nbeta);
        end
    end
    if mz > p
        ZZ = Z(ip, :);
    else
        ZZ = Z;
    end
    if mg > p
        GG = G(ip, :);
    else
        GG = G;
    end
    if mw > nalpha
        WW = W(ia, :);
    else
        WW = W;
        if nbeta > 0 & isempty(W) == 1
            WW = zeros(nalpha, nbeta);
        end
    end
    if mt > nalpha
        TT = T(ia, :);
    else
        TT = T;
    end
    if mh > nalpha
        HH = H(ia, :);
    else
        HH = H;
    end
    V = [zeros(p, ndelta), full(XX)] - ZZ * A(:, 1:end-1);
    F = ZZ * P * ZZ' + GG * GG';
    if ~isempty(g)
        alpr(:, i) = A(:, end) - A(:, 1:end-1) * g;
        malpr(:, :, i) = P + A(:, 1:end-1) * M * A(:, 1:end-1)';
        ypr(:, i) = V * g + ZZ * A(:, end);
        mypr(:, :, i) = F + V * M * V';
    else
        alpr(:, i) = A;
        malpr(:, :, i) = P;
        ypr(:, i) = ZZ * A(:, end);
        mypr(:, :, i) = F;
    end
    
    %next forecast
    if i < n
        A = [zeros(nalpha, ndelta), -full(WW), zeros(nalpha, 1)] + TT * A;
        P = TT * P * TT' + HH * HH';
    end
end
