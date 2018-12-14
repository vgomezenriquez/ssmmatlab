function [P, K, Sigma, U, iU] = ss2if(F, G, H, J, Q, S, R)
%
% This function transforms a general SS model into an innovations form
% solving the DARE. The state space form is
%
%     x_{t+1}  =  Fx_{t} + Gu_{t}
%     Y_{t}    =  Hx_{t} + Jv_t,
%
%   E[u_t]
%    [v_t][u'_s, v'_s]
%    =[Q  S]
%     [S' R]\delta_{ts}
%
%   P  = Solution of the DARE (discrete algebraic Riccati equation) =
%   MSE(\hat{x}_{T}).
%   K  = Kalman gain in the innovations model
%     \hat{x}_{t+1}  =  F\hat{x}_{t} + Ka_{t}
%     Y_{t}          =  H\hat{x}_{t} + a_{t},
%  Sigma = covariance matrix of the innovations = Var(a_{t})
%   U = Sigma^{1/2}
%   iU = Sigma^{-1/2}
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


l = size(F, 1);
m = size(H, 1);
N = J * R * J';
iN = pinv(J*R*J');
GSJt = G * S * J';
M = G * Q * G' - GSJt * iN * GSJt';
tol = max(size(M)) * norm(M) * eps(double(M));
if abs(trace(M)) < tol
    P = zeros(l);
else
    A = F - GSJt * iN * H;
    
    [U, S, V] = svd([N; -H']);
    T = [U(m+1:l+m, m+1:l+m)' * A', zeros(l); -M, eye(l)];
    Y = [U(m+1:l+m, m+1:l+m)', U(1:m, m+1:l+m)' * H; zeros(l), A];
    
    
    [W, d] = eig(T, Y);
    d = diag(d);
    [e, index] = sort(abs(d));
    WW = W(:, index(1:l));
    P = real(WW(l+1:2*l, :)*pinv(WW(1:l, :)));
end

if nargout > 1
    Sigma = H * P * H' + N;
    [Uf, SVf, Vf] = svd(Sigma);
    SSVf = diag(sqrt(diag(SVf)));
    U = Uf * SSVf * Vf'; % \Sigma^{1/2}
    iU = Vf * pinv(SSVf) * Uf'; % \Sigma^{-1/2}
    K = (F * P * H' + GSJt) * (iU' * iU);
end
