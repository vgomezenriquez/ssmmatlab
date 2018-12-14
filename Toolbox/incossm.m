function [ins, i, ferror, P, V, XX] = incossm(T, H, ndelta)
%
%
%        This function obtains the initial conditions corresponding to the
%        model
%
%        y_t = X*beta + Z*alpha_t + G*epsilon_t
%        alpha_{t+1}= W*beta + T*alpha_t + H*epsilon_t
%
%        where epsilon_t is (0,sigma^2I),
%
%        with initial state
%
%        alpha_1= c + W_0*beta + a_1 + A_1*delta
%
%        where c is (0,Omega) and delta is (0,kI) (diffuse). It is assumed
%        that W_0=[] and, therefore, cw0 below is zero.
%
%        Input parameters:
%        T    : an (nalpha x nalpha) matrix
%        H    : an (nalpha x nepsilon) matrix
%      ndelta : a positive integer, the number of diffuse components in
%               alpha_t (= number of eigenvalues of unit module in T)
%
%        Output parameters:
%        ins: an nalpha x (cc+cw0+ca1+cca1) matrix containing the initial
%             state information, according to array i below
%        i    : a  1 x 4 array containing 4 integers, i=[cc cw0 ca1 cca1],
%             where
%             cc   = nalpha if c is not missing (0 if c missing)
%             cw0  = number of columns in W_0 (0 if W_0 missing)
%             ca1  = 1 if a_1 is not missing (0 if a_1 missing)
%             cca1 = number of columns in A_1 (0 if A_1 missing)
%      ferror : flag for errors
%        P    : first output matrix of SortSchur function
%        V    : solution of the Lyapunov equation corresponding to the
%               stationary part, V = Fs*V*Fs' + Hs*Hs'
%       XX    : solution of the continuous time Lyapunov equation
%               Fn*XX - XX*Fs = A
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

i = [];
ins = [];
V = [];
ferror = 0;
nalpha = size(T, 1);
if (nalpha == ndelta)
    i = [0, 0, 0, ndelta];
    ins = eye(ndelta);
else
    [P, S] = schur(T);
    if nalpha > 1
        z = 0.;
        [P, S, ~] = SortSchur(P, S, z);
    end
    % E=ordeig(S);  [P,S]=ordschur(P,S,abs(E)>.99);
    Fn = S(1:ndelta, 1:ndelta);
    Fs = S(ndelta+1:end, ndelta+1:end);
    A = S(1:ndelta, ndelta+1:end);
    [XX, ferror] = mclyapunov(Fn, Fs, A);
    if ferror > 0
        return
    end
    PpH = P' * H;
    Hs = PpH(ndelta+1:end, :);
    Pn = P(:, 1:ndelta);
    Ps = P(:, ndelta+1:end);
    Rn = Pn;
    Rs = Ps - Pn * XX;
    A = Rn;
    V = dlyapsq(Fs, Hs);
    V = V' * V;
    Vs = Rs * V * Rs';
    i = [nalpha, 0, 0, ndelta];
    ins = [Vs, A];
end
