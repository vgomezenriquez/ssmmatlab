function [phi, theta, ierror] = ss2lcvarmaxf(H, F, K, G)
%
% This function computes a left coprime VARMAX model corresponding to
% a system (H,F,K,G),
%   x(t+1) = F*x(t) + K*u(t)
%    y(t)  = H*x(t) + G*u(t),
% where H is a k x n vector, F is a square n x n matrix, K is an n x m matrix
% and G is a k x m. The system (H,F,K,G) is not necessarily minimal.
%
% Input:  H      = a k x n matrix
%         F      = an n x n matrix
%         K      = an n x m matrix
%         G      = a k x m matrix
% If G is not given, it is assumed that m=k and G = I_k.
% Output: phi    = a k x k x l matrix containing the AR matrix polynomial
%         theta  = a k x m x l matrix containing the MA matrix polynomial,
%                  l = max(n_i:i=1,...,k), where n_i are the Kronecker
%                  indices.
%      ierror    = 0,1  a flag for errors in dimensions
% The method uses the algorithm of Chen to pass from right to left MFD. It
% is based on the decomposition
%
%   Y(t) = [zH(I-zF)^{-1}K + G]u(t) = [S(z)R^{-1}(z)K + G]u(t)
%        = [Rb^{-1}(z)Sb(z)K + G]u(t),
%
% that implies
%
%  Rb(z)Y(t) = [Sb(z)K + Rb(z)G]u(t)
%
% or
%
%  phi(z)Y(t) = theta(z)u(t),
%
% where phi(z)=Rb(z) and theta(z)=Sb(z)K + Rb(z)G. This method is used in
% Kucera (1991), p. 226.
%
% Copyright (c) 21 July 2003 by Victor Gomez
% Ministerio de Hacienda, Direccion Gral. de Presupuestos,
% Subdireccion Gral. de Analisis y P.E.,
% Alberto Alcocer 2, 1-P, D-34, 28046, Madrid, SPAIN.
% Phone : +34-915835439
% E-mail: vgomez@sepg.minhac.es
%
% The author assumes no responsibility for errors or damage resulting from the use
% of this code. Usage of this code in applications and/or alterations of it should
% be referenced. This code may be redistributed if nothing has been added or
% removed and no money is charged. Positive or negative feedback would be appreciated.
%
%
[nh, mh] = size(H);
[nf, mf] = size(F);
[nk, mk] = size(K);
phi = [];
theta = [];
ierror = 0;
iG = 0;
if nargin < 4
    if (mh ~= nf) | (nf ~= mf) | (nf ~= nk) | (nh ~= mk)
        disp('wrong dimensions of H, K or F in ss2lcvarmaxf');
        ierror = 1;
        return
    end
    m = nh;
    G = eye(m);
else
    iG = 1;
    [ng, mg] = size(G);
    if (mk ~= mg) | (nh ~= ng) | (mh ~= nf) | (nf ~= mf) | (nf ~= nk)
        disp('wrong dimensions of H, K, F or G in ss2lcvarmaxf');
        ierror = 1;
        return
    end
end

n = nf;

S(:, :, 1) = zeros(size(H));
S(:, :, 2) = H;

R(:, :, 1) = eye(n);
R(:, :, 2) = -F;


%pass from right to left coprime MFD
[phie, thetae, kro, ierror] = pright2leftcmfd(R, S, n+1);
if ierror > 0
    disp('error in function pright2leftcmfd')
    ierror
    return
end

%obtain phi and theta polynomial matrices
%phi
phi = phie;
%theta
Ks(:, :, 1) = K;
Gs(:, :, 1) = G;
[theta, ierror] = pmatmul(thetae, Ks);
theta = theta + pmatmul(phie, Gs);