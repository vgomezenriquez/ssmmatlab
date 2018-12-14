function [pj, qj] = distnj(n, j, p, q)
%      given p and q of a signature matrix J=diag(I_p,-I_q) such that n=p+q and j<=n,
%      this function returns integers pj and qj such that diag(I_pj,-I_qj)
%      is the submatrix of J formed with the last n-j+1 rows and columns. Any
%      of p, q, pj or qj can be zero.
%
%---------------------------------------------------
% USAGE: [pj,qj]=distnj(n,j,p,q)
% where:    n = integer
%           j = integer <= n
%           p,q= integers such that J = diag(I_p,-I-q) is a signature
%           matrix with n=p+q
%
%---------------------------------------------------
% RETURNS: pj and qj,
%---------------------------------------------------
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
if n ~= p + q
    error('n nonequal to p + q in distnj')
end
if j - 1 <= p
    pj = p - j + 1;
    qj = q;
else
    pj = 0;
    qj = n - j + 1;
end