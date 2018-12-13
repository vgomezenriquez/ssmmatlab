function [abar, bbar, cbar, t, k] = mobsvf(a, b, c, tol)
%OBSVF  Observability staircase form.
%
%   [ABAR,BBAR,CBAR,T,K] = MOBSVF(A,B,C) returns a decomposition
%   into the observable/unobservable subspaces.
%
%   If the observability matrix, Ob=OBSV(A,C), has rank r <= n =
%   SIZE(A,1), then there is a similarity transformation T such that
%
%      Abar = T * A * T' ,  Bbar = T * B  ,  Cbar = C * T' .
%
%   and the transformed system has the form
%
%          | Ano   A12|           |Bno|
%   Abar =  ----------  ,  Bbar =  ---  ,  Cbar = [ 0 | Co].
%          |  0    Ao |           |Bo |
%
%                                              -1           -1
%   where (Ao,Bo) is controllable, and Co(sI-Ao) Bo = C(sI-A) B.
%
%   The last output K is a vector of length n containing the
%   number of observable states identified at each iteration
%   of the algorithm.  The number of observable states is SUM(K).
%

if nargin < 4
    tol = [];
end

[aa, bb, cc, t, k] = mctrbf(a', c', b', tol);
abar = aa';
bbar = cc';
cbar = bb';
