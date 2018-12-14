function [abar, bbar, cbar, t, k] = mctrbf(a, b, c, tol)
%CTRBF  Controllability staircase form.
%
%   [ABAR,BBAR,CBAR,T,K] = MCTRBF(A,B,C) returns a decomposition
%   into the controllable/uncontrollable subspaces.
%
%   If the controllability matrix, Co=CTRB(A,B), has rank r <= n
%   = SIZE(A,1), then there is a similarity transformation T such
%   that
%
%      Abar = T * A * T' ,  Bbar = T * B ,  Cbar = C * T'
%
%   and the transformed system has the form
%
%             | Anc    0 |           | 0 |
%      Abar =  ----------  ,  Bbar =  ---  ,  Cbar = [Cnc| Cc].
%             | A21   Ac |           |Bc |
%                                              -1          -1
%   where (Ac,Bc) is controllable, and Cc(sI-Ac)Bc = C(sI-A)B.
%
%   The last output K is a vector of length n containing the
%   number of controllable states identified at each iteration
%   of the algorithm.  The number of controllable states is SUM(K).
%

if nargin < 4
    tol = [];
end

[a1, b1, k, u] = stair(a, b, tol);
t = flipdim(u', 1);
k = flipdim(k, 1);
abar = t * a * t';
bbar = t * b;
cbar = c * t';
