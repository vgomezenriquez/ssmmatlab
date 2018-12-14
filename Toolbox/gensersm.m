function B = gensersm(phi, alpha, th, stda, ctr, Ns, N, l, seed)
%
% This function obtains Ns simulated series of length N that follow the
% model defined by phi, alpha, th, stda and ctr.
%
% input arguments:
% phi: AR polynomial
% alpha: nonstationary polynomial (differencing operator)
% th: MA polynomial
% stda: standard deviation of the innovations in percentage of the levels.
% ctr: = 0 include neither a constant nor a time trend in the model
%      = 1 include a constant in the model
%      = 2 include a constant plus a time trend in the model
% Ns: number of simulations
% N: series length
% l: number of initial observations discarded in the simulated series
%
% output arguments:
% B: an Ns_N x 1 vector containing the stacked simulated series
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
% set system matrices
%
% [Z,T,H,A,Sigma,Xi]=arimam(phi,alpha,th);
phi = fliplr(phi);
p = size(phi, 2) - 1;
alpha = fliplr(alpha);
d = size(alpha, 2) - 1;
th = fliplr(th);
q = size(th, 2) - 1;
phip(:, :, p+1) = 1;
for i = 1:p + 1
    phip(:, :, i) = phi(i);
end
alphap(:, :, d+1) = 1;
for i = 1:d + 1
    alphap(:, :, i) = alpha(i);
end
phip = pmatmul(phip, alphap);
thp(:, :, q+1) = 1;
for i = 1:q + 1
    thp(:, :, i) = th(i);
end
[H, F, K, J, ierror] = qarmax2ss1(phip, thp);
[A, Sigma, Xi] = vincovma(F, K, stda);
%
% loop for the simulations
%
% randn('state',seed);
% rng(seed)                          %new matlab randon number generator
if isOctave
    randn('state', seed);
else
    rng(seed); %use the same start for random number generator
end
B = zeros(Ns*N, 1);
for i = 1:Ns
    y = varimass(l, N, H, F, K, A, Sigma, Xi, stda);
    %
    % generate a constant if ctr=1
    % generate a constant plus a time trend if ctr=2
    %
    if ~isempty(ctr)
        y = ctr + y;
    end
    B(N*(i - 1)+1:N*i) = y(1:N);
end
%  save data\simser y -ascii
