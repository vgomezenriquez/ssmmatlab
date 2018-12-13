function [e, E, rSigmat] = sqrt_ckms(y, Y, str, maxupdt, tol)
%
%
%        This function applies to the series (y,Y) the square root CKMS
%        recursions corresponding to the model
%
%        y_t = Y_t*beta + H*x_t + K*a_t
%        x_{t+1}= F*x_t + a_t,
%
%        where a_t is (0,sigmar2),
%
%        with initial state
%
%        alpha_1= (0,Omega).
%
%        Input parameters:
%        y  : an n x p matrix of observations
%        Y  : an n x (p x nbeta) matrix containing the Y_t matrices; it can be []
%        str: a structure containing the model parameters
%        maxupdt : an integer parameter equal to the maximum number of updates
%                  before passing to the steady state recursions. If
%                  maxupdt is empty, it is made equal to the number of
%                  observations.
%        tol: a real parameter that controls when to pass to the steady
%             state recursions. This happens when the difference between the
%             (1,1) element of the covariance square root of the filter and the
%             (1,1) element of the covariance square root of the model
%             innovations is less than or equal to tol.
%
%        Output parameters:
%      (e,E)   : augmented residual vector
%    rSigmat : the square roots of the residual covariance matrices
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

e = [];
E = [];
rSigmat = [];
%
% system dimensions
%
[n, p] = size(y);
[my, ny] = size(Y);
nbeta = ny;
if isempty(maxupdt)
    maxupdt = n;
end

kro = str.kro;
so = sum(kro); %system order (McMillan degree)
mkro = max(kro);
F = str.Fs;
H = str.Hs;

if isfield(str, 'cov')
    c = str.cov;
    %compare with this for steady state
    target = 1.e10; %we would need a spectral factorization to obtain target
else
    Sigma = str.sigmar2;
    phi = str.phis;
    theta = str.thetas;
    %compute covariances of y
    [c, ierror] = macgf(phi, theta, Sigma, mkro+1);
    %compare with this for steady state
    target = abs(sqrt(Sigma(1, 1)));
end

stc = []; %stack covariances
for i = 1:mkro
    stc = [stc; c(:, :, i+1)];
end
%set up matrix N
N = [];
for i = 1:p
    for k = 1:kro(i)
        N = [N; stc((k - 1)*p+i, :)];
    end
end
%initial Sigmat and barK can be obtained using covariance data only
Sigmat = c(:, :, 1); %initial innovations cov. matrix
barK = N; %initial barK matrix
urSigmat = chol(Sigmat); %initial square root of innovations cov. matrix
barKp = barK / urSigmat; %initial barKp matrix
barL = barKp; %initial barL matrix
%initial state vector
A = zeros(so, nbeta+1);
% F=sparse(F); H=sparse(H);


steady = 0;
dif = 1.e10;
for i = 1:n
    %  i
    ip = (i - 1) * p + 1:i * p;
    YY = y(i, :)';
    if my > p
        XX = Y(ip, :);
    else
        XX = Y;
    end
    %update innovation
    V = [full(XX), YY] - mulHA(H, A, kro);
    e = [e; V(:, end)'];
    E = [E; V(:, 1:end-1)'];
    %update state prediction
    A = mulFA(F, A, kro) + barKp * (urSigmat' \ V);
    rSigmat = [rSigmat; urSigmat'];
    if steady == 0
        M = [urSigmat; (mulHA(H, barL, kro))'];
        [Q, R] = jqrt(M, p, p);
        % % the expressions M'*J*M and R'*J*R should coincide in the following
        %   J=eye(2*p); J(p+1:end,p+1:end)=-eye(p);
        %   M'*J*M, R'*J*R
        %update square root of innovations covariance matrix
        nurSigmat = R(1:p, :);
        if (i <= maxupdt) && (dif > tol)
            %     format long g
            %     dif= abs(abs(nurSigmat(1,1))-target)
            %     i,target,nurSigmat(1,1)
            %     format short
            urSigmat = nurSigmat;
            Mx = [barKp'; (mulFA(F, barL, kro))'];
            QtB = qtb(Q, Mx, p, p);
            %update barKp and barL matrices
            barKp = QtB(1:p, :)';
            barL = QtB(p+1:end, :)';
        else
            steady = 1; %no more updates, steady state instead
        end
        dif = abs(abs(nurSigmat(1, 1))-target);
    end
end
