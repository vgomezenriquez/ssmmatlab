function [ff, beta, e, f, str, stx, recrs] = exactmedfvc(beta, y, x, str, Y, chb)
% PURPOSE: given a structure, it computes the functions in ff such that the
% expression ff'*ff is minimized with the Levenberg-Marquardt method.
% The state space echelon form is:
%
%  alpha_{t+1} = F*alpha_{t} + B*x_t{t} + K*a_{t}
%      y_{t}   = Y_{t}*beta + H*alpha_{t} + D*x_{t}  + a_{t}
%
%---------------------------------------------------
% USAGE: [ff,beta,e,f,str,stx,recrs]=exactmedfvc(beta,y,x,str,Y,chb)
% where:    beta   = the parameter vector
%           y      = an (nobs x neqs) matrix of y-vectors
%           x      = matrix of input variables (nobs x nx)
%           str    = a structure containing the model information
%           Y      = an (nobs x (neqs x nbeta)) regression matrix
%                    (neqs x nbeta) if it is time invariant
%           chb    = 1 compute hb and Mb in Kalman filter
%                    0 do not compute hb and Mb in Kalman filter
%                    2 compute recursive residuals updating regression p.
%                    3 compute recursive residuals with fixed regression p.
%---------------------------------------------------
% RETURNS: ff   = a vector containing the individual functions at the
%                 solution
%          beta = the parameter vector, possibly modified
%          e    = a vector containing the standardized residuals
%          f    = a scalar containing the determinantal term
%          str  = the input structure str, possibly modified
%          stx  = a structure containing the following fields
%       .X,.Z,.G,.W,.T,.H,.ins,.i are the matrices and initial conditions
%       information corresponding to the state space model for the
%       endogenous part
%         .hb   = the beta estimator
%         .Mb   = the Mse of the beta estimator
%         .A    = the estimated augmented state vector at the end of
%                 filtering
%         .P    = the Mse of A at the end of filtering
%  .B,.D, .m1   = the initial state estimate and the matrices B and D for
%                 the state space model corresponding to the exogenous part
%        .kro   = the kronecker indices for the undecoupled model
%         recrs = standardized recursive residuals
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
ff = [];
e = [];
f = [];
A = [];
P = [];
recrs = [];
Mb = [];
m1 = [];
[nobs, neqs] = size(y);
[nY, nbeta] = size(Y);

vgams = str.vgam;
bind = str.bind;
[nbind, mbind] = size(bind); %nbind includes the mean parameters
nparma = nbind - neqs; %number of parameters in ar and ma parts
Lparm = beta(nparma+1:end); %save parameters for the Cholesky factor L
for i = 1:nparma
    vgams(bind(i)) = beta(i);
end
str.vgams = vgams;
str = param2sse(str);

%check stationarity and invertibility. If necessary, change parameters
iar = chkstainv(str.Fs); %if iar >1, the model is not stationary
if iar > 1
    %  fprintf(1,'model nonstationary, iar = %2d\n',iar);
    %convert Phi(z) into Phi(lambda*z) for an appropriate lambda
    vgam = enfstab(str, 'phi  ');
    str.vgams = vgam;
    str = param2sse(str);
    for i = 1:nparma
        beta(i) = vgam(bind(i));
    end
    beta(nparma+1:end) = Lparm; %insert parameters for the Cholesky factor
end
ima = chkstainv(str.Fs-str.Ks*str.Hs); %if ima >1, the model is not invertible
if ima > 1
    %  fprintf(1,'model noninvertible, ima = %2d\n',ima);
    %convert Th(z) into Th(lambda*z) for an appropriate lambda
    vgam = enfstab(str, 'theta');
    str.vgams = vgam;
    str = param2sse(str);
    for i = 1:nparma
        beta(i) = vgam(bind(i));
    end
    beta(nparma+1:end) = Lparm; %insert parameters for the Cholesky factor
end


%compute covariance matrix of residuals using the last parameters in beta
L = zeros(neqs);
l = 0;
betam = [1, beta(nparma+1:end)];
for i = 1:neqs
    cont = neqs - i + 1;
    ind = l + 1:l + cont;
    L(i:end, i) = betam(ind);
    l = l + cont;
end
sigmar = L * L';
str.sigmar2 = sigmar;

p = neqs;
[nx, r] = size(x);
T = str.Fs;
Z = str.Hs;
HH = str.Ks;
HH = HH * L;
GG = L;
[nalpha, mf] = size(T);
if r > 0
    Bs = str.Bs;
    Ds = str.Ds;
    inc = 1;
    kro = str.kro;
    [z, rx1] = varmafil(x, T, Z, Bs, Ds, kro, inc);
    %subtract filtered inputs from output
    yf = y - z;
    %estimate regression parameters and initial state
    if nY > p
        YY = Y;
    elseif nY == p
        YY = repmat(Y, nobs, 1);
    else
        YY = Y;
    end
    Yx1 = [YY, rx1]; %design matrix
    yfvec = vec(yf'); %put the observations in a vector
    bet = Yx1 \ yfvec;
    m1 = bet(nbeta+1:end); %estimated initial state
    yfvec = yfvec - rx1 * m1; %subtract estimated initial state
    yx = reshape(yfvec, neqs, nobs);
    y = yx'; %transform from vector to matrix
else
    Bs = [];
    Ds = [];
end

%set up regression matrices
W = [];
X = Y;
%set up initial state
Sigma = dlyapsq(T, HH);
Sigma = Sigma' * Sigma;
ins = Sigma;
i = [nalpha, 0, 0, 0];

%We have used y_t = V_t + U_t, where V_t is the exogenous part. Store in a
%structure the model for U_t, the regression estimator and its mse if Y is
%not empty, the state space at the end of filtering and its mse, and the
%state space model corresponding to the exogenous part, V_t.
stx.X = X;
stx.Z = Z;
stx.G = GG;
stx.W = W;
stx.T = T;
stx.H = HH;
stx.ins = ins;
stx.i = i;

if (chb == 0) || (chb == 1)
    %likelihood evaluation
    [e, f, hb, Mb, A, P] = scakfle2(y, X, Z, GG, W, T, HH, ins, i, chb);
    ff = e * f;
elseif (chb == 2)
    if ~isempty(X)
        [Xt, Pt, hb, Mb, initf, recrs] = scakff(y, X, Z, GG, W, T, HH, ins, i);
    else
        hb = [];
        [Xt, Pt, recrs] = scakfff(y, X, Z, GG, W, T, HH, ins, i, hb);
    end
elseif (chb == 3)
    if ~isempty(X)
        [Xt, Pt, hb] = scakff(y, X, Z, GG, W, T, HH, ins, i);
    else
        hb = [];
    end
    [Xt, Pt, recrs] = scakfff(y, X, Z, GG, W, T, HH, ins, i, hb);
end

stx.hb = hb;
stx.Mb = Mb;
stx.A = A;
stx.P = P;
stx.B = Bs;
stx.D = Ds;
stx.m1 = m1;
stx.kro = str.kro;
