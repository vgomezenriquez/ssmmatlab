function str = armaxe2sse(str)
% PURPOSE: given a VARMAX model in echelon form, it computes the state
% space echelon form
%---------------------------------------------------
% USAGE: str = armaxe2sse(str)
% where:    str    = a structure containing the structure of the VARMAX in
%                    echelon form
%---------------------------------------------------
% RETURNS: str = a structure containing the previous structure plus
%                the state space echelon form matrices according to the
%                following state space model:
%
%  alpha_{t+1} = Fs*alpha_{t} + Bs*x_t{t} + Ks*a_{t}
%      y_{t}   = Hs*alpha_{t} + Ds*x_{t}  + a_{t}
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

phis = str.phis;
thetas = str.thetas;
gammas = str.gammas;
npar = str.npar;
s = str.s;
kro = str.kro;
nx = str.m;

nlag = max(kro);

% phist=str.phist; thetast=str.thetast;
% gammast=str.gammast;

% compute state space matrices corresponding to VARMAX echelon form
% First compute VARMAX polynomials with \Phi_0=I_s=\Theta_0. They are necessary
% for the state space form
% phist=phis; thetast=thetas; gammast=gammas;
% phi0=phis(:,:,1); inv=0;
% if any(any(phi0 - eye(s)))
%  inv=1;
% end
% if inv == 1
%  if nx > 0
%   gammast(:,:,1)=phi0\gammast(:,:,1);
%  end
%  for i=2:nlag+1
%   phist(:,:,i)=phi0\phis(:,:,i);
%   thetast(:,:,i)=phi0\thetas(:,:,i);
%   if (nx > 0)
%    gammast(:,:,i)=phi0\gammas(:,:,i);
%   end
%  end
% end
% phist(:,:,1)=eye(s); thetast(:,:,1)=eye(s);
[phist, thetast, gammast, ierror] = armaxe2armax(phis, thetas, gammas);
str.phist = phist;
str.thetast = thetast;
str.gammast = gammast;
% then compute state space matrices

% F=str.F; H=str.H;
sk = sum(kro);
Fs = zeros(sk);
Hs = zeros(s, sk);

ch = 0;
for i = 1:s
    if kro(i) > 0
        Hs(i, ch+1) = 1;
        Fs(ch+1:ch+kro(i)-1, ch+2:ch+kro(i)) = eye(kro(i)-1);
    end
    chb = 0;
    for j = 1:s
        if kro(i) > 0
            Fs(ch+kro(i), chb+1:chb+npar(i, j)) = ...
                -phis(i, j, kro(i)+1:-1:kro(i)-npar(i, j)+2);
            chb = chb + kro(j);
        else
            Hs(i, chb+1:chb+npar(i, j)) = -phis(i, j, 1:-1:-npar(i, j)+2);
            chb = chb + kro(j);
        end
    end
    ch = ch + kro(i);
end


str.Fs = Fs;
str.Hs = Hs;

% K=str.K; Ks=K; D=str.D; Ds=D; B=str.B; Bs=B;
Ks = zeros(sk, s);
Ds = zeros(s, nx);
Bs = zeros(sk, nx);
if nlag > 0
    Psi = thetast(:, :, 2) - phist(:, :, 2);
    if nx > 0
        Ds = gammast(:, :, 1);
        Psig = gammast(:, :, 2) - phist(:, :, 2) * Ds;
    end
else
    Psi = [];
    Psig = [];
end
for i = 2:nlag
    A = [];
    for j = i:-1:2
        A = [A, phist(:, :, j)];
    end
    Psix = thetast(:, :, i+1) - phist(:, :, i+1) - A * Psi;
    Psi = [Psi; Psix];
end
ch = 0;
for i = 1:s
    kri = kro(i);
    A = [];
    for k = i:s:s * kri
        A = [A; Psi(k, :)];
    end
    if kri > 0
        Ks(ch+1:ch+kri, :) = A;
        ch = ch + kri;
    end
end
str.Ks = Ks;
if nx > 0
    for i = 2:nlag
        A = [];
        for j = i:-1:2
            A = [A, phist(:, :, j)];
        end
        Psix = gammast(:, :, i+1) - phist(:, :, i+1) * Ds - A * Psig;
        Psig = [Psig; Psix];
    end
    ch = 0;
    for i = 1:s
        kri = kro(i);
        A = [];
        for k = i:s:s * kri
            A = [A; Psig(k, :)];
        end
        if kri > 0
            Bs(ch+1:ch+kri, :) = A;
            ch = ch + kri;
        end
    end
end
str.Bs = Bs;
str.Ds = Ds;
