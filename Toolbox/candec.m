function [comp, ierrcandec] = candec(phir, phis, thr, ths, phirst, s, dr, ds, sconp)
%************************************************************************
% PURPOSE: This function performs the canonical decomposition of the ARIMA
%          model
%          phir(B)*phis(B^s)y_t = thr(B)*ths(B^s)*a_t
%The innovation variance, sigma2a, is assumed to be unity.
%---------------------------------------------------
% USAGE: [comp,ierrcandec] = candec(phir,phis,thr,ths,phirst,s,dr,ds,sconp)
%
% Inputs: phir : a polynomial containing the regular AR part
%         phis : a polynomial containing the seasonal AR part
%         thr  : a polynomial containing the regular MA part
%         ths  : a polynomial containing the seasonal MA part
%       phirst : a polynomial containing the stationary regular AR part
%            s : a positive integer, the number of seasons
%           dr : a positive integer, the number of regular differences
%           dr : a positive integer, the number of seasonal differences
%        sconp : a positive number, standard deviation of the series model
%                innovations
% Note: all of the previous polynomials are expressed in the same variable.
% The polynomials are given by an array like [ a_n, ... a_1, a_0], where
% the polynomial is a_0 + a_1*z + ... + a_n*z^n.
%
%  Output: comp, a structure containing the following fields
%   .ptnum=thrc;     % trend-cycle numerator
%   .ptden=phir;     % trend-cycle denominator
%   .ptnur=dr+ds;    % number of nonstationary roots in phir
%   .ptvar=sigma2r;  % variance of the trend-cycle innovations (*)
%   .stnum=thsc;     % seasonal numerator
%   .stden=phis;     % seasonal denominator
%   .stnur=(s-1)*ds; % number of nonstationary roots in phis
%   .stvar=sigma2s;  % variance of the seasonal innovations (*)
%   .rt=thtc;        % transitory component (MA term)
%   .rtvar=sigma2t;  % variance of the transitory component innovations (*)
%   .itvar=sigma2i;  % variance of the irregular component (*)
%   .sigmaa=sconp;   % standard deviation of the series model innovations
%   .phi=phirst;     % stationary AR trend polynomial
%(*) in units of the series model innovations
%    ierrcandec : flag for errors
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
%*************************************************************************

thrc = [];
thsc = [];
sigma2r = [];
sigma2s = [];
thtc = [];
sigma2t = [];
ierrcandec = 0;

%form MA part (model numerator)
th = conv(thr, ths);

%compute pseudospectrum

%numerator
nnum = length(th) - 1;
num = acgf(1., th, nnum);
%denominator
ndenr = length(phir) - 1;
denr = acgf(1., phir, ndenr);
ndens = length(phis) - 1;
dens = acgf(1., phis, ndens);

%change S_n = z^n + z^(-n) variables to powers of U = z + z^(-1)
%numerator
[unum, ierr1] = sn2u(num); %total numerator as a polynomial in U
%denominator
[udenr, ierr2] = sn2u(denr);
[udens, ierr3] = sn2u(dens);
uden = conv(udenr, udens); %total denominator as a polynomial in U

%divide num by den if necessary
dunum = length(unum) - 1;
dudenr = length(udenr) - 1;
dudens = length(udens) - 1;
duden = dudenr + dudens;
if dunum >= duden
    [uq, ur] = deconv(unum, uden);
    dur = dunum - 1;
    ur = ur(2:end);
else
    uq = 0.;
    ur = unum;
    dur = dunum;
end

%perform partial fraction expansion
%polynomials udenrm and udensm are assumed to be coprime
%first, obtain gcd and polynomials unr and uns such that
%unr(x)*udenrm(x) + uns(x)*udensm(x) = gcd(x)
udenrm = zeros(1, 1, dudenr+1);
udenrm(:, :, 1:end) = fliplr(udenr);
udensm = zeros(1, 1, dudens+1);
udensm(:, :, 1:end) = fliplr(udens);
[G, U, Dr, Nr, ierrglcd] = glcd(udenrm, udensm);
ng = length(G); %G = gcd should be a constant
A = G(:);
na = length(A);
if na > 1
    normaa = norm(A);
    %  tol=double(na)*eps(normaa);  %Matlab's default
    %  tol=tol*10.;
    tol = 1d-14;
    for i = 2:ng
        if abs(G(:, :, i)) > tol
            disp(G(:, :, i)), disp(tol)
            disp('udenrm and udensm are not coprime in candec')
            return
        end
    end
    clear G;
    G = A(1);
end
[nu, mu, pu] = size(U);
unr = zeros(1, pu);
uns = zeros(1, pu);
unr(1:end) = U(1, 1, :);
uns(1:end) = U(2, 1, :);
unr = fliplr(unr) ./ G;
uns = fliplr(uns) ./ G; %polynomials are normalized
%the following sum should be equal to gcd (one)
% [C, ierrsumpol] = sumpol(conv(unr,udenr),conv(uns,udens))
%then, obtain the numerators of the fractions in the expansion
a = conv(unr, ur);
b = conv(uns, ur);
%in the following, rs and rr are the numerators of the regular and seasonal
%spectra. We clean them from trailing zeros.
tol = 1d-12;
[q, rs] = deconv(a, udens);
rs = cleanpol(rs, tol); %rs is the numerator of udens
[qq, rr] = deconv(b, udenr);
rr = cleanpol(rr, tol); %rr is the numerator of udenr
%the following sum should be equal to zero
% [C, ierrsumpol] = sumpol(q,qq)

%minimization using the MATLAB function fminbnd
%the variable is U=z + z^(-1). If z=exp(-ix), then U=2*cos(x).
%Since x \in [0, \pi], U \in [-2, 2].

[xr, fvalr1, exitflagr] = mifmin(rr, udenr, -2., 2.); %fvalr is the regular minimum
fvalr2 = polyval(rr, -2.) / polyval(udenr, -2.); %try also at the border (pi)
if fvalr2 < fvalr1
    fvalr = fvalr2; %xr=-2.;
else
    fvalr = fvalr1;
end
[xs, fvals1, exitflags] = mifmin(rs, udens, -2., 2.); %fvals is the seasonal minimum
fvals2 = polyval(rs, 2.) / polyval(udens, 2.); %try also at the border (0)
if fvals2 < fvals1
    fvals = fvals2; %xs=2.;
else
    fvals = fvals1;
end

%subtract minimum values from spectra to make them canonical
[rrmin, irrmin] = sumpol(rr, -udenr.*fvalr);
[rsmin, irsmin] = sumpol(rs, -udens.*fvals);

%add these minimum values to the irregular spectrum (uq)
[irregspct, irspct] = sumpol(uq, fvalr);
[irregspct, irspct] = sumpol(irregspct, fvals); %irregular spectrum

%solve the two covariance factorization problems to obtain the MA parts
%of the canonical components
[thrc, sigma2r, ierrpu2mar] = pu2ma(rrmin);
[thsc, sigma2s, ierrpu2mas] = pu2ma(rsmin);
ierrcandec = ierrpu2mar + ierrpu2mas;

%if irregular is not white noise, make it canonical too
ni = length(irregspct);
if ni > 1
    [xi, fvali1, exitflagr] = mifmin(irregspct, 1., -2., 2.);
    fvali2 = polyval(irregspct, 2.);
    if fvali2 < fvali1
        fvali1 = fvali2; %xi=2.;
    end
    fvali2 = polyval(irregspct, -2.);
    if fvali2 < fvali1
        fvali = fvali2; %xi=-2.;
    else
        fvali = fvali1;
    end
    [irmin, irrmin] = sumpol(irregspct, -fvali);
    [thtc, sigma2t, ierrpu2mai] = pu2ma(irmin);
    sigma2i = fvali;
else
    sigma2i = irregspct;
end

comp.ptnum = thrc;
comp.ptden = phir;
comp.ptnur = dr + ds;
comp.ptvar = sigma2r;
comp.stnum = thsc;
comp.stden = phis;
comp.stnur = (s - 1) * ds;
comp.stvar = sigma2s;
comp.rt = thtc;
comp.rtvar = sigma2t;
comp.itvar = sigma2i;
comp.sigmaa = sconp;
comp.phi = phirst;
