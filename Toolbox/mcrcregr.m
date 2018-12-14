function [D, nr, yd, DA, ferror] = mcrcregr(y, x)
%
% This function applies the CRC criterion to the multivariate y series to
% obtain the number of unit roots. This criterion is the generalization to
% the multivariate case of the one described for univariate time series in
% the paper "A Strongly Consistent Criterion to Decide Between I(1) and
% I(0) Processes Based on Different Convergence Rates" by Víctor Gómez,
% (2013), Communications in Statistics - Simulation and Computation, 42,
% pp. 1848-1864.
% Maximum regular differencing considered is one.
%
% Inputs: y: matrix containing the output series
%         x: matrix containing the input series
%  Output: D: an (ny x ny x seas+1) matrix 'differencing' matrix polynomial
%         nr: an integer, number of unit roots
%         yd= matrix containing the 'differenced' series
%         DA= matrix of the form [DAr Indxr], where DAr is the regular
%             parameterization of the differencing polynomial, and Indxr is
%             an index vector to identify the l.i. rows of DAr.
%         ferror= a flag for errors
%
% Copyright (c) 21 July 2014 by Victor Gomez
% Ministerio de Hacienda y A.P., Direccion Gral. de Presupuestos,
% Subdireccion Gral. de Analisis y P.E.,
% Alberto Alcocer 2, 1-P, D-34, 28046, Madrid, SPAIN.
% Phone : +34-915835439
% E-mail: vgomez@sepg.minhap.es
%
% The author assumes no responsibility for errors or damage resulting from the use
% of this code. Usage of this code in applications and/or alterations of it should
% be referenced. This code may be redistributed if nothing has been added or
% removed and no money is charged. Positive or negative feedback would be appreciated.
%

ferror = 0;
yd = [];
nr = 0;
D = [];
DA = [];

if nargin < 2
    x = [];
end

[n, s] = size(y);
[nx, m] = size(x);
if ~isempty(x)
    if (nx ~= n)
        ferror = 1;
        disp('mcrcregr: nobs in x-matrix not the same as y-matrix');
        return
    end;
end
%         maxr: maximum regular differencing order considered
maxr = 1;

%correction for x
if ~isempty(x)
    beta = mulols(y, x);
    y = y - x * beta;
    x = [];
end

%differenced series
yd = y;

%maximum orders
maxrs = maxr * s;

if (s == 1)
    %there is only one series
    [nr1, nr] = crcregr(y, maxr);
    if (nr > 0)
        yd = diferm(yd, 1);
    end
    return
else
    %preliminary univariate unit root analysis
    nru = zeros(1, s);
    for i = 1:s
        [nr11, nr1] = crcregr(y(:, i), maxr);
        %  plot(y(:,i)),i,nr1, pause
        nru(i) = nru(i) + nr1;
    end
    nrsu = sum(nru);
    if nrsu == 0
        return
    end
end

can = .13; %.11
seas = 0;
alphamax = .499;
alpha1 = min(alphamax, .5-1/n);
alpha2 = min(alphamax, .5-1/(n^.60));
% alpha1
% alpha2
% pause

DAr = [];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% First step: AR(p)(1)_seas model
%
hm = max(n^(-alphamax), n^(-alpha1));
hm1 = 1 - hm;

p = 6;
P = 0;
hr3 = 1; %perform only the first two stages of the HR method
[str, ferror] = estvarmaxpqrPQR(y, x, seas, [p, 0, 0], [P, 0, 0], hr3);
%copy in a the residuals estimated with the long VARX of the original
%series
a = str.residv;


%'differencing' polynomial
Dz(:, :, 1) = eye(s);
Ds(:, :, 1) = eye(s);
Dr(:, :, 1) = eye(s);

%obtain matrices containing the unit roots
%if p > 1, we should estimate unit roots using the companion form, not
%using Pi(1)-eye(s).
Pir = zeros(s);
%matrices Pi(1) and dPi(1) of regular part
phir = str.phis(:, :, 1:p+1);
% Pi1=-eye(s); Pid1=zeros(s);
for i = 1:p
    phirip1 = phir(:, :, i+1);
    Pir = Pir + phirip1;
end
% e=eig(Pir),abs(e),hm1


%obtain unit roots
[nr] = eurpi(Pir, hm1);
if nr > 0
    r = s - nr; %rank of Pi1r
    t = s - r; %rank of alphaor'*Pid1r*betaor
    %  r,t,nr
    %estimate the model again with the rank r imposed
    ir = 1;
    [Pi1r] = evarma11rurimp(y, x, a, ir, r);
    Pid1r = -Pi1r - eye(s);
    %obtain differencing polynomial
    [Dr1, Fd, ferror] = mdifpol(r, t, Pi1r, Pid1r);
    F = -Dr1(:, :, 2);
    DAr = Fd;
    Dr(:, :, 1) = Dr1(:, :, 1);
    Dr(:, :, 2) = -F;
    %  [yd,ferror]=mdifer(yd,F,1,y);
    %  [zd,ferror]=mdifer(zd,F,1);
    Dz = pmatmul(Dr, Dz);
end

%obtain ''differenced series''
if (nr > 0)
    [yd, ferror] = dtimesy(Dz, y);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Second step: VARMA(1,1) model
%
hms = max(n^(-alphamax), n^(-alpha2));
hm1s = 1 - hms;
mflag = 0;
if (nr < maxrs)
    mflag = 1;
end
while mflag == 1
    mflag = 0;
    %estimate VARMA(1,1) model
    hr3 = 0; %run third stage of HR method
    [str, ferror] = estvarmaxpqrPQR(yd, x, seas, [1, 1, 0], [0, 0, 0], hr3);
    if isfield(str, 'phis3')
        if str.nonst3 == 0
            Pir = str.phis3(:, :, 2);
            Thr = str.thetas3(:, :, 2);
        else
            Pir = str.phis(:, :, 2);
            Thr = str.thetas(:, :, 2);
        end
        %    abs(eig(Pir)),abs(eig(Thr)),hm1s
    else
        Pir = str.phis(:, :, 2);
        Thr = str.thetas(:, :, 2);
        %    abs(eig(Pir)),abs(eig(Thr)),hm1s
    end
    %obtain unit root in regular part.
    [nrp, ferror] = findurpir(Pir, -Thr, hm1s, can);
    if (nrsu > 0) && (nr == 0) && (nrp == 0)
        nrp = nrp + 1;
    end
    if (nrp > 0) && (nr < maxrs) && (nrp < maxrs)
        nr = nr + 1;
        mflag = 1;
        %obtain corrected and 'differenced' series
        r = s - nr; %rank of Pi1
        t = s - r; %rank of alphaor'*Pid1 *betaor
        %    r,t,nr
        %estimate the model again with the rank r imposed
        ir = 1;
        [Pi1r] = evarma11rurimp(y, x, a, ir, r);
        Pid1r = -Pi1r - eye(s);
        %obtain differencing polynomial
        [Dr1, Fd, ferror] = mdifpol(r, t, Pi1r, Pid1r);
        F = -Dr1(:, :, 2);
        DAr = Fd;
        Dr(:, :, 1) = Dr1(:, :, 1);
        Dr(:, :, 2) = -F;
        Dz = pmatmul(Dr, Ds);
        [yd, ferror] = dtimesy(Dz, y);
    end
end
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% 'differencing' polynomial
D = Dz;

%parameterization
[DA, ferror] = parambeta(DAr);
