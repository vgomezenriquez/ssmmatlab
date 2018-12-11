function [compf, ferror] = dbptanbut(D, Omegap1, Omegap2, Omegas2, Di, Thetac, ...
    Lambda)
%
% This function obtains the band-pass filter based on the Butterworth
% tangent filter corresponding to the parameters D(1), D(2), Omegap1,
% Omegap2, Omegas2. See "The Use of Bitterworth Filters for Trend and Cycle
% Estimation in Economic Time Series", Gómez, V. (2001), Journal of
% Business and Economic Statistics, 19, 365-373.
% The filter model is
%
%           z_t = s_t + n_t,
%   Alpha(B)s_t = num(B)b_t,    Var(b_t)=1
%
% where Alpha(z) = (1 - 2*Alph*z + z^2)^Di, num(z) = (1 - z^2)^Di, and n_t
% and b_t are independent white noises.
%
% The filter numerator is num*(1/sa). The filter denominator is den. Thus,
%
% H(z) = (1/sa)*(num(z)/den(z))
%
% The other filter is
%
% G(z) = (sqrt(Lambda)/sa)*(Alpha(z)/den(z))
%
% Input parameters:
%     D       : a (1 x 2) array containing the design tolerances D1 and D2.
%               It can be empty.
%     Omegap1 : a number, design frequency Omegap1 divided by pi. Required.
%     Omegap2 : a number, design frequency Omegap2 divided by pi. Required.
%     Omegas2 : a number, design frequency Omegas2 divided by pi. It can be
%               empty
%     Di      : a number, the exponent in Alpha(z) and num(z). It can be
%               empty.
%     Thetac  : a number, the frequency, divided by pi, of gain .5 in the
%               But. tan. filter. It can be empty.
%     Lambda  : a number, the signal to noise ratio (sigma^2_n/sigma^2_b)
%               in the But. tangent filter. It can be empty.
% Note: The usual specification is D, Omegap1, Omegap2 and Omegas2 (Di,
%       Thetac and Lambda empty). Alternatively, the user can enter
%       Omegap1, Omegap2, Di and Thetac (D, Omegas2 and Lambda empty) or
%       Omegap1, Omegap2, Di and Lambda (D, Omegas2 and Thetac empty).
%
% Output parameters: compf, a structure containing the following fields
%     .nterm1 : number of factors of degree 2 in den polynomial below
%     .nterm2 : number of factors of degree 4 in den polynomial below
%     .term1  : factors of degree 2 in den polynomial below
%     .term2  : factors of degree 4 in den polynomial below
%     .num    : a polynomial of degree 2*Di, (1 - z^2)^Di
%     .den    : a polynomial of degree 2*Di
%     .sa     : a positive number (sigmaa)
%     .Alpha  : a polynomial of degree 2*Di, (1 - 2*Alph*z + z^2)^Di
%     .Alph   : a number, in Alpha(z) = (1 - 2*Alph*z + z^2)^Di
%     .Di     : a positive integer
%     .Thetac : a number, the frequency of gain .5 in the But. tan. filter
%     .Lambda : a positive number, the square root of the noise to signal
%               ratio (sigma^2_n/sigma^2_b) in the But. tangent filter.
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

ferror = 0;
if isempty(Omegap1) || isempty(Omegap2)
    disp('Omegap1 and Omegap2 must be entered in dbptanbut')
    ferror = 1;
    return
end
Omegap1 = Omegap1 * pi;
Omegap2 = Omegap2 * pi;

if (nargin == 6)
    %specification Omegap1, Omegap2, Di and Thetac
    if isempty(Di) || isempty(Thetac)
        disp('Di and Thetac must be entered in dbptanbut when there are')
        disp('six arguments')
        ferror = 2;
        return
    end
    Thetac = Thetac * pi;
    dd = double(Di);
    Lambda = 1.D0 / ((tan(Thetac/2.D0))^Di);
elseif (nargin == 7)
    %specification Omegap1, Omegap2, Di and Lambda
    if isempty(Di) || isempty(Lambda)
        disp('Di and Lambda must be entered in dbptanbut when there are')
        disp('seven arguments')
        ferror = 3;
        return
    end
    dd = double(Di);
    Thetac = 2.D0 * atan(1.D0/exp(log(Lambda)/(dd + dd)));
    Lambda = sqrt(Lambda);
else
    %specification D, Omegap1, Omegap2 and Omegas2
    if isempty(D) || isempty(Omegas2)
        disp('D and Omegas2 must be entered in dbptanbut when there are less')
        disp('than six arguments')
        return
    end
    Omegas2 = Omegas2 * pi;
    % transformation in the frequency domain: THETA = OMEGA-OMEGA_P1
    thetap = Omegap2 - Omegap1;
    thetas = Omegas2 - Omegap1;
    sum = 0.D0;
    for i = 1:2
        sum = sum + log(D(i)) - log(1.D0-D(i));
    end
    deno = 2.D0 * (log(tan(thetap/2.D0)) - log(tan(thetas/2.D0)));
    dd = sum / deno;
    Di = round(dd);
    dd = double(Di);
    sum = (log(D(1)) - log(1.D0-D(1))) / (2.D0 * dd);
    sum = exp(sum);
    sum = (2.D0 * tan(thetap/2.D0)) / sum;
    sum = 2.D0 * atan(sum/2.D0);
    Thetac = sum;
    Lambda = 1.D0 / ((tan(Thetac/2.D0))^Di);
end
% now transformation in the time domain
Alph = cos((Omegap1 + Omegap2)/2.D0) / cos((Omegap2 - Omegap1)/2.D0);
nterm2 = floor(Di/2);
if (mod(Di, 2) == 0)
    nterm1 = 0;
else
    nterm1 = 1;
end
compf.nterm1 = nterm1;
compf.nterm2 = nterm2;
compf.term1 = [];
compf.term2 = [];
sa = Lambda;
% set up moving average polynomial for the aggregate series
b = tan(Thetac/2.D0);
% den vector
den = 1.D0;
for i = 1:nterm2
    alp0 = b^2 + 1.D0 + b * sqrt(2.D0*(1.D0 - cos(double(2*i-1)*pi/dd)));
    alp2 = b^2 + 1.D0 - b * sqrt(2.D0*(1.D0 - cos(double(2*i-1)*pi/dd)));
    alp1 = 2.D0 * (b^2 - 1.D0);
    sa = sa * abs(alp0);
    delta(1) = 1.D0;
    alp1 = alp1 / alp0;
    alp2 = alp2 / alp0;
    % do transformation
    delta(2) = Alph * (alp1 - 2.D0);
    delta(3) = Alph * Alph * (1.D0 - alp1 + alp2) - alp1;
    delta(4) = Alph * (alp1 - 2.D0 * alp2);
    delta(5) = alp2;
    delta = fliplr(delta);
    compf.term2{i} = delta;
    den = conv(delta, den);
    clear delta;
end
if (nterm1 == 1)
    alp0 = b + 1.D0;
    alp1 = b - 1.D0;
    sa = sa * abs(alp0);
    alp1 = alp1 / alp0;
    pol = [-alp1, Alph * (alp1 - 1.D0), 1.D0];
    compf.term1 = pol;
    den = conv(pol, den);
end
delta = [-1.D0, 0.D0, 1.D0];
num(1) = 1.D0;
for i = 1:Di
    num = conv(delta, num);
end
delta = [1.D0, -2.D0 * Alph, 1.D0];
Alpha(1) = 1.D0;
for i = 1:Di
    Alpha = conv(delta, Alpha);
end

compf.num = num;
compf.den = den;
compf.sa = sa;
compf.Alpha = Alpha;
compf.Alph = Alph;
compf.Di = Di;
compf.Thetac = Thetac;
compf.Lambda = Lambda;
