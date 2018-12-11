function [compf, ferror] = dtanbut(D, Thetap, Thetas, Di, Thetac, Lambda)
%
% This function obtains the tangent Butterworth filter corresponding to the
% parameters d (differencing degree) and xc (frequency where gain is 1/2).
% See "The Use of Bitterworth Filters for Trend and Cycle Estimation in
% Economic Time Series", Gómez, V. (2001), Journal of Business and Economic
% Statistics, 19, 365-373.
% The filter model is
%
%           z_t = s_t + n_t,
%   Alpha(B)s_t = num(B) b_t,    Var(b_t)=1
%
% where Alpha(z) = (1 - z)^Di, num(z) = (1 + z)^Di, and n_t and b_t are
% independent white noises.
%
% The filter numerator is (1/sa). The filter denominator is den. Thus,
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
%     Thetap  : a number, design frequency Thetap divided by pi. It can be
%               empty.
%     Thetas  : a number, design frequency Thetas divided by pi. It can be
%               empty.
%     Di      : a number, the exponent in Alpha(z). It can be empty.
%     Thetac  : a number, the frequency, divided by pi, of gain .5 in the
%               But. tangent filter. It can be empty.
%     Lambda  : a number, the signal to noise ratio (sigma^2_n/sigma^2_b)
%               in the But. tangent filter. It can be empty.
% Note: The usual specification is D, Thetap and Thetas (Di, Thetac and
%       Lambda empty). Alternatively, the user can enter Di and Thetac (D,
%       Thetap, Thetas and Lambda empty) or Di and Lambda (D, Thetap,
%       Thetas and Thetac empty).
%
% Output parameters: compf, a structure containing the following fields
%     .num    : a polynomial of degree Di, (1 + z)^Di (filter numerator)
%     .den    : a polynomial of degree Di  (filter denominator)
%     .sa     : a positive number
%     .Alpha  : a polynomial, (1 - z)^Di
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
if (nargin == 5)
    %specification Di and Thetac
    if isempty(Di) || isempty(Thetac)
        disp('Di and Thetac must be entered in dtanbut when there are')
        disp('five arguments')
        ferror = 1;
        return
    end
    dd = double(Di);
    Thetac = Thetac * pi;
    Lambda = 1.D0 / ((tan(Thetac/2.D0))^Di);
elseif (nargin == 6)
    %specification Di and Lambda
    if isempty(Di) || isempty(Lambda)
        disp('Di and Lambda must be entered in dtanbut when there are')
        disp('six arguments')
        ferror = 2;
        return
    end
    dd = double(Di);
    Thetac = 2.D0 * atan(1.D0/exp(log(Lambda)/(dd + dd)));
    Lambda = sqrt(Lambda);
else
    %specification D, Thetap and Thetas
    if isempty(D) || isempty(Thetap) || isempty(Thetas)
        disp('D , Thetap and Thetas2 must be entered in dtanbut when')
        disp('there are less than five arguments')
        ferror = 3;
        return
    end
    sum = 0.D0;
    for i = 1:2
        sum = sum + log(D(i)) - log(1.D0-D(i));
    end
    Thetap = Thetap * pi;
    Thetas = Thetas * pi;
    deno = 2.D0 * (log(tan(Thetap/2.D0)) - log(tan(Thetas/2.D0)));
    dd = sum / deno;
    Di = round(dd);
    dd = double(Di);
    sum = (log(D(1)) - log(1.D0-D(1))) / (2.D0 * dd);
    sum = exp(sum);
    sum = (2.D0 * tan(Thetap/2.D0)) / sum;
    sum = 2.D0 * atan(sum/2.D0);
    Thetac = sum;
    Lambda = 1.D0 / ((tan(Thetac/2.D0))^Di);
end
nterm2 = floor(Di/2);
if (mod(Di, 2) == 0)
    nterm1 = 0;
else
    nterm1 = 1;
end
den = 1.D0;
sa = Lambda;
%form moving average polynomial for the aggregate series
b = tan(Thetac/2.D0);
%form den vector
for i = 1:nterm2
    alp0 = b^2 + 1.D0 + b * sqrt(2.D0*(1.D0 - cos(double(2*i-1)*pi/dd)));
    alp2 = b^2 + 1.D0 - b * sqrt(2.D0*(1.D0 - cos(double(2*i-1)*pi/dd)));
    alp1 = 2.D0 * (b^2 - 1.D0);
    sa = sa * abs(alp0);
    den = conv([alp2 / alp0, alp1 / alp0, 1.D0], den);
end
if (nterm1 == 1)
    alp0 = b + 1.D0;
    alp1 = b - 1.D0;
    sa = sa * abs(alp0);
    den = conv([alp1 / alp0, 1.D0], den);
end
num = 1.D0;
delta = [1., 1.];
for i = 1:Di
    num = conv(delta, num);
end
delta = [-1., 1.];
Alpha(1) = 1.D0;
for i = 1:Di
    Alpha = conv(delta, Alpha);
end

compf.num = num;
compf.den = den;
compf.sa = sa;
compf.Alpha = Alpha;
compf.Di = Di;
compf.Thetac = Thetac;
compf.Lambda = Lambda;
