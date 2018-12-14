function ggbptanbut(D, omp1, omp2, oms2, d, alph, lambda)
%
% This function plots the gain function corresponding to a
% band-pass filter based on the tangent Butterworth filter.
% Input parameters:
%     D       : a (1 x 2) array containing the design tolerances D1 and D2.
%               It can be empty.
%     ompp1   : a number, design frequency Omegap1 divided by pi. It can be
%               empty
%     omp2    : a number, design frequency Omegap2 divided by pi. It can be
%               empty
%     oms2    : a number, design frequency Omegas2 divided by pi. It can be
%               empty
%     d       : a number, the exponent in Alpha(z) and num(z). Required
%     alph    : alph parameter in 1 - 2*alph*z + z^2. Required
%     lambda  : a number, the square root of the signal to noise ratio
%              (sigma^2_n/sigma^2_b) in the But. tangent filter. Required
% Note: The parameters d, alph and lambda define the filter. If D, omp1,
%       omp2 and oms2 are entered by the user, the program will draw the
%       toleration lines.
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

if isempty(d) || isempty(alph) || isempty(lambda)
    disp('d, alph and lambda should be specified in ggbptanbut')
    return
end

lambda = lambda^2;
t = 0:.01:pi;
a = sin(t).^(2 * d);
b = (cos(t) - alph).^(2 * d);
c1 = (a ./ (a + lambda * b)).^2;
%
% line for gain equal to one half
%
c2 = ones(size(t)) .* .5;

if ~isempty(D) && ~isempty(omp1) && ~isempty(omp2) && ~isempty(oms2)
    omp1 = omp1 * pi;
    omp2 = omp2 * pi;
    oms2 = oms2 * pi;
    %
    % specification lines
    %
    s1 = omp1:.01:omp2;
    c3 = ones(size(s1));
    c4 = c3 * (1 - D(1));
    s2 = 0:.01:1 - D(1);
    n = length(s2);
    c5 = ones(size(s2));
    c6 = ones(size(s2));
    c7 = ones(size(s2));
    for j = 1:n,
        c5(j) = omp1;
        c6(j) = omp2;
        c7(j) = s2(j);
    end
    s4 = oms2:.01:pi;
    c8 = ones(size(s4));
    c8 = D(2) * c8;
    plot(t, c1, '-', t, c2, '--', s1, c3, '-', s1, c4, '-', c5, c7, '-', c6, c7, '-', s4, c8, '-')
else
    plot(t, c1, '-', t, c2, '--')
end
title('Band-pass filter');
axis([0, pi, 0, 1.25]);
