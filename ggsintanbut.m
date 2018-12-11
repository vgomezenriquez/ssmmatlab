function ggsintanbut(D, Thetap, Thetas, d, thc)
%
% This function plots the gain function corresponding to a sine or
% tangent Butterworth filter.
% Input parameters:
%     D       : a (1 x 2) array containing the design tolerances D1 and D2.
%               It can be empty.
%     Thetap  : a number, design frequency Thetap divided by pi. It can be
%               empty
%     Thetas  : a number, design frequency Thetas divided by pi. It can be
%               empty
%     d       : a number, the exponent in Alpha(z) and num(z). Required
%     thc     : a number, the frequency for gain .5, divided by pi.
%               Required
% Note: The parameters d, and thc define the filter. If D, Thetap and
%       Thetas are entered by the user, the program will draw the
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

if isempty(d) || isempty(thc)
    disp('d and thc should be specified in ggsintanbut')
    return
end

thc = thc * pi;
t = 0:.01:pi;
a = tan(t/2).^(2 * d);
b = tan(thc/2)^(2 * d);
a = a / b;
c = ones(size(t));
c1 = (c ./ (1 + a)).^2; %tangent
a = sin(t/2).^(2 * d);
b = sin(thc/2)^(2 * d);
a = a / b;
c2 = (c ./ (1 + a)).^2; %sine

%
% line for gain equal to one half
%
c3 = ones(size(t)) .* .5;

if ~isempty(D) && ~isempty(Thetap) && ~isempty(Thetas)
    Thetap = Thetap * pi;
    Thetas = Thetas * pi;
    %
    % specification lines
    %
    tt = 0:.01:1 - D(1);
    n = length(tt);
    c4 = ones(size(tt));
    c5 = ones(size(tt));
    for j = 1:n,
        c4(j) = Thetap;
        c5(j) = tt(j);
    end
    s = 0:.01:Thetap;
    c6 = ones(size(s));
    c7 = c6 - D(1) * c6;
    r = Thetas:.01:pi;
    c8 = ones(size(r));
    c8 = D(2) * c8;
    plot(t, c1, '-', t, c2, '--', t, c3, '-.', c4, c5, '-', s, c6, '-', s, c7, '-', r, c8, '-')
else
    plot(t, c1, '-', t, c2, '--', t, c3, '-.')
end
title('Butterworth filters');
axis([0, pi, 0, 1.25]);
legend('BFT ', 'BFS ');
