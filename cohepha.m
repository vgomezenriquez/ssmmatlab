function [co, ph, ga] = cohepha(cxy, qxy, fxx, fyy)
%
%        This function computes the coherence, gain and phase
%        (see Granger's book)
%
%     INPUTS:
%------------
%       cxy : cospectrum
%       qxy : quadrature spectrum
%       fxx : (smoothed) periodogram of x
%       fyy : (smoothed) periodogram of y
%
%
%    OUTPUTS:
%------------
%    co     : coherence
%    ph     : phase angle
%    ga     : gain
%
% Copyright (c) 21 July 2003 by Victor Gomez
% Ministerio de Hacienda, Direccion Gral. de Presupuestos,
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
alpha = cxy.^2 + qxy.^2;
co = alpha ./ (fxx .* fyy); % coherence
ga = ((alpha).^(.5)) ./ fyy; % gain
l = length(cxy);
ph = zeros(l, 1);
for i = 1:l
    ph(i) = atan2(qxy(i), cxy(i)); % phase angle
end
