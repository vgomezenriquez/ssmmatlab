function plotspcd(outa)
%
% function to plot the spectra of the canonical decomposition of an ARIMA
% model previously identified with function arimaestos.
%
% phi(B)*phi_s(B^s)*(delta*delta_s*y_t -mu) =
% th(B)*th_s(B^s)*a_t
%
% Copyright (c) 21 July 2015 by Victor Gomez
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
thrc = outa.compcd.ptnum; % trend-cycle numerator
phir = outa.compcd.ptden; % trend-cycle denominator
%   outa.compcd.ptnur;    % number of nonstationary roots in phir
sigma2r = outa.compcd.ptvar; % variance of the trend-cycle innovations (*)
thsc = outa.compcd.stnum; % seasonal numerator
phis = outa.compcd.stden; % seasonal denominator
%   outa.compcd.stnur; % number of nonstationary roots in phis
sigma2s = outa.compcd.stvar; % variance of the seasonal innovations (*)
thtc = outa.compcd.rt; % transitory outa.compcdonent (MA term)
sigma2t = outa.compcd.rtvar; % variance of the transitory outa.compcdonent innovations (*)
%(*) in units of the series model innovations

%plot canonical trend spectrum
specgraph(thrc, phir, sigma2r)
title('canonical trend spectrum')
pause
%plot canonical seasonal spectrum
specgraph(thsc, phis, sigma2s)
title('canonical seasonal spectrum')
pause
%plot canonical transitory spectrum
if ~isempty(thtc)
    phit = 1.;
    specgraph(thtc, phit, sigma2t)
    title('canonical transitory spectrum')
    pause
end
close all
