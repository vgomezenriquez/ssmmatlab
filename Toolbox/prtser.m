function prtser(fid, fname, yor, y, ny, datei, inft, lam)
%*************************************************************************
% This function prints time series into text file
%
%    INPUTS:
%      fid : file identifier, needed for writing the output into text file
%    fname : name of the series
%      yor : original time series
%        y : time series used in the estimation etc.
%       ny : number of observations
%    datei : calendar structure
%     inft : structure containing printing options
%      .fh :  flag for header and years
%      .wd : format width
%      .nd : number of decimal points
%   .scale = 1 : scale data if necessary
%          = 0 : do not scale data
%      lam = 0 : compute logs of y
%          = 1 : do not compute logs
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
%**************************************************************************

fprintf(fid, 'Series:%s\n', fname);
fprintf(fid, 'Number of observations:%5i\n', ny);
my = min(abs(yor));
if my <= 1d-9, inft.scale = 0;
end %do not scale if numbers are too small
tabla(yor, datei, inft);
fprintf(fid, '\n');
if lam == 0 %print series in logs
    fprintf(fid, 'Logs of the series:\n');
    fprintf(fid, 'Number of observations:%5i\n', ny);
    my = min(abs(y));
    if my <= 1d-9, inft.scale = 0;
    end %do not scale if numbers are too small
    tabla(y, datei, inft);
    fprintf(fid, '\n');
end
