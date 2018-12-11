function inft = minft(fid, fh, wd, nd, scale)
%************************************************************************
% This function puts the printing options into structure inft
%
%          INPUTS:
%            .fid: the device on which the table will be written
%             .fh: flag for header and years
%             .wd: format width
%             .nd: number of decimal points
%          .scale: =1 scale data if necessary
%                  =0 do not scale data
%
%          OUPTUT:
%            inft: structure with printing options given by the inputs
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

inft.fid = fid;
inft.fh = fh;
inft.wd = wd;
inft.nd = nd;
inft.scale = scale;
