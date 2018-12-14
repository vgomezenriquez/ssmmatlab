function [olsres] = OLSres(out)
%
%
% This function obtains the OLS residuals after having used function
% arimaestos, arimaestni or arimaestwi
%
% input arguments:
% out: a structure, the output of function arimaestos, arimaestni or
% arinamestwi
%
% output arguments:
% res: a vector containing the OLS residuals
%
% Copyright (c) 21 July 2003 by Victor Gomez
% Ministerio de Hacienda, Direccion Gral. de Presupuestos,
% Subdireccion Gral. de Analisis y P.E.,
% Alberto Alcocer 2, 1-P, D-34, 28046, Madrid, SPAIN.
% Phone : +34-915835439
% E-mail: vgomez@sepg.minhac.es
%
% The author assumes no responsibility for errors or damage resulting from the use
% of this code. Usage of this code in applications and/or alterations of it should
% be referenced. This code may be redistributed if nothing has been added or
% removed and no money is charged. Positive or negative feedback would be appreciated.
%
if isfield(out.model, 'hb')
    %   Ee=out.model.matsis.OLSMatrix;
    %   beta=out.model.hb;
    %   olsres=Ee(:,end) - Ee(:,1:end-1)*beta;
    olsres = out.model.matsis.olsres;
else
    olsres = [];
end
