function str = param2sse(str)
% PURPOSE: given a vector of Hannan-Rissanen estimates, it computes the
% state space echelon form
%---------------------------------------------------
% USAGE: str = param2sse(str)
% where:    str    = a structure containing the vector of second step
%                    estimates
%---------------------------------------------------
% RETURNS: str = a structure containing the previous structure plus
%                the matrices of the VARMAX echelon form
%---------------------------------------------------
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

% compute polynomial matrices corresponding to the second step estimates
str = param2armaxe(str);

str = armaxe2sse(str);
