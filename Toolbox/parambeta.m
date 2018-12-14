function [DA, ferror] = parambeta(beta)
%
% Given an (s x r) matrix of rank r, this function parameterizes beta by
% finding r linearly independent rows. It returns an (s x r+1) matrix of
% the form [beta Idx], where beta is parameterized and Idx is an index such
% that Idx(i) = 0 if the i-th row is linearly independent and Idx(i) = 1 if
% the i-th row is linearly dependent.
%
% Inputs  :    betap : an (s x r) matrix
%  Output :     DA   : an (s x r+1) matrix such that DA=[beta Idx]
%
%
% Copyright (c) 21 July 2014 by Victor Gomez
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
DAr = beta;
Indxr = [];
if ~isempty(beta)
    [s, r] = size(beta);
    betap = beta';
    %Indx gives the linearly independent rows of DAr (those equal to zero)
    [Q, R, Indxr, ferror] = housref(betap);
    Br = [];
    for i = 1:s
        if Indxr(i) == 0
            Br = [Br; beta(i, :)];
        end
    end
    Br = pinv(Br);
    DAr = DAr * Br;
    Indxr = Indxr';
end
DA = [DAr, Indxr];
