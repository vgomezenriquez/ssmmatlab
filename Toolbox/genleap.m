function Y = genleap(Iy, Im, N, Mq)

% this function generates the leap year variable.  It works until 2100.
%
% this function generates the leap year variable.  It works until 2100.
%
% input variables       Iy      : the initial year
%                       Im      : the initial period
%                       N       : the length of the desired vector
%                       Mq      : the series frequency (=12 for monthly,
%                                                       =4 for quarterly)
%
% output variables      Y       : N x 1 vector containing the leap year
%                                 variable%
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

YY = trade(Iy, Im, N, 2, 0, [], Mq);
Y = YY(:, 2);