function ser = lbma
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
% --    THIS IS A SIMULATED MA(1) SERIES USED BY G.M. LJUNG AND
% --    G.E.P. BOX (1976) IN THE PAPER "MAXIMUM LIKELIHOOD ESTIMATION
% --    IN THE ARMA MODEL" (TR 476, STATISTICS DEPARTMENT, UNIVERSITY
% --    OF WISCONSIN, ALSO APPEARED IN BIOMETRIKA
% --    THE MODEL IS
% --          Z(T) = (1-0.9B)A(T)     WITH VAR(A(T))=1.0
% --    THERE ARE 75 OBSERVATIONS

yor = load(fullfile('data', 'lbma.dat'));
bg_year = 1970;
bg_per = 1;
freq = 1;
Y = [];
ser.yor = yor;
ser.Y = Y;
ser.bg_year = bg_year;
ser.bg_per = bg_per;
ser.freq = freq;
