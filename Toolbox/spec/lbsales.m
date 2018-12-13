function ser = lbsales
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
% --    MONTHLY SALES OF A COMPANY (JANUARY 1965 - MAY 1971)
% --    THIS SERIES WAS ANALYZED BY CHATFIELD AND PROTHERO (1973)
% --    IN THE PAPER "BOX-JENKINS SEASONAL FORECASTING:  PROBLEM
% --    IN A CASE STUDY" (JRSS A 136: 295-336).  G.M. LJUNG AND
% --    G.E.P. BOX (1976) USED IT AS AN EXAMPLE IN THE PAPER
% --    "MAXIMUM LIKELIHOOD ESTIMATION IN THE ARMA MODEL" (TR 476,
% --    STATISTICS DEPARTMENT, UNIVERSITY OF WISCONSIN, ALSO APPEARED
% --    IN BIOMETRIKA.)
% --    THERE ARE 77 OBSERVATIONS

yor = load(fullfile('data', 'lbsales.dat'));
bg_year = 1965;
bg_per = 1;
freq = 12;
Y = [];
ser.yor = yor;
ser.Y = Y;
ser.bg_year = bg_year;
ser.bg_per = bg_per;
ser.freq = freq;
