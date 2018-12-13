function g = gammln(x)
%
% auxiliary function called in functions gacf and gser
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


cof = [76.18009172947146; -86.50532032941677; 24.01409824083091; -1.231739572450155; ...
    .1208650973866179e-2; -.5395239384953e-5];
stp = 2.5066282746310005;
y = x;
tmp = x + 5.5;
tmp = (x + 0.5) * log(tmp) - tmp;
ser = 1.000000000190015;
for j = 1:6
    y = y + 1;
    ser = ser + cof(j) / y;
end
g = tmp + log(stp*ser/x);
