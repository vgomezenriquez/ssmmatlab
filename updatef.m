function [f,fc]=updatef(ff,ffc)
%
% function to update the vector whose nonlinear sum of squares is
% minimized. It is stored in the form f=(f^(1/n))*(2^(fc/n) ) to avoid
% underflow and overflow.
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
     while abs(ff) >= 1.
       ff=ff*0.0625;
       ffc=ffc+4;
     end
     while abs(ff) <0.0625
       ff=ff*16.;
       ffc=ffc-4;
     end
     f=ff;
     fc=ffc;
