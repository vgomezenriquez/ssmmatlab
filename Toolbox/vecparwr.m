function str = vecparwr(str)
% PURPOSE: given a structure of a VARMAX model in echelon form, forms a
% vector with the nonrestricted parameters and restricted parameters
%---------------------------------------------------
% USAGE: str = vecparwr(str)
% where:    str    = a structure containing the structure of the VARMAX 
%                    model in echelon form
%---------------------------------------------------
% RETURNS: str = a structure containing the previous structure plus 
%                a vector with the the nonrestricted parameters and 
%                restricted parameters
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
phi=str.phi; theta=str.theta; gamma=str.gamma; neqs=str.s; 
nx=str.m; kro=str.kro; nlag=max(kro); 
vgam=[];
for i=2:nlag+1
 vgam=[vgam; vec(phi(:,:,i))];
end
for i=1:nlag+1
 vgam=[vgam; vec(theta(:,:,i))];
end
if (nx > 0)
 for i=1:nlag+1
  vgam=[vgam; vec(gamma(:,:,i))];
 end
end
vgam=[vgam; NaN(neqs,1)];   %parameters for the mean
%  vgam
str.vgam=vgam;
