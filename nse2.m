function [ct2, str] = nse2(y, residv, x, tsig2, str)
% PURPOSE:
% Eliminates nonsignificant parameter after second step of HR method.
% The Kronecker indices should be preserved. This means
% that in each row of the VARMAX model the maximum degree is the Kronecker
% index. For example, if s=1, p or q is equal to the Kronecker index. If
% s=2, in the first row, p_1 or q_1 is equal to k_1, in the second row, p_2
% or q_2 is equal to k_2, etc.
% Put nsig2=1 and choose tsig2 for insignificant t-value
%---------------------------------------------------
% USAGE: [ct2,str] = nse2(y,residv,x,tsig2,str)
% where:    str    = a structure containing the structure of the VARMAX
% model
%---------------------------------------------------
% RETURNS:
%         ct2 = the number of parameters eliminated
%          str = a structure containing the inverted model
%---------------------------------------------------
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ct2 = 0;
nsig2 = 1;
while nsig2
    [mp, minc, mint] = nselimhr2(y, x, str);
    if mint >= tsig2
        nsig2 = 0;
    else
        %   mp, minc, mint
        if strcmp(mp, 'ph')
            str.phi(minc(1), minc(2), minc(3)) = 0;
            str.nparm = str.nparm - 1;
            ct2 = ct2 + 1;
        end
        if strcmp(mp, 'th')
            str.theta(minc(1), minc(2), minc(3)) = 0;
            str.nparm = str.nparm - 1;
            ct2 = ct2 + 1;
            if minc(3) == 1
                str.phi(minc(1), minc(2), minc(3)) = 0;
            end
        end
        if strcmp(mp, 'ga')
            str.gamma(minc(1), minc(2), minc(3)) = 0;
            str.nparm = str.nparm - 1;
            ct2 = ct2 + 1;
        end
        str = mhanris2(y, residv, x, str);
    end
end
str = vecparwr(str); %fill in vgam with new model parameters
