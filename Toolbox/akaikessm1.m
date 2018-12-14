function [T, H, Z, ferror] = akaikessm1(phip, thp)
%
%        This function obtains Akaike's state space representation of
%        minimal dimension corresponding to the ARMA model
%
%        y_t = [thp(z)/phip(z)] a_t
%
%        where a_t is (0,1). The model is
%
%        x_{t+1} = T x_t + H*a_t
%           y_t    = Z x_t + a_t,
%
%        where
%            [0 1 0  ... ....   0]        [ Psi_1    ]
%            [0 0 1  ... ....   0]        [ Psi_2    ]
%     T =    [ ...   ...  ...    ],   H = [ ...      ],
%            [0 0 0  ... ....   1]        [ Psi_{r-1}]
%            [-phi_r ...   -phi_1]        [ Psi_r    ]
%
%     Z =    [1 0 0  ...  ... 0],
%
%     r = degree(phip) = degree(thp)  and phi^{1}(z)*theta(z) =
%     Psi_0 + Psi_1*z + Psi_2*z^2+ ...
%
%        Input parameters:
%        phip   : a (1 x np+1) array
%        thp    : a (1 x nt+1) array
%
%        Output parameters:
%        T    : a (r x r) matrix
%        H    : a (r x 1) matrix
%        Z    : a (1 x r) matrix
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
ferror = 0;
np = length(phip) - 1;
nt = length(thp) - 1;
if np ~= nt
    disp('dimension mismatch in akaikessm1')
    ferror = 1;
    return
end

npm1 = np - 1;
T = zeros(np);
T(1:npm1, 2:np) = eye(npm1);
T(np, :) = -phip(1:end-1);
psip = poldiv(fliplr(thp), fliplr(phip), np);
H = psip(2:end)';
Z = zeros(1, np);
Z(1) = 1.;
