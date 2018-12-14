function [Tsp, Hsp, Zsp, ferror] = cascade(den, Alpha, phip, thp, sigma)
%
%        This function obtains a cascade implementation of the state space
%        form corresponding to the product of filters in the ARMA model
%
%        y_t = [Alpha(z)/den(z)]*[thp(z)/phip(z)] a_t
%
%        where a_t is (0,sigma^2). In each factor the degree of the
%        numerator has to be equal to that of the denominator.
%        If z_t = [thp(z)/phip(z)] a_t has a state space form
%
%        x^z_{t+1} = T_z x^z_t + H_z*sigma*e_t                   (1)
%           z_t    = Z_z x^z_t + sigma*e_t,
%
%        where Var(e_t)=1,  and if y_t = [Alpha(z)/den(z)] z_t has a state
%        space form
%
%        x^y_{t} = T_y x^y_{t-1} + H_y z_t                       (2)
%           y_t  = Z_y x^y_t,
%
%        then the following state space form for y_t is a cascade
%        implementation
%
%        x_{t}   = [ T_y   H_y*Z_z]x_{t-1}  + [ H_y*sigma ]e_t
%                  [  0    T_z    ]           [ H_z*sigma ]
%          y_t   = [ Z_y     0    ]x_t ,
%
%        where x_t = [x^y_t'  x^z_{t+1}']'.
%        Here, T_z is a square matrix with dimension equal to the degree of
%        phip(z), T_y is a square matrix with dimension equal to the degree
%        of den(z) plus one, and the representations (1) and (2) are
%        Akaike's representations. Note that (1) is minimal while (2) is
%        not.
%
%        Input parameters:
%        den    : a (1 x nbp) array
%        Alpha  : a (1 x nal) array
%        phip   : a (1 x np+1) array
%        thp    : a (1 x nt+1) array
%        sigma  : a positive constant
%
%        Output parameters:
%        Tsp    : an (nalpha x nalpha) matrix
%        Hsp    : an (nalpha x nepsilon) matrix
%        Zsp    : an (p x nalpha) matrix
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

[Tp, Hp, Zp, ferror1] = akaikessm1(phip, thp);
if ferror1 > 0
    ferror = ferror1;
    return
end

[Ts, Hs, Zs, ferror2] = akaikessm2(den, Alpha);
if ferror2 > 0
    ferror = ferror2;
    return
end

np = length(phip) - 1;
nbp = length(den);
nsp = np + nbp;

Tsp = zeros(nsp);
Tsp(1:nbp, 1:nbp) = Ts;
Tsp(nbp+1:end, nbp+1:end) = Tp;

Hsp = zeros(nsp, 1);
Hsp(1:nbp) = Hs;
Hsp(nbp+1:end) = Hp;

Tsp(1:nbp, nbp+1) = Hsp(1:nbp);
Hsp = Hsp * sigma;

Zsp = zeros(1, nsp);
Zsp(1:nbp) = Zs;
