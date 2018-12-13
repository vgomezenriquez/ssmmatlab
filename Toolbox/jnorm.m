function jnr = jnorm(x, p, q)
%      given an n-column vector x and a signature matrix J=diag(I_p,-I_q), this function
%      calculates the J-norm of x, that is, the quantity sqrt(x'*J*x), assuming that
%      x'*J*x > 0.
%
%---------------------------------------------------
% USAGE: jnr = jnorm(x,p,q)
% where:    x = an n-column vector
%           p,q= integers such that J = diag(I_p,-I-q) is a signature
%           matrix
%---------------------------------------------------
% RETURNS: the J-norm of x
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
n = length(x);
if n ~= p + q
    error('n nonequal to p + q in jnorm')
end
jnr = 0.d0;
if p > 0
    jnr = jnr + mnorm(x(1:p));
end
if q > 0
    jnr = jnr - mnorm(x(p+1:p+q));
end
if jnr < 0
    %  jnr
    error('j-norm negative in jnorm')
else
    jnr = sqrt(jnr);
end
