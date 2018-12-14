function vgam = enfstab(stre, pol)
%
% This function multiplies the parameters of a matrix polynomial by some
% factor of the form .95^n until it becomes stable (all roots outside the
% unit circle).
%
% Input arguments:
% stre: a structure, containing model information. In particular,
% stre.vgams: vector containing the parameters of the polynomial matrix
%       in the form vec(phi_1),...,vec(phi_r), vec(theta_0),vec(theta_1),
%       ..., vec(theta_r), vec(gamma_0),...,vec(gamma_r)
% pol: 'phi' or 'theta'
%
% Output arguments:
% vgam: the parameter vector containing corresponding to the stable
%       polynomial matrix
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

vgam = stre.vgams;
kro = stre.kro;
s = stre.s;
% str = matechelon(kro,s,m);  %create structure
% str = vecparwr(str);        %create vector of parameters
r = max(kro);
s2 = s * s;
s2r = s2 * r;

%First option: transform Phi(z) into Phi(lambda*z)
if pol == 'phi  '
    l = s2r;
    vgamp = vgam(1:l);
    maeval = max(abs(eig(stre.Fs)));
elseif pol == 'theta'
    f = s2r + s2 + 1;
    l = s2r + s2 + s2r;
    vgamp = vgam(f:l);
    maeval = max(abs(eig(stre.Fs-stre.Ks*stre.Hs)));
end
epsilon = 0.015;
c = 1. + epsilon;
lambda = 1 / (c * maeval);

cont = 0;
li = lambda;
for i = 1:r
    for j = cont + 1:cont + s2
        vgamp(j) = vgamp(j) * li;
    end
    cont = cont + s2;
    li = li * lambda;
end
if pol == 'phi  '
    l = s2r;
    vgam(1:l) = vgamp;
elseif pol == 'theta'
    f = s2r + s2 + 1;
    l = s2r + s2 + s2r;
    vgam(f:l) = vgamp;
end

% %Second option: Hannan and Kavalieris, (1984), "Multivariate Linear Time
% %Series Models".
% if pol == 'phi  '
%  l=s2r; vgamp=vgam(1:l);
% elseif pol == 'theta'
%  f=s2r+s2+1; l=s2r+s2+s2r;
%  vgamp=vgam(f:l);
% end
% const=.995; fac=const;
% for i=1:500
% %  i,fac,pol
%  vgamp=vgamp.*fac;
%  if pol == 'phi  '
%   l=s2r; vgam(1:l)=vgamp;
%  elseif pol == 'theta'
%   f=s2r+s2+1; l=s2r+s2+s2r;
%   vgam(f:l)=vgamp;
%  end
%  str.vgams=vgam; str.vgamtv=str.vgam;
%  str = param2sse(str);
%  if pol == 'phi  '
%   iar=chkstainv(str.Fs);
%   if iar == 0, break, end
%  elseif pol == 'theta'
%   ima=chkstainv(str.Fs-str.Ks*str.Hs);
%   if ima == 0, break, end
%  end
%  fac=fac*const;
% end
