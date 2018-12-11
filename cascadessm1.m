function [Fsp, Gsp, Hsp, Jsp, ferror] = cascadessm1(Fp, Gp, Hp, Fs, Gs, Hs, Js)
%
%        This function obtains a cascade implementation of the state space
%        form corresponding to the product of VARMA filters
%
%        y_t = As(z)*Ap(z) a_t
%
%        where a_t is (0,Sigma).
%        If z_t = Ap(z) a_t has a state space form
%
%        x^p_{t+1} = F_p x^p_t + G_p*a_t                   (1)
%           z_t    = H_p x^p_t + a_t,
%
%        and if y_t = As(z) z_t has a state space form
%
%        x^s_{t+1} = F_s x^s_t + G_s*z_t                       (2)
%           y_t    = H_s x^s_t + J_s*z_t,
%
%        then the following state space form for y_t is a cascade
%        implementation
%
%        x_{t+1} = [ F_s   G_s*H_p]x_t  + [ G_s ]a_t
%                  [  0      F_p  ]       [ G_p ]
%          y_t   = [ H_s   J_s*H_p]x_t  + [ J_s ]a_t,
%
%        where x_t = [x^s_t'  x^p_t']'.
%        The representations (1) and (2) are
%        Akaike's representations. Note that (1) and (2) are minimal.
%
%        Input parameters:
%        Fp    : an (np x np) matrix
%        Gp    : an (np x n) matrix
%        Hp    : an (n x np) matrix
%        Fs    : an (ns x ns) matrix
%        Gs    : an (ns x m) matrix
%        Hs    : an (m x ns) matrix
%        Js    : an (m x ns) matrix
%
%        OuFput parameters:
%        Fsp    : an (nalpha x nalpha) matrix, nalpha = np + ns.
%        Gsp    : an (nalpha x 1) matrix
%        Hsp    : an (m x nalpha) matrix
%        Jsp    : an (m x ns) matrix
%
% Copyright (c) 23 March 2018 by Victor Gomez
% Ministerio de Hacienda, Direccion Gral. de Presupuestos,
% % Subdireccion Gral. de Analisis y P.E.,
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
Fsp = [];
Gsp = [];
Hsp = [];
Jsp = [];

[np, mp] = size(Fp);
[nGp, mGp] = size(Gp);
[nHp, mHp] = size(Hp);
if (np ~= mp) || (nGp ~= np) || (mHp ~= np)
    disp('dimension mismath in cascadessm2')
    ferror = 1;
    return
end
[ns, ms] = size(Fs);
[nGs, mGs] = size(Gs);
[nHs, mHs] = size(Hs);
[nJs, mJs] = size(Js);
if (ns ~= ms) || (nGs ~= ns) || (mGs ~= mJs) || (nHs ~= nJs) || ...
        (mHs ~= ns) || (mGs ~= nHp) || (mGp ~= mJs)
    disp('dimension mismath in cascadessm2')
    ferror = 1;
    return
end

nsp = np + ns;

Fsp = zeros(nsp);
Fsp(1:ns, 1:ns) = Fs;
Fsp(ns+1:end, ns+1:end) = Fp;

Gsp = zeros(nsp, mGp);
Gsp(1:ns, :) = Gs;
Gsp(ns+1:end, :) = Gp;

Fsp(1:ns, ns+1:nsp) = Gs * Hp;

Hsp = zeros(nHs, nsp);
Hsp(:, 1:mHs) = Hs;
Hsp(:, mHs+1:end) = Js * Hp;
Jsp = Js;
