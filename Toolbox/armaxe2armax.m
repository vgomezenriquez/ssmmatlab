function [phi, theta, gamma, ierror] = armaxe2armax(phie, thetae, gammae)
%
% This function computes VARMAX polynomials with \Phi_0=I_s=\Theta_0.
%---------------------------------------------------
% USAGE: [phi,theta,gamma] = armaxe2armax(phie,thetae,gammae)
% where:    phie   = a k x k polynomial matrix with phi(0) nonsingular
%           thetae = a k x k polynomial matrix
%          gammae  = a k x m polynomial matrix
%---------------------------------------------------
% RETURNS:
%           phi  = the AR polynomial matrix
%          theta = the MA polynomial matrix
%          gamma = the input polynomial matrix
%        ierror =1, dimension mismatch in phi and theta
%               =0, there are no errors on input
%---------------------------------------------------
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

phi = [];
theta = [];
gamma = [];
ierror = 0;
[p1, p2, p3] = size(phie);
[t1, t2, t3] = size(thetae);
[s1, s2, s3] = size(gammae);
if (p1 ~= p2) | (p1 ~= t1)
    disp('wrong dimensions of phi or theta in armaxe2armax');
    ierror = 1;
    return
end


phi = phie;
theta = thetae;
gamma = gammae;
phi0 = phie(:, :, 1);
inv = 0;
if any(any(phi0-eye(p1)))
    inv = 1;
end
if inv == 1
    if s2 > 0
        gamma(:, :, 1) = phi0 \ gamma(:, :, 1);
    end
    for i = 2:p3
        phi(:, :, i) = phi0 \ phie(:, :, i);
        theta(:, :, i) = phi0 \ thetae(:, :, i);
        if (s2 > 0)
            gamma(:, :, i) = phi0 \ gammae(:, :, i);
        end
    end
end
phi(:, :, 1) = eye(p1);
theta(:, :, 1) = eye(t1);
