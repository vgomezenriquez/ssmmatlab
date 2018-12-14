function [str, ierror] = enfstap(str)
% PURPOSE: enforces stationarity in a VARMAX model in echelon form
%          using polynomial methods
%---------------------------------------------------
% USAGE:  str = enfstap(str)
% where:    str is a structure containing all the information about the
%           model
%---------------------------------------------------
% RETURNS: the same structure with the stationary model
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
ierror = 0;
H = str.H;
Sigma = str.Sigma;
kro = str.kro;

% Kb=zeros(size(K)); Hb=zeros(size(H)); Sigmab=zeros(size(Sigma));
% phi=str.phi; theta=str.theta; npar=str.npar;
% phib=zeros(size(phi));
% vgams=str.vgams;

% Fp=F-K*H;
Sigmai = pinv(Sigma); %matrices for the inverse process
[s, mh] = size(H);
%
% % obtain arma process corresponding to inverse process
% [thetatp,phitp,ierror]=ss2lcvarmax(-K',Fp',H',G);

phi = str.phi;
theta = str.theta;
%pass from left to right MFD
[thetat, phit, krolr, ierrorp] = pleft2rightcmfd(theta, phi, sum(kro)+1);
phitp = pmattrans(phit); %transpose phit
[np, mp, pp1] = size(phitp);

% %Kronecker indices of the inverse process
% [phiei,thetaei,kroi,ierrori] = pecheform(transpmat(thetat),phitp);
% kroi
% [Hm,Fm,Gm,Jm,ierror] = qarmax2ss2(transpmat(thetat),phitp);
% [Hro,Fro,Kro,u,krob,ierror]=pkroindrev(Hm,Fm,Gm);krob


% %check that the eigenvalues of thetatp and phitp are correct
% [Hm,Fm,Gm,Jm,ierror] = qarmax2ss2(thetat,phit);
% FF=abs(eig(F))
% FFGHm=abs(eig(Fm-Gm*Hm))
% FFp=abs(eig(Fp))
% FFm=abs(eig(Fm))

% obtain invertible inverse process
% left divide phitp(z) by phitp(0)
A = phitp(:, :, 1);
if any(any(A-eye(s)))
    phitp(:, :, 1) = eye(s);
    for i = 2:pp1
        phitp(:, :, i) = A \ phitp(:, :, i);
    end
end

phis = phitp;
if (sum(sum(Sigmai-eye(s))) ~= 0)
    for i = 1:pp1
        phis(:, :, i) = phitp(:, :, i) * Sigmai;
    end
end
% phis
% phitp
G = pmmulbf(phis, phitp); %generating function of MA part
Gm = G(:, :, pp1:end);
[Sigmab, Theta, ierrorpm, iter, normdif] = pmspectfac(Gm, 20, 1e-3);

%adjustment for A different from zero because Theta(0)=I
phib(:, :, 1) = A;
if any(any(A-eye(s)))
    for i = 2:pp1
        phib(:, :, i) = A * Theta(:, :, i-1);
    end
else
    for i = 2:pp1
        phib(:, :, i) = Theta(:, :, i-1);
    end
end


% %pass invertible arma model to state space form
% [Hm,Fm,Gm,Jm,ierror] = qarmax2ss2(thetat,transpmat(phib));
% %check that the eigenvalues are correct
% FFGHm=abs(eig(Fm-Gm*Hm)),FFm=abs(eig(Fm))
% Hb=Gm'; Kb=-Hm'; Fb=Fm'+Kb*Hb;


% %obtain minimal state space form
% [Hro,Fro,Kro,u,kroo,ierror]=pkroindrev(Hm,Fm,Gm);
% kroo,ierror
% %obtain stationary state space model from invertible inverse process
% Hb=Kro'; Kb=-Hro'; Fb=Fro'+Kb*Hb;
% % aa=eig(F)
% % bb=eig(Fb)
% % cc=eig(Fp)
% % dd=eig(Fb-Kb*Hb)
% % % check whether Kronecker indices have changed
% % [Hro,Fro,Kro,u,kroo,ierror]=pkroindrev(Hb,Fb,Kb);
% % kroo,ierror


% %obtain echelon form. Here we impose the original Kronecker indices, that
% %may not coincide with the Kronecker indices that correspond to the new
% %system. For this reason, the procedure is not stable.
% [nh,mh]=size(Hb); Gb=eye(nh);
% [phib,thetab,Hb,Fb,Kb,ierror]=echeform(Hb,Fb,Kb,Gb,kro);
% % Es=abs(eig(Fb))
% % Gs=abs(eig(Fb-Kb*Hb))
%
%
% str.Kb=Kb; str.Sigmab=pinv(Sigmab); str.phib=phib; str.Fb=Fb; str.thetab=thetab;
% str.Hb=Hb;
% % str.Kb
% % str.Sigmab
% % str.phib


%alternative polynomial method to obtain the echelon form. The problem is
%that the Kronecker indices may change and function pecheform needs the exact
%Kronecker indices. The procedure works at least for the generic neighborhood.


% pass from right to left MFD
phibt = pmattrans(phib);
[thetar, phir, krorl, ierrorp] = pright2leftcmfd(thetat, phibt, sum(kro));
if ierrorp == 3
    nn = floor(sum(kro)*2);
    [thetar, phir, krorl, ierrorp] = pright2leftcmfd(thetat, phibt, nn);
end


%assuming the Kronecker indices have not changed
% [phie,thetae,kro,ierror] = pecheform(phir,thetar,[],kro);
%compute Kronecker indices
[phieb, thetaeb, krob, ierrorp] = pecheform(phir, thetar);
if any(krob-kro)
    phie = phi;
    thetae = theta;
    Sigmab = Sigmai;
    disp('Kronecker indices changed in enfstap from ')
    disp(kro)
    disp('to:')
    disp(krob)
    ierror = 1;
else
    phie = phieb;
    thetae = thetaeb;
end

stre = str;
stre.phis = phie;
stre.thetas = thetae;
stre.gammas = str.gamma;
stre.npar = str.npar;
stre.s = str.s;
stre.kro = str.kro;
stre.m = str.m;
stre = armaxe2sse(stre);
Kb = stre.Ks;
phib = stre.phis;
Fb = stre.Fs;
Hb = stre.Hs;

str.Kb = Kb;
str.Sigmab = pinv(Sigmab);
str.phib = phie;
str.Fb = Fb;
str.thetab = thetae;
str.Hb = Hb;
% str.Kb
% str.Sigmab
% % str.phib
