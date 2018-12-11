function str = param2armaxe(str)
% PURPOSE: given a vector of Hannan-Rissanen estimates, it computes the VARMAX
% echelon form
%---------------------------------------------------
% USAGE: str = param2armaxe(str)
% where:    str    = a structure containing the vector of second step
%                    estimates
%---------------------------------------------------
% RETURNS: str = a structure containing the previous structure plus
%                the matrices of the VARMAX echelon form
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
phis = str.phi;
thetas = str.theta;
gammas = str.gamma;
s = str.s;
kro = str.kro;
vgams = str.vgams;
nx = str.m;

vgamtv = str.vgamtv;
phitv = NaN(size(phis));
thetatv = NaN(size(thetas));
gammatv = NaN(size(gammas));

neqs2 = s * s;
nlag = max(kro);
for i = 2:nlag + 1
    A = [];
    Atv = [];
    for j = 1:s
        A = [A, vgams((i - 2)*neqs2+(j - 1)*s+1:(i - 2)*neqs2+j*s)];
        Atv = [Atv, vgamtv((i - 2)*neqs2+(j - 1)*s+1:(i - 2)*neqs2+j*s)];
    end
    phis(:, :, i) = A;
    phitv(:, :, i) = Atv;
end
ii = (i - 1) * neqs2;
for i = 1:nlag + 1
    A = [];
    Atv = [];
    for j = 1:s
        A = [A, vgams(ii+(i - 1)*neqs2+(j - 1)*s+1:ii+(i - 1)*neqs2+j*s)];
        Atv = [Atv, vgamtv(ii+(i - 1)*neqs2+(j - 1)*s+1:ii+(i - 1)*neqs2+j*s)];
    end
    thetas(:, :, i) = A;
    thetatv(:, :, i) = Atv;
end
% for i=1:nlag+1
%  Atv=[];
%  for j=1:s
%   Atv=[Atv vgamtv(ii+(i-1)*neqs2+(j-1)*s+1:ii+(i-1)*neqs2+j*s)];
%  end
%  thetatv(:,:,i)=Atv;
% end
phis(:, :, 1) = thetas(:, :, 1);
phitv(:, :, 1) = thetatv(:, :, 1);
ii = ii + i * neqs2;
if (nx > 0)
    for i = 1:nlag + 1
        A = [];
        Atv = [];
        for j = 1:nx
            A = [A, vgams(ii+(i - 1)*s*nx+(j - 1)*s+1:ii+(i - 1)*s*nx+j*s)];
            Atv = [Atv, vgamtv(ii+(i - 1)*s*nx+(j - 1)*s+1:ii+(i - 1)*s*nx+j*s)];
        end
        gammas(:, :, i) = A;
        gammatv(:, :, i) = Atv;
    end
end
mu = vgams(end-s+1:end); %parameters for the mean
mutv = vgamtv(end-s+1:end);

if nx == 0
    Phi = zeros(s);
    for i = 1:nlag + 1
        Phi = Phi + phis(:, :, i);
    end
    musers = Phi \ mu;
    str.musers = musers; % mean value obtained from the constant
end

str.phis = phis;
str.phitv = phitv;
str.thetas = thetas;
str.thetatv = thetatv;
str.gammas = gammas;
str.gammatv = gammatv;
str.mu = mu;
str.mutv = mutv;
