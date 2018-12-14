function [phir, phis, thr, ths, phirst] = arima2rspol(phi, Phi, th, Th, freq, dr, ds)
%*************************************************************************
%
% This function returns the trend, seasonal and stationary polynomials,
% phir, thr, phis, ths and phirst, corresponding to the canonical
% decomposition of an ARIMA model,
%
%        y_t = [thr(B)/phir(B)]b_t + [ths(B)/phis(B)]c_t + u_t,
%
% where phirst, if it exists, is a stationary factor of phir.
%
% The original ARIMA model is given by the regular and seasonal
% polynomilas, phi, Phi, th and Th, in matrix polynomial format, such that
%
%        phi(B)*Phi(B^s)y_t = th(B)*Th(B^s)*a_t
%
% For example, phi(:,:,1)=1; phi(:,:,2)=-1, etc.
%
%*************************************************************************
%     INPUTS :
%        phi : regular autoregressive polynomial
%        Phi : seasonal autoregressive polynomial
%         th : regular moving average polynomial
%         Th : seasonal moving average polynomial
%       freq : frequency of the data
%         dr : number of regular differences
%         ds : number of seasonal differences
%
%    OUTPUTS :
%       phir : autoregressive trend polynomial
%       phis : autoregressive seasonal polynomial
%        thr : moving average trend polynomial
%        ths : moving average seasonal polynomial
%      phirst: stationary autoregressive trend polynomial (factor of phir)
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
%*************************************************************************

s = freq;

if ~isempty(phi)
    phi = squeeze(phi)';
    % Find positive and complex roots and assign them to phir
    % The negative roots are assigned to phis
    rphi = roots(phi);
    if (any(real(rphi) > 0)) || (~isreal(rphi))
        j = (real(rphi) > 0) & imag(rphi) == 0 | (imag(rphi) ~= 0);
        prphi = rphi(j);
        lprphi = length(prphi);
        phir = 1.;
        for i = 1:lprphi
            phir1 = [-prphi(i), 1];
            phir = conv(phir, phir1);
        end
        phi = fliplr(phi);
        phis = deconv(phi, phir);
    else
        phis = fliplr(phi);
        phir = 1.;
    end
else
    phir = 1.;
    phis = 1.;
end

if ~isempty(Phi)
    Phi = squeeze(Phi)';
    [nP, mP, lP] = size(Phi); %mP - 1 = deg(Phi) <= 1
    if mP > 1
        Phix = zeros(1, s+1);
        Phix(1) = 1.;
        Phix(end) = Phi(end);
        % Find positive roots and assign them to phir
        % The rest of the roots are assigned to phis
        rPhi = roots(Phix);
        if any(real(rPhi) > 0)
            j = (real(rPhi) > 0) & (imag(rPhi) == 0);
            nrPhi = rPhi(j);
            lnrPhi = length(nrPhi);
            phirx = 1.;
            for i = 1:lnrPhi
                phir1 = [-nrPhi(i), 1];
                phirx = conv(phirx, phir1);
            end
            phir = conv(phir, phirx);
            Phix = fliplr(Phix);
            phisx = deconv(Phix, phirx);
            phis = conv(phis, phisx);
        else
            phis = conv(phis, fliplr(Phix));
        end
    end
end

if ~isempty(th)
    th = squeeze(th)';
    th = fliplr(th);
else
    th = 1.;
end

if ~isempty(Th)
    Th = squeeze(Th)';
    Th = fliplr(Th);
else
    Th = 1;
end

% Regular differencing

if dr >= 0 && ds >= 0 && (ceil(dr) == floor(dr)) && (ceil(ds) == floor(ds))
    rd = dr + ds;
    if rd > 0
        phird1 = [-1, 1];
        phird = phird1;
        for i = 1:rd - 1
            phird = conv(phird, phird1);
        end
    else
        phird = 1;
    end
end

phirst = phir;
phir = conv(phir, phird);


% Seasonal differencing

if ds > 0 && (ceil(ds) == floor(ds))
    phisd = ones(1, s);
    phisd1 = phisd;
    for i = 1:ds - 1
        phisd = conv(phisd, phisd1);
    end
elseif ds == 0
    phisd = 1;
end
phis = conv(phis, phisd);


% Regular polynomial thr
thr = th;

% Seasonal polynomial ths
if isempty(Th)
    ths = 1;
else
    sma = length(Th) - 1;
    ths = zeros(1, 1+sma*s);
    for i = 0:sma
        ths(1+i*s) = Th(i+1);
    end
end

end