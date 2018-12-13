function ser = W12
%Wei's example, 'Time Series Analysis, 1990, pp. 302-308'
M = load(fullfile('data', 'W12.dat')); % 2 series
% n=size(M,1);% Longitud de las series
% M(5:7,1)=NaN(3,1);   %we add some missing values
% M(:,1)=M(:,1) + (100.).*(1:size(M,1))'; we add a mean for the dif. series
ser.bg_year = 1000;
ser.bg_per = 1;
ser.freq = 1; % Fecha ficticia 4 cifras
ser.ninput = 1; % Hay 1 serie explicativa
ser.Yin = M(:, 2); % Serie explicativa
ser.yor = M(:, 1); % Serie explicada
ser.nlagtf = -1; % Nro. retardos para ident. automatica determinado por el programa
ser.tfident = 1; % Se hace identificacion automatica
ser.pr = 1; %Se crea informe-resumen
ser.gft = 1; %Se hacen graficos
ser.dr = 1;
ser.fixdif = 0;
ser.prelivar = 0;
%input fictitious forecasts
npr = 12;
ser.npr = npr;
twxn = 2 * M(end, 2);
ser.modpred(1).pred = zeros(npr, 1);
ser.modinput.mod = 0;
for i = 1:npr
    ser.modpred(1).pred(i) = twxn - M(end-i, 2); %use symmetry to obtain the forecasts of x1
end
