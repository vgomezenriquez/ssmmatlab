function X = pmmulbf(A, B);
%  multiplies matrix pol. in z (A0+A1*z+...)by matrix pol.in z^{-1}
% (B0+B1*z^{-1}+...)'
% in output:X(1) is the coeff. of z^{-1}(max)
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

[na, ma, p] = size(A);
[nb, mb, q] = size(B);

%the following works only for A and B having the same size
% C=B(:,:,1);
% for i=1:floor(q/2)
%  B(:,:,i)=B(:,:,q-i+1)'; B(:,:,q-i+1)=C'; C=B(:,:,i+1);
% end
% if mod(q,2) ~= 0
%  i=floor(q/2)+1;
%  B(:,:,i)=B(:,:,i)';
% end

%this works for the general case
C = B(:, :, 1);
D = zeros(mb, nb, q);
for i = 1:floor(q/2)
    D(:, :, i) = B(:, :, q-i+1)';
    D(:, :, q-i+1) = C';
    C = B(:, :, i+1);
end
if mod(q, 2) ~= 0
    i = floor(q/2) + 1;
    D(:, :, i) = B(:, :, i)';
end


X = pmatmul(A, D);
