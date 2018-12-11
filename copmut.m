function [U, T] = copmut(Um, Tm, ColsUT, da, du, n, m, full)
% this function computes the Unimodular polynomial matrices U and T
% Author Felix Aparicio-Perez, Instituto Nacional de Estadistica, Spain
ngradT = da + du;
T = zeros(n, m, ngradT+1);
U = zeros(m, m, du+1);
if (full == 0)
    for k = 1:du + 1
        U(:, :, du+2-k) = Um(m*(k - 1)+1:m*k, ColsUT);
    end
end
Tm = permat(Tm, da+du+1, n);
%Tm = P'*Tm
for k = 1:ngradT + 1
    T(:, :, ngradT+2-k) = Tm(n*(k - 1)+1:n*k, ColsUT);
end
if (full)
    U = T(n-m+1:end, :, 1:du+1); % We extract the U matrix from the lower part of the T matrix
    % (see (2) in Henrion and Sebek) but U has degree du only
    T = T(1:n-m, :, :); % We eliminate from T the U matrix that we have extracted
end
