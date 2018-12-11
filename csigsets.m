function [Sig, ColsUT, Cs, Ncs] = csigsets(SigBar, dtot, du, n, m)
% This routine extracts the Ci sets and, if possible, a triangular shape of
% size m
% Author Felix Aparicio-Perez, Instituto Nacional de Estadistica, Spain

Ncs = zeros(n, 1);
ColsUT = zeros(m, 1);
Cs = zeros(n, m*(du + 1));
for i = 1:n
    icual = find((SigBar <= i * (dtot + 1)) & (SigBar >= i * (dtot + 1) - dtot));
    Ncs(i) = size(icual, 2);
    Cs(i, 1:Ncs(i)) = icual;
end
icual2 = find(Ncs ~= 0);
if (max(size(icual2)) >= m) %En este caso hay solución
    for j = 1:m
        Sig(j) = SigBar(Cs(icual2(j), Ncs(icual2(j))));
        ColsUT(j) = Cs(icual2(j), Ncs(icual2(j)));
    end
else
    Sig = [];
end
