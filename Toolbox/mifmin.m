function [y, fval, exitflag] = mifmin(r, den, x1, x2)
%
% minimization function called in candec

options = optimset('TolX', 1e-12, 'Display', 'off'); % Turn off Display
% options = optimset('TolX',1e-12,'Display', 'iter');
if isOctave
    polyW = @(x) poly(x, r, den);
    [y, fval, exitflag] = fminbnd(polyW, x1, x2, options);
else
    [y, fval, exitflag] = fminbnd(@poly, x1, x2, options);
end
    function y = poly(x) % Compute the rational expression.
        y = polyval(r, x) / polyval(den, x);
    end
end
