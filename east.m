function Y = east(Iy, Im, N, Idur, Mq, Yd)
%
% this function generates the variable used to correct for Easter
% effect. Given the initial year Iy, the initial month Im, the desired
% length of the series N, and the duration of the effect Idur, the function
% computes for each month the proportion of the period Idur before Easter
% which falls in that month. It works for the period 1901-2099.
%
% input variables       Iy      : the initial year
%                       Im      : the initial period
%                       N       : the length of the desired vector
%                       Idur    : the length of the period before Easter
%                                 that the effect is thought to prevail
%                       Mq      : the series frequency (=12 for monthly,
%                                                       =4 for quarterly)
%                       Yd      : 199 x 2 array containing the dates of
%                                 Easter
% output variables      Y       : N x 1 array containing the variable%
% Copyright (c) 21 July 2015 by Victor Gomez
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


if (Mq == 4)
    iim = Im;
    in = N;
    Im = (Im - 1) * 3 + 1;
    N = N * 3;
end
Y = zeros(N, 1);
ian = Iy - 1900;
icount = 13 - Im;
ic = 13 - Im;
ik = 0;
istop = 0;
dur = double(Idur);
while 1
    if (icount >= N), istop = 1;
        Im = 12 + N - icount;
    end
    imes = Yd(ian, 1);
    idia = Yd(ian, 2);
    if (istop == 1)
        if (Im >= imes - 1), Im = 1;
        else, Im = 12;
        end
    end
    ind = imes - Im + 1;
    %    if (imes == 3), text='MARCH'; else, text='APRIL'; end
    if (ik > 0), ind = ind + ic + (ik - 1) * 12;
    end
    if (ind <= N)
        idia1 = idia - 1;
        if (imes > Im)
            if (Idur <= idia1)
                Y(ind) = .5D0;
                if (imes == 4), Y(ind-1) = -.5D0;
                else, ...
                        if (ind + 1 <= N), Y(ind+1) = -.5D0;
                        end, end
            else
                Y(ind) = (Yd(ian, 2) - 1.D0) / dur - .5D0;
                Y(ind-1) = .5D0 - (Yd(ian, 2) - 1.D0) / dur;
            end
        elseif (imes == Im)
            if (Idur <= idia1)
                Y(ind) = .5D0;
                if (imes == 3) & (ind + 1 <= N), Y(ind+1) = -.5D0;
                end
            else
                Y(ind) = (Yd(ian, 2) - 1.D0) / dur - .5D0;
            end
        end
        ian = ian + 1;
        icount = icount + 12;
        Im = 1;
        ik = ik + 1;
    else
        istop = 1;
    end
    if (istop == 0), continue, end
    if (Mq == 4)
        Im = iim;
        N = in;
        for i = 1:N
            ii = (i - 1) * 3;
            Y(i) = sum(Y(ii+1:ii+3));
        end
        Y(N+1:3*in) = [];
    end
    break
end
