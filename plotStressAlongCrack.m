function plotStressAlongCrack()

    period = 2;
    hl = 0.5;

    function res = syy(x)
        S = sqrt(sin(pi*x/period).^2 - sin(pi*hl/period)^2);

        res = imag(1i*(sin(pi*x/period)./S - 1));
    end

    st_func = @() st_collinear_crack_row(period);
%     st_func = @() st_single_crack();
    [st_func_mode1, st_func_mode2] = st_func();
    
    %stFuncAlongCrack = @(x) sqrt((hl+x)/(hl-x)) * [0, 1] * syy(x) * [0; 1];
    stFuncAlongCrack = @(x) cos(x)/sqrt(pi)/sin(hl) * sqrt(tan(hl)) * sqrt((sin(hl)+sin(x))/(sin(hl)-sin(x)));


    stMean = integral(stFuncAlongCrack, -hl, hl, 'ArrayValued', true,'AbsTol',1e-17);
        
    xvals = -hl : 0.00001*hl : hl;
    yvals = syy(xvals);

    plot(xvals, yvals);
end