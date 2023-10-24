function avgTr = evaluateAverageTractions()
    function mean = functionMean(func, from, to)
        mean = integral(func, from, to, 'ArrayValued', true, 'RelTol',1e-6,'AbsTol',1e-10) / (to - from);
    end

    period = 2;
    hl = 0.5;
    [stFuncMode1, stFuncMode2] = st_collinear_crack_row(period);

    c1Center = [ 0; 0.5 ];
    c2Center = [ 0;-0.5 ];
    norm = [ 0; 1 ];
    tang = [ 1; 0 ];
    extStress = [ 1; 0 ];

    B = [ extStress; extStress ];
    A = eye(4, 4);

    sigC1OnC2Mode1 = functionMean(@(s) stFuncMode1(s, c2Center(2), c1Center(1), c1Center(2), 0, hl), -hl, hl);
    sigC1OnC2Mode2 = functionMean(@(s) stFuncMode2(s, c2Center(2), c1Center(1), c1Center(2), 0, hl), -hl, hl);
    sigC2OnC1Mode1 = functionMean(@(s) stFuncMode1(s, c1Center(2), c2Center(1), c2Center(2), 0, hl), -hl, hl);
    sigC2OnC1Mode2 = functionMean(@(s) stFuncMode2(s, c1Center(2), c2Center(1), c2Center(2), 0, hl), -hl, hl);

    A(1, 3) = - norm' * sigC2OnC1Mode1 * norm;
    A(1, 4) = - norm' * sigC2OnC1Mode2 * norm;
    A(2, 3) = - norm' * sigC2OnC1Mode1 * tang;
    A(2, 4) = - norm' * sigC2OnC1Mode2 * tang;
    A(3, 1) = - norm' * sigC1OnC2Mode1 * norm;
    A(3, 2) = - norm' * sigC1OnC2Mode2 * norm;
    A(4, 1) = - norm' * sigC1OnC2Mode1 * tang;
    A(4, 2) = - norm' * sigC1OnC2Mode2 * tang;

    avgTr = A \ B;
end