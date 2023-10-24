function plotStressCompAlongCrack()
    period = 2;
    hl = 0.5;
    x0 = 0;
    y0 = 0;
    phi0 = 0;
   %st_func = @st_single_crack;
    st_func = @() st_collinear_crack_row(period);
    [st_func_mode1, st_func_mode2] = st_func();
    
    xList = -hl : 0.0001*hl : hl;
    numPoints = numel(xList);
    stMode1List = zeros(numPoints, 3);
    stMode2List = zeros(numPoints, 3);

    for i = 1:numPoints
        x = xList(i);
        y = 0.001;

        st_mode1 = st_func_mode1(x, y, x0, y0, phi0, hl);
        st_mode2 = st_func_mode2(x, y, x0, y0, phi0, hl);

        stMode1List(i, :) = [ st_mode1(1, 1), st_mode1(2, 2), st_mode1(1, 2) ];
        stMode2List(i, :) = [ st_mode2(1, 1), st_mode2(2, 2), st_mode2(1, 2) ];
    end

    figure("Name", "Mode I");
    plot(xList, stMode1List(:, 1));    

    figure("Name", "Mode II");
    plot(xList, stMode2List(:, 2));
end