function plotSifByDefinition()  % i.e. via limits
    period = 2;
    hl = 0.5;
    x0 = 0;
    y0 = 0;
    phi0 = 0;
    st_func = @st_single_crack;
   %st_func = @() st_collinear_crack_row(period);
    [st_func_mode1, st_func_mode2] = st_func();
    
    initDist = (period - 2*hl);
    distList = initDist/10 : -initDist/100 : initDist/100;
    numPoints = numel(distList);
    sif1 = zeros(numPoints, 1);
    sif2 = zeros(numPoints, 1);

    for i = 1:numPoints
        x = hl + distList(i);
        y = 0;

        st_mode1 = st_func_mode1(x, y, x0, y0, phi0, hl);
        st_mode2 = st_func_mode2(x, y, x0, y0, phi0, hl);

        sif1(i) = sqrt(2*pi*distList(i)) * st_mode1(2, 2);
        sif2(i) = sqrt(2*pi*distList(i)) * st_mode1(1, 2);
    end

    figure("Name", "SIF I");
    plot(distList, sif1);
    

    figure("Name", "SIF II");
    plot(distList, sif2);    
end