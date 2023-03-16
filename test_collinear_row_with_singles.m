periods = [ 2 ];
hls = 0.5;
numCracksVert = 2;
vertSpacing = 1;
bstLst = zeros(2, 2, 3);
    bstLst(2, 2, 1) = 1;
    bstLst(1, 2, 2) = 1;
    bstLst(2, 1, 2) = 1;
    bstLst(1, 1, 3) = 1;
    bstLst(2, 2, 3) = 1;

resultsDir = sprintf("results/collinear_vs_singles/hl=%.3f_numVertCracks=%i_vertSpacing=%0.3f/", ...
    hls, ...
    numCracksVert, ...
    vertSpacing);
if ~exist(resultsDir, "dir")
    mkdir(resultsDir);
else
    %error("Error! Directory %s already exists!", resultsDir);
end

for iPeriod = 1:numel(periods)
    period = periods(iPeriod);
    fnDetails = sprintf("period=%0.3f", period);

    [sifColSol, avgTrColSol] = evaluateUsingCollinearSt(periods(1), hls(1), numCracksVert, vertSpacing, bstLst);
    save(strcat(resultsDir, "collinear-row-st-sol_", fnDetails, ".mat"), "sifColSol", "avgTrColSol");
end





function [sifMat, avgTract] = evaluateUsingCollinearSt(period, hl, numCracksVert, vertSpacing, bstLst)

    ccLst = zeros(numCracksVert, 2);
    cpLst = zeros(numCracksVert, 1);
    chlLst = hl*ones(numCracksVert, 1);
    sifCrackEvalLst = 1:numCracksVert;
    st_func = @() st_collinear_crack_row(period);

    for i = 2:numCracksVert
        ccLst(i, 2) = ccLst(i - 1, 2) - (-1)^i * (i - 1) * vertSpacing;
    end

    if false
        figure;    
        plot(ccLst(:, 1), ccLst(:, 2), 'o');
    end

    [sifMat, avgTract] = crack_interact(ccLst, cpLst, chlLst, bstLst, sifCrackEvalLst, st_func);
end

function [sifMat, avgTract] = evaluateUsingSingleCrackSt(period, hl, numCracksVert, vertSpacing, bstLst)

    numHorCracksOneSide = 5; % how many horizontal cracks are added for each crack in the stack from ONE SIDE (i.e. the total number is double)
    numCracks = (2*numHorCracksOneSide + 1)*numCracksVert;

    ccLst = zeros(numCracks, 2);
    cpLst = zeros(numCracks, 1);
    chlLst = hl*ones(numCracks, 1);
    sifCrackEvalLst = 1:numCracksVert;
    st_func = @st_single_crack;

    % CONTINUE HERE! Stopped at 13:00, 16/03/23

    for i = 2:numCracksVert
        ccLst(i, 2) = ccLst(i - 1, 2) - (-1)^i * (i - 1) * vertSpacing;
    end

    if false
        figure;    
        plot(ccLst(:, 1), ccLst(:, 2), 'o');
    end

    [sifMat, avgTract] = crack_interact(ccLst, cpLst, chlLst, bstLst, sifCrackEvalLst, st_func);
end