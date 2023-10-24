function kachCalculate()
                
    resFn = "";
    
    bstLst = zeros(2, 2, 1);
        bstLst(2, 2, 1) = 1;
    st_func = @() st_single_crack();

    bigCrLen = 1; % m

    smCrLenLst = [ bigCrLen/10, bigCrLen/2, bigCrLen ];

    for iSmCrLen = 1:length(smCrLenLst)
        smCrLen = smCrLenLst(iSmCrLen);

        smCrCoordXLst = 0 : 1/4*((bigCrLen+smCrLen)/2) : (bigCrLen+smCrLen)/2;
        vertSpaceLst = [ smCrLen/4, smCrLen/2, smCrLen/1, bigCrLen/3, 2*bigCrLen/3 ];

        for iSmCrCoordX = 1:length(smCrCoordXLst)
            smCrCoordX = smCrCoordXLst(iSmCrCoordX);

            for iVertSpaceLst = 1:length(vertSpaceLst)
                vertSpace = vertSpaceLst(iVertSpaceLst);
                
                smCrCoordY = vertSpace;
                
                ccLst = [ 0, 0;
                          smCrCoordX, smCrCoordY ];
                hlLst = [ bigCrLen, smCrLen ];
            
                [sifM, avgM] = crack_interact(ccLst, zeros(length(hlLst), 1), hlLst, bstLst, 1:2, st_func);

            end
        end
    end



end