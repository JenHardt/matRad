function resultGUI = matRad_TopasRecalcRBEfromLET(resultGUI,modelName,alphaX,betaX)

if nargin <=2
    alphaX = 0.1;
    betaX = 0.05;
end

if isfield(resultGUI,'LET')

else
    numOfBeams = sum(contains(fieldnames(resultGUI),'LET_beam'));
end

if isfield(resultGUI,'w')
    wBeam = resultGUI.w;
else
    wBeam = ones(resultGUI.totalNumOfBixels);
end
if isfield(resultGUI,'LET')
    numOfVoxel = numel(resultGUI.LET);
    numOfBeams = sum(contains(fieldnames(resultGUI),'LET_beam'));
    doseGrid = size(resultGUI.LET);
    dij = 0;
elseif isfield(resultGUI,'mLETDose')
    numOfVoxel = size(resultGUI.mLETDose{1},1);
    numOfBeams = resultGUI.numOfBeams;
    numOfBeamlets = size(resultGUI.mLETDose{1},2);
    dij = 1;
end

ABratio = alphaX / betaX;

modelParam.p0_MCN   = 0.999064;     % according to https://www.ncbi.nlm.nih.gov/pubmed/26459756
modelParam.p1_MCN   = 0.35605;
modelParam.p2_MCN   = 1.1012;
modelParam.p3_MCN   = -0.0038703;

modelParam.p0_WED   = 1; % https://www.ncbi.nlm.nih.gov/pubmed/22909391
modelParam.p1_WED   = 0.434;
modelParam.p2_WED   = 1;

if ~dij
    mAlphaDose = spalloc(numOfVoxel,numOfBeams,1);
    mSqrtBetaDose = spalloc(numOfVoxel,numOfBeams,1);

    for i = 1:numOfBeams
        beamInfo(i).suffix = ['_beam', num2str(i)];
        beamInfo(i).logIx  = ([1:numOfBeams]' == i);
    end
    beamInfo(numOfBeams+1).suffix = '';
    beamInfo(numOfBeams+1).logIx  = true(numOfBeams);

    for i = 1:length(beamInfo)-1
        LET = resultGUI.(['LET', beamInfo(i).suffix]);
        LET(isnan(LET)) = 0;

        switch modelName
            case 'MCN'

                RBEmax     = modelParam.p0_MCN + ((modelParam.p1_MCN * LET )./ ABratio);
                RBEmin     = modelParam.p2_MCN + (modelParam.p3_MCN  * sqrt(ABratio) .* LET);

                resultGUI.(['alpha_recalc_MCN', beamInfo(i).suffix]) = RBEmax    .* alphaX;
                resultGUI.(['beta_recalc_MCN', beamInfo(i).suffix])  = RBEmin.^2 .* betaX;

            case 'WED'

                RBEmax     = modelParam.p0_WED + ((modelParam.p1_WED * LET )./ ABratio);
                RBEmin     = modelParam.p2_WED;

                resultGUI.(['alpha_recalc_WED', beamInfo(i).suffix]) = RBEmax    .* alphaX;
                resultGUI.(['beta_recalc_WED', beamInfo(i).suffix])  = RBEmin.^2 .* betaX;

        end

        mAlphaDose(:,i)         = reshape(resultGUI.(['alpha_recalc_' modelName, beamInfo(i).suffix]) .* resultGUI.(['physicalDose', beamInfo(i).suffix]),[],1);
        mSqrtBetaDose(:,i)      = reshape(sqrt(resultGUI.(['beta_recalc_' modelName, beamInfo(i).suffix])) .* resultGUI.(['physicalDose', beamInfo(i).suffix]),[],1);

    end

    for i = 1:length(beamInfo)
        resultGUI.(['effect_recalc_' modelName, beamInfo(i).suffix])       = full(mAlphaDose * beamInfo(i).logIx + (mSqrtBetaDose * beamInfo(i).logIx).^2);
        resultGUI.(['effect_recalc_' modelName, beamInfo(i).suffix])       = reshape(resultGUI.(['effect_recalc_' modelName, beamInfo(i).suffix]),doseGrid);

        resultGUI.(['RBExD_recalc_' modelName, beamInfo(i).suffix])        = zeros(size(resultGUI.(['effect_recalc_' modelName, beamInfo(i).suffix])));
        resultGUI.(['RBExD_recalc_' modelName, beamInfo(i).suffix])        = (sqrt(alphaX.^2 + 4 .* betaX .* resultGUI.(['effect_recalc_' modelName, beamInfo(i).suffix])) - alphaX)./(2.*betaX);

        resultGUI.(['RBE_recalc_' modelName, beamInfo(i).suffix])          = resultGUI.(['RBExD_recalc_' modelName, beamInfo(i).suffix])./resultGUI.(['physicalDose', beamInfo(i).suffix]);
    end

else
%     resultGUI.mLETDose{1}(isnan(resultGUI.mLETDose{1})) = 0;

    for beamLet = 1:numOfBeamlets
        switch modelName
            case 'MCN'

                RBEmax     = modelParam.p0_MCN + ((modelParam.p1_MCN * full(resultGUI.mLETDose{1}(:,beamLet)) )./ ABratio);
                RBEmin     = modelParam.p2_MCN + (modelParam.p3_MCN  * sqrt(ABratio) .* full(resultGUI.mLETDose{1}(:,beamLet)));

                resultGUI.alpha_recalc{1}(:,beamLet) = RBEmax    .* alphaX;
                resultGUI.beta_recalc{1}(:,beamLet)  = RBEmin.^2 .* betaX;

            case 'WED'

                RBEmax     = modelParam.p0_WED + ((modelParam.p1_WED * full(resultGUI.mLETDose{1}(:,beamLet)) )./ ABratio);
                RBEmin     = modelParam.p2_WED;

                resultGUI.alpha_recalc{1}(:,beamLet) = RBEmax    .* alphaX;
                resultGUI.beta_recalc{1}(:,beamLet)  = RBEmin.^2 .* betaX;

        end
    end

    resultGUI.mAlphaDose_recalc{1}      = resultGUI.alpha_recalc{1} .* resultGUI.physicalDose{1};
    resultGUI.mSqrtBetaDose_recalc{1}   = sqrt(resultGUI.beta_recalc{1}) .* resultGUI.physicalDose{1};

    resultGUI = rmfield(resultGUI,{'alpha_recalc','beta_recalc'});
end

resultGUI = orderfields(resultGUI);
end