function [resultGUI,pln] = matRad_calcDoseModulated(ct,stf,pln,cst,weights)

global matRad_cfg;
matRad_cfg =  MatRad_Config.instance();

pln = matRad_cfg.loadDefaultParam(pln);


samples = pln.propHeterogeneity.sampling.numOfSamples;

switch pln.propHeterogeneity.sampling.mode
    case 'TOPAS'
        %         pln.propMC.materialConverter.densityCorrection = 'TOPAS2'; %'default','TOPAS1','TOPAS2'
        pln.propMC.materialConverter.HUSection = 'advanced'; %'default','advanced'
        pln.propMC.materialConverter.HUToMaterial = 'default'; %'default','simpleLung','advanced'
        if ~isfield(pln.propMC.materialConverter,'addSection')
            if strcmp(pln.propHeterogeneity.sampling.method,'poisson')
                pln.propMC.materialConverter.addSection = 'poisson';
            else
                pln.propMC.materialConverter.addSection = 'sampledDensities'; %'none','lung','poisson','sampledDensities' (the last 2 only with modulation)
            end
        end

        histories = pln.propMC.histories;
        calcExternal = pln.propMC.externalCalculation;
    case 'matRad'

    otherwise
        matRad_cfg.dispError('No sampling mode other than TOPAS and matRad implemented');
end

% parallelComputationTOPAS =1;
% if parallelComputationTOPAS
%
%     for i = 1:samples
%         ct_mod{i} = matRad_modulateDensity(ct,cst,pln,Pmod,modulation);
%     end
%
%     pln.propMC.proton_engine = 'TOPAS';
%         if strcmp(modulation,'poisson')
%             pln.propMC.materialConverter = 'HUToWaterSchneider_mod';
%         else
%             if ~isfield(pln.propMC,'materialConverter')
%                 pln.propMC.materialConverter = 'HUToWaterSchneider';
%             end
%         end
%
%         resultGUI_mod = matRad_calcDoseDirectMC(ct_mod,stf,pln,cst,weights,mode{2}/samples,mode{3});
%
%         if ~mode{3}
%             %     resultGUI.(['physicalDose',num2str(s)]) = resultGUI.(['physicalDose',num2str(s)]) + resultGUI_mod.physicalDose/s;
%             if strcmp(pln.bioParam.quantityOpt,'RBExD')
%                 resultGUI.RBExD = resultGUI.RBExD + resultGUI_mod.RBExD/samples;
%             end
%             resultGUI.physicalDose = resultGUI.physicalDose + resultGUI_mod.physicalDose/samples;
%             std{i} = resultGUI_mod.std;
%         end
%
%
%
% else

% Turn info and warning messages off for modulation
logLevel = matRad_cfg.logLevel;
matRad_cfg.logLevel = 1;

switch pln.propHeterogeneity.sampling.mode
    case 'TOPAS'
        [ctR,cstR] = matRad_resampleTopasGrid(ct,cst,pln,stf);
    case 'matRad'
        ctR = ct;
        cstR = cst;
end
resultGUI.physicalDose = zeros(ct.cubeDim);
if strcmp(pln.bioParam.quantityOpt,'RBExD')
    resultGUI.RBExD = zeros(ct.cubeDim);
end

% set this flag so that the modulated cube is not overwritten in matRad_calcDoseInit
pln.propDoseCalc.useGivenEqDensityCube = true;
data = cell(samples,1);
std = cell(samples,1);

for i = 1:samples
    fprintf([pln.propHeterogeneity.sampling.mode,': Dose calculation for CT %i/%i \n'],i,samples)
    ct_mod = matRad_modulateDensity(ctR,cstR,pln);
    ct_mod.sampleIdx = i;
    %%
    switch pln.propHeterogeneity.sampling.mode
        case 'TOPAS'
            pln.propMC.proton_engine = 'TOPAS';
            pln.propMC.numOfRuns = 1;
            pln.propHeterogeneity.sampling.histories = histories/samples;
            resultGUI_mod = matRad_calcDoseDirectMC(ct_mod,stf,pln,cstR,weights);

            if ~calcExternal
                %     resultGUI.(['physicalDose',num2str(s)]) = resultGUI.(['physicalDose',num2str(s)]) + resultGUI_mod.physicalDose/s;
                if strcmp(pln.bioParam.quantityOpt,'RBExD')
                    resultGUI.RBExD = resultGUI.RBExD + resultGUI_mod.RBExD/samples;
                end
                resultGUI.physicalDose = resultGUI.physicalDose + resultGUI_mod.physicalDose/samples;
                std{i} = resultGUI_mod.physicalDose_std;
            end
        case 'matRad'
            resultGUI_mod = matRad_calcDoseDirect(ct_mod,stf,pln,cstR,weights);

            if strcmp(pln.bioParam.quantityOpt,'RBExD')
                resultGUI.RBExD = resultGUI.RBExD + resultGUI_mod.RBExD/samples;
            end
            resultGUI.physicalDose = resultGUI.physicalDose + resultGUI_mod.physicalDose/samples;
    end

    data{i} = resultGUI_mod.physicalDose;

end

% Calculate Standard deviation between samples
meanDiff = 0;
for k = 1:samples
    meanDiff = meanDiff + (data{k} - resultGUI.physicalDose).^2;
end
varMean = meanDiff./(samples - 1)./samples;
stdMean = sqrt(varMean);

stdSum = stdMean * samples;
varSum = stdSum.^2;

resultGUI.physicalDose_std = sqrt(varSum);

if strcmp(pln.propHeterogeneity.sampling.mode,'TOPAS') && ~calcExternal
    resultGUI.physicalDose_std_individual = std;
end

%Change loglevel back to default;
matRad_cfg.logLevel = logLevel;

end

