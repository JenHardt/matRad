function [resultGUI,pln] = matRad_calcDoseModulated(ct,stf,pln,cst,weights)

% Instance of matRad configuration class
matRad_cfg =  MatRad_Config.instance();

% load default parameters in case they haven't been set yet
pln = matRad_cfg.getDefaultProperties(pln,{'propDoseCalc','propHeterogeneity','propMC'});

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
    case 'MCsquare'
        calcExternal = false;
        histories = pln.propMC.histories;
        if ~isfield(pln.propMC,'materialConverter') || ~isfield(pln.propMC.materialConverter,'addSection')
            pln.propMC.materialConverter.addSection = 'sampledDensities';
        end
    case 'matRad'
        calcExternal = false;
    otherwise
        matRad_cfg.dispError('No sampling mode other than TOPAS and matRad implemented');
end

% parallelComputationTOPAS =1;
% if parallelComputationTOPAS
%
%     for i = 1:samples
%         ct_mod{i} = heterogeneityConfig.modulateDensity(ct,cst,pln,Pmod,modulation);
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

% Instance of HeterogeneityCorrection configuration class
heterogeneityConfig = MatRad_HeterogeneityConfig.instance();

% Perform resampling to dose grid if necessary (modulation is performed on the resampled grid)
switch pln.propHeterogeneity.sampling.mode
    case 'TOPAS'
        topasConfig = MatRad_TopasConfig();
        pln.propMC.proton_engine = 'TOPAS';
        [ctR,cstR,stfR] = topasConfig.resampleGrid(ct,cst,pln,stf);
    case 'MCsquare'
        pln.propMC.proton_engine = 'MCsquare';
        ctR = ct;
        cstR = cst;
        stfR = stf;
    case 'matRad'
        ctR = ct;
        cstR = cst;
end

% set this flag so that the modulated cube is not overwritten in matRad_calcDoseInit
pln.propDoseCalc.useGivenEqDensityCube = true;

% Allocate empty resultGUI and space for individual physical doses to calculate their standard deviation
resultGUI = struct;
data = cell(samples,1);

for i = 1:samples
    % Write current sample to the console
    matRad_cfg.dispInfo([pln.propHeterogeneity.sampling.mode,': Dose calculation for CT ' num2str(i) '/' num2str(samples)])

    % Modulate density of ct cube
    ct_mod = heterogeneityConfig.modulateDensity(ctR,cstR,pln);

    % Save number of samples in modulated CT (e.g. used for TOPAS folder generation)
    ct_mod.sampleIdx = i;

    % Switch between the different modes of calculation currently implemented
    % WARNING: Implementation of MCsquare is currently not finished
    switch pln.propHeterogeneity.sampling.mode
        case {'TOPAS','MCsquare'}
            % Set TOPAS parameters
            pln.propMC.numOfRuns = 1;
            pln.propHeterogeneity.sampling.histories = histories/samples;

            % Calculate dose with modulated CT
            resultGUI_mod = matRad_calcDoseDirectMC(ct_mod,stfR,pln,cstR,weights);

            if ~calcExternal
                % Accumulate averaged results
                resultGUI = heterogeneityConfig.accumulateOverSamples(resultGUI,resultGUI_mod,samples);

                % Save individual standard deviation
                if isfield(resultGUI_mod,'physicalDose_std')
                    resultGUI.physicalDose_std_individual{i} = resultGUI_mod.physicalDose_std;
                end
            end
        case 'matRad'
            % Calculate dose with modulated CT
            resultGUI_mod = matRad_calcDoseDirect(ct_mod,stf,pln,cstR,weights);

            % Accumulate averaged results
            resultGUI = heterogeneityConfig.accumulateOverSamples(resultGUI,resultGUI_mod,samples);
    end

    % Save individual physical doses to calculate standard deviation
    if ~calcExternal
        data{i} = resultGUI_mod.physicalDose;
    end

end

if ~calcExternal
    % Calculate standard deviation between samples
    resultGUI.physicalDose_std = heterogeneityConfig.calcSampleStd(data,resultGUI.physicalDose);
end

% Change loglevel back to default;
matRad_cfg.logLevel = logLevel;

end

