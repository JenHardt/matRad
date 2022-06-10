function ct = matRad_modulateDensity(ct,cst,pln)
% matRad density modulation function to calculate sampled ct struct
%
% call
%   ct = matRad_modulateDensity(ct,cst,Pmod,mode)
%
% input
%   ct:             ct struct
%   cst:            matRad cst struct
%   Pmod:           Modulation power according to which the modulation will
%                   be created
%   mode:           mode for density modulation ('binominal','poisson')
%                   note: poisson only available for Pmod = 250,450 and 800
%
% output
%   ct:             ct struct with modulated density cube
%
% References
%   [1] Poisson sampling: https://iopscience.iop.org/article/10.1088/1361-6560/aa641f
%   [2] Detwiler et al., 2021, Compendium of Material Composition Data for Radiation Transport Modeling
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Copyright 2020 the matRad development team.
%
% This file is part of the matRad project. It is subject to the license
% terms in the LICENSE file found in the top-level directory of this
% distribution and at https://github.com/e0404/matRad/LICENSES.txt. No part
% of the matRad project, including this file, may be copied, modified,
% propagated, or distributed except according to the terms contained in the
% LICENSE file.
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Instance of MatRad_Config class
matRad_cfg = MatRad_Config.instance();

% get all unique lung indices from lung segmentations
idx = contains(cst(:,2),'Lung');
if sum(idx)==0
    matRad_cfg.dispError('No lung segmentation found in cst.\n');
end
lungIdx = [cst{idx,4}];
lungIdx = unique(vertcat(lungIdx{:}));

% calculate ct cube from cubeHU if not specified
if ~isfield(ct,'cube')
    ct = matRad_calcWaterEqD(ct,pln);
end

switch pln.propHeterogeneity.sampling.method
    case 'binomial'
        % set lung tissue density
        rhoLung = 1.05; % [2] [g/cm^3]

        % isolate lung densities and normalize to lung tissue density
        pLung = ct.cube{1}(lungIdx) / rhoLung;

        % scrap lung densities larger than 1 (cannot be used in sampling, since distributions are defined over [0,1])
        if any(pLung > 1)
            lungIdx = lungIdx(pLung <= 1);
            pLung = ct.cube{1}(lungIdx) / rhoLung;
        end

        % calculate size of substructures
        d = pln.propHeterogeneity.modPower/1000 ./ (1-pLung) / rhoLung; % [1] eq.8: Pmod = d*(1-pLung) * rhoLung

        % length of a voxel (uniform voxels are assumed)
        D = ct.resolution.y;

        % calculate number of substructures inside of a single voxel (round in case of discrete distribution)
        if pln.propHeterogeneity.sampling.continuous
            n = D./d;
            largeEnough = n > 1;
        else
            n = round(D./d);
            largeEnough = n >= 1;
        end

        % Don't modulate voxel with less than 1 substructures
        lungIdx = lungIdx(largeEnough);
        pLung = pLung(largeEnough);
        n = n(largeEnough);
        
        % get samples from the binomial distribution (discrete or continuous approximation)
        samples = matRad_sampleBino(n,pLung,length(lungIdx),pln.propHeterogeneity.sampling.continuous);
        
        % revert normalization to get values between [0,rhoLung]
        samples = samples * rhoLung;
        
        % write samples to CT and convert to Hounsfield Units
        ct.cube{1}(lungIdx) = samples;
        ct.cubeHU{1}(lungIdx) = 1024*(ct.cube{1}(lungIdx)-1);

    case 'poisson'
        % read density modulation look-up table for supported modulation powers (Pmod = 250, 450, 800 mu)
        DensMod = matRad_loadModDist(pln.propHeterogeneity.modPower);

        % Connect each density-probability-couple with a number that will later be transformed the HU-Value:
        % The maximum HU-Value of the HU set by Schneider etal is 2995: So the HU of the modulated density must be at least 2995+1=2996 ; This values is prerocessed with the later used RescaleIntercept and RescaleSlope. See also calculation for Threshold_HU_Value_to_double.
        Min_HU_for_DensMod_in_double=((2995+1000)/1);

        % Real HU values in row 3
        DensMod(:,3) = (Min_HU_for_DensMod_in_double+1:Min_HU_for_DensMod_in_double+length(DensMod))' - 1000;

        % Use 50 times as many samples since samples are being scrapped in the next step
        numOfSamples = 50 * length(lungIdx);

        % Sample density modulation
        zz1 = ceil(rand(numOfSamples,1)*length(DensMod));
        zz2 = rand(numOfSamples,1);
        ix = zz2 <= DensMod(zz1,2);
        newDist = zz1(ix);

        % Write samples in HU to CT
        ct.cubeHU{1}(lungIdx) = DensMod(newDist(1:numel(lungIdx)),3);

        % Write density cube in CT
        %     ct.cube{1}(lungIdx) = (ct.cubeHU{1}(lungIdx)+1)/1024;

        % Descrete sampling of the density distribution
        %     P = [0; cumsum(DensMod(:,2))];
        %     samples = discretize(rand(numel(lungIdx),1),P);
        %     ct.cube{1}(lungIdx) = samples / max(samples);
end

% In case of Monte Carlo simulations, the Hounsfield Units have to be adjusted to work with the respective
% density to material converters.
if any(contains(pln.propHeterogeneity.sampling.mode,{'TOPAS','MCsquare'}))
    % Only include different densities that are significantly different (1e-3)
    % This is done to significantly increase the computation time
    lung = round(ct.cube{1}(lungIdx),3);

    % Set minimum air density to 0.001225 to avoid deviding by 0
    lung(lung < 0.001225) = 0.001225;

    % Extract distinct densities and respective number of occurences
    [numOfOccurences,sampledDensities] = groupcounts(lung);

    % 
    switch pln.propMC.materialConverter.addSection
        % Modes used to include samples in material converter in case of discrete sampling (either 0 or 1.05)
        case 'lung'
            ct.cubeHU{1}(lungIdx(lung == 1.05)) = 0;
            ct.cubeHU{1}(lungIdx(lung == 0))    = -999;
        case 'poisson'
            ct.cubeHU{1}(lungIdx(lung == 1.05)) = 3020;
            ct.cubeHU{1}(lungIdx(lung == 0))    = 2997;
        % Main mode used by TOPAS and MCsquare pipeline to include samples in material converter
        case 'sampledDensities'
            % Sort densities
            [~,sortIdx] = sort(lung);
            
            % Set individual "virtual" Hounsfield Units (beginning at 6000) for each density and repeat for their number of occurence 
            lungDensitiesNewSorted = repelem(6000:5999+numel(numOfOccurences),1,numOfOccurences);
            
            % Write sorted lung densities to HU cube (using also sorted indices)
            ct.cubeHU{1}(lungIdx(sortIdx)) = lungDensitiesNewSorted;
            
            % Save sampled densities and indices to CT to write in material converter
            ct.sampledDensities   = sampledDensities;
            ct.sampledLungIndices = lungIdx(sortIdx);
        otherwise
            matRad_cfg.dispWarning('Lung modulation should be used with a separate section in the Schneider converter.\n');
    end
end

% Set flag to indicate that the CT has been modulated
ct.modulated = 1;

% plot histogram of the the lung density distribution
% figure, histogram(ct.cube{1}(lungIdx))
end