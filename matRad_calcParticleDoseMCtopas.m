function dij = matRad_calcParticleDoseMCtopas(ct,stf,pln,cst,calcDoseDirect)
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% matRad TOPAS Monte Carlo proton dose calculation wrapper
%   This calls a TOPAS installation (not included in matRad due to
%   licensing model of TOPAS) for MC simulation
%
% call
%   dij = matRad_calcParticleDoseMCtopas(ct,stf,pln,cst,calcDoseDirect)
%
% input
%   ct:                         matRad ct struct
%   stf:                        matRad steering information struct
%   pln:                        matRad plan meta information struct
%   cst:                        matRad cst struct
%   calcDoseDirect:             binary switch to enable forward dose
%                               calcualtion
% output
%   dij:                        matRad dij struct
%
% References
%
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Copyright 2019 the matRad development team.
%
% This file is part of the matRad project. It is subject to the license
% terms in the LICENSE file found in the top-level directory of this
% distribution and at https://github.com/e0404/matRad/LICENSES.txt. No part
% of the matRad project, including this file, may be copied, modified,
% propagated, or distributed except according to the terms contained in the
% LICENSE file.
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Instance of MatRad_Config class
matRad_cfg = MatRad_Config.instance();

% handle inputs
if nargin < 5
    calcDoseDirect = false;
end
if ~calcDoseDirect
    matRad_cfg.dispWarning('You have selected TOPAS dij calculation, this may take a while ^^');
    pln.propMC.calcDij = true;
end
if ~isfield(pln.propStf,'useRangeShifter')
    pln.propStf.useRangeShifter = false;
end

% set calcMC flag just in case it hasn't been set yet
pln.propMC.calcMC = true;

% load default parameters in case they haven't been set yet
pln = matRad_cfg.loadDefaultParam(pln);

% load TOPAS config from pln or from class
if isfield(pln,'propMC') && isfield(pln.propMC,'config')
    if isa(pln.propMC.config,'MatRad_TopasConfig')
        matRad_cfg.dispInfo('Using given Topas Configuration in pln.propMC.config!\n');
        topasConfig = pln.propMC.config;
    else
        % Create a default instance of the configuration
        topasConfig = MatRad_TopasConfig();

        % Overwrite parameters
        % mc = metaclass(topasConfig); %get metaclass information to check if we can overwrite properties
        if isstruct(pln.propMC.config)
            props = fieldnames(pln.propMC.config);
            for fIx = 1:numel(props)
                fName = props{fIx};
                if isprop(topasConfig,fName)
                    % We use a try catch block to catch errors when trying to overwrite protected/private properties
                    % instead of a metaclass approach
                    try
                        topasConfig.(fName) = pln.propMC.config.(fName);
                    catch
                        matRad_cfg.dispWarning('Property ''%s'' for MatRad_TopasConfig will be omitted due to protected/private access or invalid value.',fName);
                    end
                else
                    matRad_cfg.dispWarning('Unkown property ''%s'' for MatRad_TopasConfig will be omitted.',fName);
                end
            end
        else
            matRad_cfg.dispError('Invalid Configuration in pln.propMC.config');
        end
    end
else
    topasConfig = MatRad_TopasConfig();
end

% override default parameters from external parameters if available
if isfield(pln.propHeterogeneity,'sampling') && isfield(pln.propHeterogeneity.sampling,'histories')
    topasConfig.numHistories = pln.propHeterogeneity.sampling.histories;
else
    topasConfig.numHistories = pln.propMC.histories;
end
if isfield(pln.propMC,'numOfRuns')
    topasConfig.numOfRuns = pln.propMC.numOfRuns;
end
if isfield(pln.propMC,'materialConverter') && isfield(pln.propMC.materialConverter,'HUToMaterial')
    topasConfig.materialConverter.HUToMaterial = pln.propMC.materialConverter.HUToMaterial;
end

% set nested folder structure if external calculation is turned on (this will put new simulations in subfolders)
if pln.propMC.externalCalculation
    if isfield(pln,'patientID')
        topasConfig.workingDir = [topasConfig.workingDir pln.radiationMode filesep pln.patientID filesep pln.patientID '_'];
    end
    topasConfig.workingDir = [topasConfig.workingDir pln.radiationMode,'_',pln.machine,'_',datestr(now, 'dd-mm-yy')];
    if isfield(ct,'sampleIdx')
        topasConfig.workingDir = [topasConfig.workingDir '_' num2str(ct.sampleIdx,'%02.f') filesep];
    end
end

%% Initialize dose grid and dij

% load calcDoseInit as usual
matRad_calcDoseInit;

% for TOPAS we explicitly downsample the ct to the dose grid (might not be necessary in future versions with separated grids)
[ctR,~,~] = topasConfig.resampleGrid(ct,cst,pln,stf);

% overwrite CT grid in dij in case of modulation.
if isfield(ctR,'ctGrid')
    dij.ctGrid = ctR.ctGrid;
end

% fill bixels, rays and beams in case of dij calculation or external calculation
if ~calcDoseDirect || pln.propMC.externalCalculation
    counter = 1;
    for f = 1:dij.numOfBeams
        for r = 1:stf(f).numOfRays
            for b = 1:stf(f).numOfBixelsPerRay(r)
                dij.bixelNum(counter) = b;
                dij.rayNum(counter)   = r;
                dij.beamNum(counter)  = f;
                counter = counter + 1;
            end
        end
    end
end

%% sending data to topas

% Load and create TOPAS Base Data
load([pln.radiationMode,'_',pln.machine],'machine');
machine.data = matRad_overrideBaseData(machine.data);

%  Collect weights
if calcDoseDirect
    w = zeros(sum([stf(:).totalNumOfBixels]),ctR.numOfCtScen);
    counter = 1;
    for i = 1:length(stf)
        for j = 1:stf(i).numOfRays
            rayBix = stf(i).numOfBixelsPerRay(j);
            w(counter:counter+rayBix-1,:) = stf(i).ray(j).weight;
            counter = counter + rayBix;
        end
    end
end

if isfield(pln,'bioParam') && strcmp(pln.bioParam.quantityOpt,'RBExD')
    topasConfig.scorer.RBE = true;
    [dij.ax,dij.bx] = matRad_getPhotonLQMParameters(cst,dij.doseGrid.numOfVoxels,1,VdoseGrid);
    dij.abx(dij.bx>0) = dij.ax(dij.bx>0)./dij.bx(dij.bx>0);
end

% save current directory to revert back to later
currDir = cd;

for shiftScen = 1:pln.multScen.totNumShiftScen

    % manipulate isocenter
    for k = 1:length(stf)
        stf(k).isoCenter = stf(k).isoCenter + pln.multScen.isoShift(shiftScen,:);
    end

    % Run simulations for each scenario
    for ctScen = 1:pln.multScen.numOfCtScen
        for rangeShiftScen = 1:pln.multScen.totNumRangeScen
            if pln.multScen.scenMask(ctScen,shiftScen,rangeShiftScen)

                %Overwrite CT (TEMPORARY - we should use 4D calculation in
                %TOPAS here)
                %                 ctR.cubeHU = cubeHUresampled(ctScen);
                %                 ctR.cube = cubeResampled(ctScen);
                % Delete previous topas files so there is no mix-up
                files = dir([topasConfig.workingDir,'*']);
                files = {files(~[files.isdir]).name};
                fclose('all');
                for i = 1:length(files)
                    delete([topasConfig.workingDir,files{i}])
                end

                % actually write TOPAS files
                if calcDoseDirect
                    topasConfig.writeAllFiles(ctR,pln,stf,machine,w(:,ctScen));
                else
                    topasConfig.writeAllFiles(ctR,pln,stf,machine);
                end

                % change director back to original directory
                cd(topasConfig.workingDir);

                % save dij and weights, they are needed for later reading the data back in
                if pln.propMC.externalCalculation
                    save('dij.mat','dij')
                    save('weights.mat','w')
                    matRad_cfg.dispInfo('TOPAS simulation skipped for external calculation\n');
                else
                    for beamIx = 1:numel(stf)
                        for runIx = 1:topasConfig.numOfRuns
                            fname = sprintf('%s_field%d_run%d',topasConfig.label,beamIx,runIx);
                            if isfield(pln.propMC,'verbosity') && strcmp(pln.propMC.verbosity,'full')
                                topasCall = sprintf('%s %s.txt',topasConfig.topasExecCommand,fname);
                            else
                                topasCall = sprintf('%s %s.txt > %s.out > %s.log',topasConfig.topasExecCommand,fname,fname,fname);
                            end

                            % initiate parallel runs and delete previous files
                            if topasConfig.parallelRuns
                                finishedFiles{runIx} = sprintf('%s.finished',fname);
                                if isfile(finishedFiles{runIx})
                                    delete(finishedFiles{runIx});
                                end
                                topasCall = [topasCall '; touch ' finishedFiles{runIx} ' &'];
                            end

                            % Actual simulation happening here
                            matRad_cfg.dispInfo('Calling TOPAS: %s\n',topasCall);
                            [status,cmdout] = system(topasCall,'-echo');

                            % Process TOPAS output and potential errors
                            cout = splitlines(string(cmdout));
                            if status == 0
                                matRad_cfg.dispInfo('TOPAS simulation completed succesfully\n');
                            else
                                if status == 139
                                    matRad_cfg.dispError('TOPAS segmentation fault: might be caused from an outdated TOPAS version or Linux distribution');
                                else
                                    matRad_cfg.dispError('TOPAS simulation exited with error code %d\n "%s"',status,cout(2:end-1));
                                end
                            end
                        end

                        % wait for parallel runs to finish and process
                        if topasConfig.parallelRuns
                            runsFinished = false;
                            pause('on');
                            while ~runsFinished
                                pause(1);
                                fin = cellfun(@(f) exist(f,'file'),finishedFiles);
                                runsFinished = all(fin);
                            end
                        end

                    end

                    % revert back to original directory
                    cd(currDir);

                    %% Simulation finished - read out volume scorers from topas simulation
                    if calcDoseDirect
                        topasCubes = topasConfig.read(topasConfig.workingDir);
                    else
                        topasCubes = topasConfig.read(topasConfig.workingDir,dij);
                    end

                    % save fieldnames to dij for later use and debug
                    fnames = fieldnames(topasCubes);
                    dij.MC_tallies = fnames;

                    % Process read out topasCubes into matRad dij format
                    if calcDoseDirect
                        if any(contains(fnames,'physicalDose'))
                            for d = 1:length(stf)
                                dij.physicalDose{ctScen,1}(:,d)    = sum(w)*reshape(topasCubes.(['physicalDose_beam',num2str(d)]),[],1);
                            end
                        end
                        if any(contains(fnames,'doseToWater'))
                            for d = 1:length(stf)
                                dij.doseToWater{ctScen,1}(:,d)    = sum(w)*reshape(topasCubes.(['doseToWater_beam',num2str(d)]),[],1);
                            end
                        end
                        if any(contains(fnames,'alpha'))
                            abFields = fnames(contains(fnames,{'alpha','beta'}));
                            for ab = 1:length(abFields)
                                model = strsplit(abFields{ab},'_');
                                models{ab} = char(model(2));
                            end
                            models = unique(models);

                            for m = 1:length(models)
                                for d = 1:length(stf)
                                    dij.(['alpha_' models{m}]){ctScen,1}(:,d)           = reshape(topasCubes.(['alpha_' models{m} '_beam',num2str(d)]),[],1);
                                    dij.(['beta_' models{m}]){ctScen,1}(:,d)           = reshape(topasCubes.(['beta_' models{m} '_beam',num2str(d)]),[],1);

                                    dij.(['mAlphaDose_' models{m}]){ctScen,1}(:,d)      = dij.physicalDose{ctScen,1}(:,d) .* dij.(['alpha_' models{m}]){ctScen,1}(:,d);
                                    dij.(['mSqrtBetaDose_' models{m}]){ctScen,1}(:,d)   = sqrt(dij.physicalDose{ctScen,1}(:,d)) .* dij.(['beta_' models{m}]){ctScen,1}(:,d);
                                end
                            end
                        end
                        if any(contains(fnames,'physicalDose_std'))
                            for d = 1:length(stf)
                                dij.physicalDose_std{ctScen,1}(:,d)    = sum(w)*reshape(topasCubes.(['physicalDose_std_beam',num2str(d)]),[],1);
                            end
                        end
                        if any(contains(fnames,'LET'))
                            for d = 1:length(stf)
                                dij.LET{ctScen,1}(:,d)    = reshape(topasCubes.(['LET_beam',num2str(d)]),[],1);
                                dij.mLETDose{ctScen,1}(:,d) = dij.physicalDose{ctScen,1}(:,d) .*  dij.LET{ctScen,1}(:,d);
                            end
                        end
                    else
                        for f = 1:numel(fnames)
                            for d = 1:stf(f).totalNumOfBixels
                                dij.physicalDose{1}(:,d) = reshape(topasCubes.(fnames{f}){d},[],1);
                            end
                        end
                    end
                end
            end
        end
    end
end

% manipulate isocenter back
for k = 1:length(stf)
    stf(k).isoCenter = stf(k).isoCenter - pln.multScen.isoShift(shiftScen,:);
end
end
