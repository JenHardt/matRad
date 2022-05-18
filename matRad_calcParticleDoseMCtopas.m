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

global matRad_cfg;
matRad_cfg = MatRad_Config.instance();
pln.propMC.calcMC = true;
pln = matRad_cfg.loadDefaultParam(pln);

if nargin < 5
    calcDoseDirect = false;
end

if isfield(pln,'propMC') && isfield(pln.propMC,'config')        
    if isa(pln.propMC.config,'MatRad_TopasConfig')
        matRad_cfg.dispInfo('Using given Topas Configuration in pln.propMC.config!\n');
        topasConfig = pln.propMC.config;
    else 
        %Create a default instance of the configuration
        topasConfig = MatRad_TopasConfig();
        
        %Overwrite parameters
        %mc = metaclass(topasConfig); %get metaclass information to check if we can overwrite properties
        
        if isstruct(pln.propMC.config)
            props = fieldnames(pln.propMC.config);
            for fIx = 1:numel(props)
                fName = props{fIx};
                if isprop(topasConfig,fName)
                    %We use a try catch block to catch errors when trying
                    %to overwrite protected/private properties instead of a
                    %metaclass approach
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
        

if ~calcDoseDirect
    matRad_cfg.dispWarning('You have selected TOPAS dij calculation, this may take a while ^^');
    pln.propMC.calcDij = true;
end

if ~isfield(pln.propStf,'useRangeShifter') 
    pln.propStf.useRangeShifter = false;
end

env = matRad_getEnvironment();

%% Initialize dose Grid as usual
% for TOPAS we explicitly downsample the ct to the dose grid (might not
% be necessary in future versions with separated grids)
matRad_calcDoseInit;
[ctR,~,~] = matRad_resampleTopasGrid(ct,cst,pln,stf);
% overwrite CT grid in dij in case of modulation.
if isfield(ctR,'ctGrid')
    dij.ctGrid = ctR.ctGrid;
end

% fill bixels, rays and beams in case of dij calculation
%if ~calcDoseDirect
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
%end

%% sending data to topas

load([pln.radiationMode,'_',pln.machine],'machine');
topasConfig = MatRad_TopasConfig();
% Create Base Data
topasConfig.radiationMode = stf.radiationMode;

machine.data = matRad_overrideBaseData(machine.data);
topasBaseData = MatRad_TopasBaseData(machine,stf);%,TopasConfig);

if isfield(pln.propHeterogeneity,'sampling') && isfield(pln.propHeterogeneity.sampling,'histories')
    topasConfig.numHistories = pln.propHeterogeneity.sampling.histories;
else
    topasConfig.numHistories = pln.propMC.histories;
end
if isfield(pln.propMC,'numOfRuns')
    topasConfig.numOfRuns = pln.propMC.numOfRuns;
end
if pln.propMC.externalCalculation
    if isfield(pln,'patientID')
        topasConfig.workingDir = [topasConfig.workingDir pln.radiationMode filesep pln.patientID filesep pln.patientID '_'];
    end
    topasConfig.workingDir = [topasConfig.workingDir pln.radiationMode,'_',pln.machine,'_',datestr(now, 'dd-mm-yy')];
    if isfield(ctR,'sampleIdx')
        topasConfig.workingDir = [topasConfig.workingDir '_' num2str(ctR.sampleIdx,'%02.f') filesep];
    end
end
% topasConfig.numOfRuns = matRad_cfg.propMC.topas_defaultNumBatches;

%Collect weights
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

currDir = cd;

for shiftScen = 1:pln.multScen.totNumShiftScen
    
    % manipulate isocenter
    for k = 1:length(stf)
        stf(k).isoCenter = stf(k).isoCenter + pln.multScen.isoShift(shiftScen,:);
    end    
    
    for ctScen = 1:pln.multScen.numOfCtScen
        for rangeShiftScen = 1:pln.multScen.totNumRangeScen
            if pln.multScen.scenMask(ctScen,shiftScen,rangeShiftScen)
                
                %Overwrite CT (TEMPORARY - we should use 4D calculation in
                %TOPAS here)
%                 ctR.cubeHU = cubeHUresampled(ctScen);
%                 ctR.cube = cubeResampled(ctScen);
                %Delete previous topas files
                files = dir([topasConfig.workingDir,'*']);
                files = {files(~[files.isdir]).name};
                fclose('all');
                for i = 1:length(files)
                    delete([topasConfig.workingDir,files{i}])
                end
               
                if calcDoseDirect
                    topasConfig.writeAllFiles(ctR,pln,stf,topasBaseData,w(:,ctScen));
                else
                    topasConfig.writeAllFiles(ctR,pln,stf,topasBaseData);
                end
                
                % Run simulation for current scenario
                cd(topasConfig.workingDir);
                              
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

                            if topasConfig.parallelRuns
                                finishedFiles{runIx} = sprintf('%s.finished',fname);
                                delete(finishedFiles{runIx});
                                topasCall = [topasCall '; touch ' finishedFiles{runIx} ' &'];
                            end
                            
                            matRad_cfg.dispInfo('Calling TOPAS: %s\n',topasCall);
                            [status,cmdout] = system(topasCall,'-echo');
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
                    
                    cd(currDir);
                    
                    %% Simulation finished - read out volume scorers from topas simulation
                    if calcDoseDirect
                        topasCubes = matRad_readTopasData(topasConfig.workingDir);
                    else
                        topasCubes = matRad_readTopasData(topasConfig.workingDir,dij);
                    end
                    
                    fnames = fieldnames(topasCubes);
                    dij.MC_tallies = fnames;
                    
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

                                    [dij.ax,dij.bx] = matRad_getPhotonLQMParameters(cst,dij.doseGrid.numOfVoxels,1,VdoseGrid);
                                    dij.abx(dij.bx>0) = dij.ax(dij.bx>0)./dij.bx(dij.bx>0);

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
