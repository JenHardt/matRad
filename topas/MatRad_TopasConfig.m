classdef MatRad_TopasConfig < handle
    % MatRad_TopasConfig class definition
    %
    %
    % References
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
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
    
    properties
        matRad_cfg = MatRad_Config.instance(); %Instance of matRad configuration class
        
        topasExecCommand; %Defaults will be set during construction according to TOPAS installation instructions and used system
        
        parallelRuns = false; %Starts runs in parallel
        
        workingDir; %working directory for the simulation
        
        label = 'matRad_plan';
        
        %Simulation parameters
        numThreads = 0; %number of used threads, 0 = max number of threads (= num cores)
        numOfRuns = 5; %Default number of runs / batches
        modeHistories = 'num'; %'frac';
        fracHistories = 1e-4; %Fraction of histories to compute
        numHistories = 1e6; %Number of histories to compute
        verbosity = struct( 'timefeatures',0,...
            'cputime',true,...
            'run',0,...
            'event',0,...
            'tracking',0,...
            'material',0,...
            'maxinterruptedhistories',1000);
        
        minRelWeight = .00001; %Threshold for discarding beamlets. 0 means all weights are being considered, can otherwise be assigned to min(w)
        
        useOrigBaseData = false; % base data of the original matRad plan will be used?
        beamProfile = 'biGaussian'; %'biGaussian' (emittance); 'simple'
        useEnergySpectrum = false;
        
        %Not yet implemented
        %beamletMode = false; %In beamlet mode simulation will be performed for a dose influence matrix (i.e., each beamlet simulates numHistories beamlets)
        
        pencilBeamScanning = true; %This should be always true (enables deflection)
        
        
        
        %Image
        materialConversion = 'HUToWaterSchneider'; %'CustomWaterRSP';
        %materialConversion = 'CustomWaterRSP';
        arrayOrdering = 'F'; %'C';
        rsp_basematerial = 'Water';
        rsp_vecLength = 10000; %Bins for RSP values when using RSP conversion with custom image converter
        rsp_methodCube = 2; %1: TsBox with variable voxels, 2: TsImageCube with Density Bins and Custom Image converter
        
        
        
        
        %Scoring
        addVolumeScorers = true
        scoreDose = true;
        scoreTrackCount = false;
        scoreDij = false;
        scoreRBE = false;
        %scoreLET = true;
        %         scoreReportQuantity = {'Sum','Standard_Deviation'};
        scoreReportQuantity = 'Sum';
        outputType = 'binary'; %'csv'; 'binary';%
        bioParam = struct( 'PrescribedDose',2,...
            'AlphaX',0.1,...
            'BetaX',0.05);       
        
        %Physics
        electronProductionCut = 0.5; %in mm
        radiationMode;
        modules_protons     = {'g4em-standard_opt4','g4h-phy_QGSP_BIC_HP','g4decay','g4h-elastic_HP','g4stopping','g4ion-QMD','g4radioactivedecay'};
        modules_GenericIon  = {'g4em-standard_opt4','g4h-phy_QGSP_BIC_HP','g4decay','g4h-elastic_HP','g4stopping','g4ion-QMD','g4radioactivedecay'};
        modules_gamma       = {'g4em-standard_opt4','g4h-phy_QGSP_BIC_HP','g4decay'};
        
        %Geometry / World
        worldMaterial = 'G4_AIR';
        
        %filenames
        converterFolder = 'materialConverter';
        scorerFolder = 'scorer';
        outfilenames = struct(  'patientParam','matRad_cube.txt',...
            'patientCube','matRad_cube.dat');
        
        infilenames = struct(   'geometry','TOPAS_matRad_geometry.txt.in',...
            'matConv_Schneider_mod','TOPAS_materialConverter_SchneiderModulation.txt.in',...
            'matConv_Schneider_custom','TOPAS_materialConverter_Schneider_CustomLung.txt.in',...
            'surfaceScorer','TOPAS_scorer_surfaceIC.txt.in',...
            'doseScorer','TOPAS_scorer_dose.txt.in',...
            'doseScorerRBE','TOPAS_scorer_doseRBE_LEM1.txt.in',...
            'doseScorerRBE_MCN','TOPAS_scorer_doseRBE_McNamara.txt.in');
        
        
    end
    
    properties (SetAccess = private)
        thisFolder;
        
        MCparam; %Struct with parameters of last simulation to be saved to file
    end
    
    methods
        function obj = MatRad_TopasConfig()
            %MatRad_MCsquareConfig Construct configuration Class for TOPAS
            %   Default execution paths are set here
            
            obj.thisFolder = fileparts(mfilename('fullpath'));
            obj.workingDir = ['E:/Paper/results/patients/narrowSOPB/' filesep];
%             obj.workingDir = [obj.thisFolder filesep 'MCrun' filesep];
            
            %Let's set some default commands taken from topas installation
            %instructions for mac & debain/ubuntu
            if ispc %We assume topas is installed in wsl (since no windows version)
                obj.topasExecCommand = 'wsl export TOPAS_G4_DATA_DIR=~/G4Data; ~/topas/bin/topas';
            elseif ismac
                obj.topasExecCommand = 'export TOPAS_G4_DATA_DIR=/Applications/G4Data; export QT_QPA_PLATFORM_PLUGIN_PATH=/Applications/topas/Frameworks; /Applications/topas/bin/topas';
            elseif isunix
                obj.topasExecCommand = 'export TOPAS_G4_DATA_DIR=~/G4Data; ~/topas/bin/topas';
            else
                obj.topasExecCommand = '';
            end
            
        end
        
        
        
        function writeAllFiles(obj,ct,pln,stf,topasBaseData,w)
            
            %Reset MCparam structure
            obj.MCparam = struct();
            obj.MCparam.tallies = {};
            
            if isfield(pln,'prescribedDose')
                obj.bioParam.PrescribedDose = pln.prescribedDose;
            end
            if obj.scoreRBE
                obj.radiationMode = stf.radiationMode;
                if all(isfield(pln.propMC,{'AlphaX','BetaX'}))
                   obj.bioParam.AlphaX = pln.propMC.AlphaX;
                   obj.bioParam.BetaX = pln.propMC.BetaX;                
                else
                   for i = 1:length(pln.bioParam.AvailableAlphaXBetaX)
                       if contains(pln.bioParam.AvailableAlphaXBetaX{i,2},'default')
                           break
                       end
                   end
                   obj.bioParam.AlphaX = pln.bioParam.AvailableAlphaXBetaX{5,1}(1);
                   obj.bioParam.BetaX = pln.bioParam.AvailableAlphaXBetaX{5,1}(2);
                end
            end
            if isfield(pln.propMC,'calcDij') && pln.propMC.calcDij
                obj.scoreDij = true;
            end

            obj.MCparam.nbRuns = obj.numOfRuns;
            obj.MCparam.simLabel = obj.label;
            obj.MCparam.scoreReportQuantity = obj.scoreReportQuantity;
            
            
            if ~exist(obj.workingDir,'dir')
                mkdir(obj.workingDir);
                obj.matRad_cfg.dispInfo('Created TOPAS working directory %s\n',obj.workingDir);
            end
            
            obj.MCparam.workingDir = obj.workingDir;
            
            obj.matRad_cfg.dispInfo('Writing parameter files to %s\n',obj.workingDir);
            
            obj.writePatient(ct,pln);
            if ~exist('w','var')
                numBixels = sum([stf(:).totalNumOfBixels]);
                w = ones(numBixels,1);
            end
            obj.writeStfFields(ct,stf,topasBaseData,w);
            
            obj.matRad_cfg.dispInfo('Successfully written TOPAS setup files!\n')
            
            obj.writeMCparam();
        end
        
    end
    
    %Private sub functions for writing (private so the state of the configuration
    %can not be corrupted)
    methods (Access = private)
        
        function writeRunHeader(obj,fID,fieldIx,runIx)
            
            fprintf(fID,'s:Sim/PlanLabel = "%s"\n',obj.label);
            fprintf(fID,'s:Sim/ScoreLabel = "score_%s_field%d_run%d"\n',obj.label,fieldIx,runIx);
            fprintf(fID,'\n');
            
            logicalString = {'"False"', '"True"'};
            
            fprintf(fID,'i:Ma/Verbosity = %d\n',obj.verbosity.material);
            fprintf(fID,'i:Ts/TrackingVerbosity = %d\n',obj.verbosity.tracking);
            fprintf(fID,'i:Ts/EventVerbosity = %d\n',obj.verbosity.event);
            fprintf(fID,'i:Ts/RunVerbosity = %d\n',obj.verbosity.run);
            fprintf(fID,'b:Ts/ShowCPUTime = %s\n',logicalString{obj.verbosity.cputime + 1});
            fprintf(fID,'i:Tf/Verbosity = %d\n',obj.verbosity.timefeatures);
            fprintf(fID,'i:Ts/MaxInterruptedHistories = %d\n',obj.verbosity.maxinterruptedhistories);
            fprintf(fID,'Ts/NumberOfThreads = %d\n',obj.numThreads);
            fprintf(fID,'i:Ts/ShowHistoryCountAtInterval = %d\n',10^(floor(log10(1/obj.numOfRuns * obj.numHistories))-1));
            fprintf(fID,'\n');
            
            
            fprintf(fID,'s:Sim/DoseScorerOutputType = "%s"\n',obj.outputType);
            if iscell(obj.scoreReportQuantity)
                fprintf(fID,'sv:Sim/DoseScorerReport = %i ',length(obj.scoreReportQuantity));
                fprintf(fID,'"%s" ',obj.scoreReportQuantity{:});
                fprintf(fID,'\n');
            else
                fprintf(fID,'sv:Sim/DoseScorerReport = 1 "%s"\n',obj.scoreReportQuantity);
            end
            fprintf(fID,'\n');
            
            %fprintf(fID,'includeFile = %s/TOPAS_Simulation_Setup.txt\n',obj.thisFolder);
            %fprintf(fID,'includeFile = %s/TOPAS_matRad_geometry.txt\n',obj.thisFolder);
            %fprintf(fID,'includeFile = %s/TOPAS_scorer_surfaceIC.txt\n',obj.thisFolder);
        end
        
        function writeFieldHeader(obj,fID,fieldIx)
            fprintf(fID,'u:Sim/HalfValue = %d\n',0.5);
            fprintf(fID,'u:Sim/SIGMA2FWHM = %d\n',2.354818);
            fprintf(fID,'u:Sim/FWHM2SIGMA = %d\n',0.424661);
            fprintf(fID,'\n');
            
            fprintf(fID,'d:Sim/ElectronProductionCut = %f mm\n',obj.electronProductionCut);
            fprintf(fID,'s:Sim/WorldMaterial = "%s"\n',obj.worldMaterial);
            fprintf(fID,'\n');
            
            fprintf(fID,'includeFile = %s\n',obj.outfilenames.patientParam);
            fprintf(fID,'\n');
            
            fname = fullfile(obj.thisFolder,obj.infilenames.geometry);
            obj.matRad_cfg.dispInfo('Reading Geometry from %s\n',fname);
            world = fileread(fname);
            fprintf(fID,'%s\n\n',world);
        end
        
        function writeScorers(obj,fID)
            
            obj.MCparam.outputType = obj.outputType;
            
            if obj.scoreDose
                if obj.scoreRBE
                    if strcmp(obj.radiationMode,'protons')
                        fname = fullfile(obj.thisFolder,filesep,obj.scorerFolder,filesep,obj.infilenames.doseScorerRBE_MCN);
                    elseif strcmp(obj.radiationMode,'carbon') || strcmp(obj.radiationMode,'helium')
                        fname = fullfile(obj.thisFolder,filesep,obj.scorerFolder,filesep,obj.infilenames.doseScorerRBE);
                    else
                        obj.matRad_cfg.dispWarning('No specific RBE model available, using custom LEM scorer.');
                        fname = fullfile(obj.thisFolder,filesep,obj.scorerFolder,filesep,obj.infilenames.doseScorerRBE);
                        
                    end
                else
                    fname = fullfile(obj.thisFolder,filesep,obj.scorerFolder,filesep,obj.infilenames.doseScorer);
                end
                
                obj.matRad_cfg.dispDebug('Writing Biologial Scorer components.\n');
                fprintf(fID,'d:Sc/PrescribedDose = %.4f Gy\n',obj.bioParam.PrescribedDose);
                fprintf(fID,'d:Sc/AlphaX = %.4f /Gy\n',obj.bioParam.AlphaX);
                fprintf(fID,'d:Sc/BetaX = %.4f /Gy2\n',obj.bioParam.BetaX);
                fprintf(fID,'d:Sc/AlphaBetaX = %.4f Gy\n',obj.bioParam.AlphaX/obj.bioParam.BetaX);
                
                obj.matRad_cfg.dispDebug('Reading Dose Scorer from %s\n',fname);
                scorer = fileread(fname);
                fprintf(fID,'%s\n',scorer);
                if obj.scoreRBE
                    %                    obj.MCparam.tallies = {'physicalDose','alpha','beta','RBE'};
                    obj.MCparam.tallies = {'physicalDose','alpha','beta'};
                end
                if ~any(strcmp(obj.MCparam.tallies,'physicalDose'))
                    obj.MCparam.tallies{end+1} = 'physicalDose';
                end               
            end
            
            if obj.addVolumeScorers
                fileList = dir(fullfile(obj.thisFolder,filesep,obj.scorerFolder,filesep,'TOPAS_scorer_volume_*.in'));
                for fileIx=1:length(fileList)
                    fname = fullfile(obj.thisFolder,fileList(fileIx).name);
                    obj.matRad_cfg.dispDebug('Reading Volume Scorer from %s\n',fname);
                    scorer = fileread(fname);
                    fprintf(fID,'%s\n',scorer);
                    tallyLabel = regexprep(fileList(fileIx).name,'TOPAS_scorer_volume_','');
                    tallyLabel = regexprep(tallyLabel,'.txt.in','');
                    obj.MCparam.tallies{end+1} = tallyLabel;
                end
            end
            
            if obj.scoreTrackCount
                fname = fullfile(obj.thisFolder,filesep,obj.scorerFolder,filesep,obj.infilenames.surfaceScorer);
                obj.matRad_cfg.dispDebug('Reading surface scorer from %s\n',fname);
                scorer = fileread(fname);
                fprintf(fID,'%s\n',scorer);
            end
        end
        
        function writeStfFields(obj,ct,stf,baseData,w)
            
            %Bookkeeping
            obj.MCparam.nbFields = length(stf);
            
            %Sanity check
            if numel(w) ~= sum([stf(:).totalNumOfBixels])
                obj.matRad_cfg.dispError('Given number of weights (#%d) doesn''t match bixel count in stf (#%d)',numel(w), sum([stf(:).totalNumOfBixels]));
            end
            
            nozzleToAxisDistance = baseData.nozzleToIso;
            
            nParticlesTotalBixel = round(1e6*w);
            %nParticlesTotal = sum(nParticlesTotalBixel);
            maxParticlesBixel = 1e6*max(w(:));
            minParticlesBixel = round(max([obj.minRelWeight*maxParticlesBixel,1]));
            
            switch obj.modeHistories
                case 'num'
                    obj.fracHistories = obj.numHistories ./ sum(nParticlesTotalBixel);
                case 'frac'
                    obj.numHistories = sum(nParticlesTotalBixel);
                otherwise
                    obj.matRad_cfg.dispError('Invalid history setting!');
            end
            
            nParticlesTotal = 0;
            
            for beamIx = 1:length(stf)
                
                SAD = stf(beamIx).SAD;
                
                sourceToNozzleDistance = SAD - nozzleToAxisDistance;
                
                %Selection of base data given the energies
                if obj.useOrigBaseData
                    [~,ixTmp,~] = intersect([ baseData.machine.data.energy], [stf.ray.energy]);
                    for i = 1:length(ixTmp)
                        selectedData(i) =  baseData.machine.data(ixTmp(i));
                    end
                    energies = [selectedData.energy];
                else
                    selectedData = [];
                    focusIndex = baseData.selectedFocus(baseData.energyIndex);
                    for i = 1:numel(focusIndex)
                        selectedData = [selectedData, baseData.monteCarloData(focusIndex(i), i)];
                    end
                    energies = [selectedData.NominalEnergy];
                end
                
                %Get Range Shifters in field ifpresent
                allRays = [stf(beamIx).ray];
                raShis = [allRays.rangeShifter];
                [~,ix] =  unique(cell2mat(squeeze(struct2cell(raShis))'),'rows');
                
                raShis = raShis(ix);
                ix = [raShis.ID] == 0;
                raShis = raShis(~ix);
                
                %Convert ID into readable string
                for r = 1:numel(raShis)
                    if isnumeric(raShis(r).ID)
                        raShis(r).topasID = ['RangeShifter' num2str(raShis(r).ID)];
                    else
                        raShis(r).topasID = ['RangeShifter' raShis(r).ID];
                    end
                end
                
                
                
                %get beamlet properties for each bixel in the stf and write
                %it into dataTOPAS
                currentBixel = 1;
                cutNumOfBixel = 0;
                nBeamParticlesTotal(beamIx) = 0;
                
                collectBixelIdx = [];
                
                %Loop over rays and then over spots on ray
                for rayIx = 1:stf(beamIx).numOfRays
                    for bixelIx = 1:stf(beamIx).numOfBixelsPerRay(rayIx)
                        
                        nCurrentParticles = nParticlesTotalBixel(currentBixel);
                        
                        % check whether there are (enough) particles for beam delivery
                        if (nCurrentParticles>minParticlesBixel)
                            
                            collectBixelIdx(end+1) = bixelIx;
                            cutNumOfBixel = cutNumOfBixel + 1;
                            bixelEnergy = stf(beamIx).ray(rayIx).energy(bixelIx);
                            [~,ixTmp,~] = intersect(energies, bixelEnergy);
                            
                            voxel_x = -stf(beamIx).ray(rayIx).rayPos_bev(3);
                            voxel_y = stf(beamIx).ray(rayIx).rayPos_bev(1);
                            
                            dataTOPAS(cutNumOfBixel).posX = -1.*voxel_x;
                            dataTOPAS(cutNumOfBixel).posY = voxel_y;
                            
                            if obj.scoreDij
                                % write particles directly to every beamlet for dij calculation (each bixel
                                % calculated separately with full numParticles
                                dataTOPAS(cutNumOfBixel).current = uint32(nCurrentParticles / obj.numOfRuns);
                            else
                                dataTOPAS(cutNumOfBixel).current = uint32(obj.fracHistories*nCurrentParticles / obj.numOfRuns);
                            end
                            
                            if obj.pencilBeamScanning
                                % angleX corresponds to the rotation around the X axis necessary to move the spot in the Y direction
                                % angleY corresponds to the rotation around the Y' axis necessary to move the spot in the X direction
                                % note that Y' corresponds to the Y axis after the rotation of angleX around X axis
                                dataTOPAS(cutNumOfBixel).angleX = atan(dataTOPAS(cutNumOfBixel).posY / SAD);
                                dataTOPAS(cutNumOfBixel).angleY = atan(-dataTOPAS(cutNumOfBixel).posX ./ (SAD ./ cos(dataTOPAS(cutNumOfBixel).angleX)));
                                dataTOPAS(cutNumOfBixel).posX = (dataTOPAS(cutNumOfBixel).posX / SAD)*(SAD-nozzleToAxisDistance);
                                dataTOPAS(cutNumOfBixel).posY = (dataTOPAS(cutNumOfBixel).posY / SAD)*(SAD-nozzleToAxisDistance);
                            end
                            
                            if obj.useOrigBaseData
                                dataTOPAS(cutNumOfBixel).energy = selectedData(ixTmp).energy;
                                dataTOPAS(cutNumOfBixel).focusFWHM = selectedData(ixTmp).initFocus.SisFWHMAtIso(stf(beamIx).ray(rayIx).focusIx(bixelIx));
                                
                            else
                                dataTOPAS(cutNumOfBixel).energy = selectedData(ixTmp).MeanEnergy;
                                dataTOPAS(cutNumOfBixel).nominalEnergy = selectedData(ixTmp).NominalEnergy;
                                dataTOPAS(cutNumOfBixel).energySpread = selectedData(ixTmp).EnergySpread;
                                dataTOPAS(cutNumOfBixel).spotSize = selectedData(ixTmp).SpotSize1x;
                                dataTOPAS(cutNumOfBixel).divergence = selectedData(ixTmp).Divergence1x;
                                dataTOPAS(cutNumOfBixel).correlation = selectedData(ixTmp).Correlation1x;
                                dataTOPAS(cutNumOfBixel).focusFWHM = selectedData(ixTmp).FWHMatIso;
                            end
                            
                            if obj.scoreDij
                                % remember beam and bixel number
                                dataTOPAS(cutNumOfBixel).beam           = beamIx;
                                dataTOPAS(cutNumOfBixel).ray            = rayIx;
                                dataTOPAS(cutNumOfBixel).bixel          = bixelIx;
                                dataTOPAS(cutNumOfBixel).totalBixel     = currentBixel;
                            end
                            
                            %Add RangeShifterState
                            if ~isempty(raShis)
                                for r = 1:length(raShis)
                                    if stf(beamIx).ray(rayIx).rangeShifter(bixelIx).ID == raShis(r).ID
                                        raShiOut(r) = 0; %Range shifter is in beam path
                                    else
                                        raShiOut(r) = 1; %Range shifter is out of beam path / not used
                                    end
                                end
                                dataTOPAS(cutNumOfBixel).raShiOut = raShiOut;
                            end
                            
                            nBeamParticlesTotal(beamIx) = nBeamParticlesTotal(beamIx) + nCurrentParticles;
                            
                            
                        end
                        
                        currentBixel = currentBixel + 1;
                        
                    end
                end
                
                nParticlesTotal = nParticlesTotal + nBeamParticlesTotal(beamIx);
                
                % discard data if the current has unphysical values
                idx = find([dataTOPAS.current] < 1);
                dataTOPAS(idx) = [];
                cutNumOfBixel = length(dataTOPAS(:));
                if obj.scoreDij
                    historyCount(beamIx) = uint32(nBeamParticlesTotal(beamIx) / obj.numOfRuns);
                else
                    historyCount(beamIx) = uint32(obj.fracHistories * nBeamParticlesTotal(beamIx) / obj.numOfRuns);
                end
                
                if historyCount(beamIx) < cutNumOfBixel || cutNumOfBixel == 0
                    obj.matRad_cfg.dispError('Insufficient number of histories!')
                end
                
%                 while sum([dataTOPAS.current]) ~= historyCount(beamIx)
%                     randIx = randi([1 length(dataTOPAS)],1,abs(sum([dataTOPAS(:).current]) - historyCount(beamIx)));
%                     
%                     if (sum([dataTOPAS(:).current]) > historyCount(beamIx))
%                         for i = randIx
%                             dataTOPAS(i).current = dataTOPAS(i).current - 1;
%                         end
%                     else
%                         for i = randIx
%                             dataTOPAS(i).current = dataTOPAS(i).current + 1;
%                         end
%                     end
%                 end
%                 
                %adjust current to actual histories (by adding/subtracting
                %from random rays)
                while sum([dataTOPAS(:).current]) ~= historyCount(beamIx)
                    
                    diff = sum([dataTOPAS.current]) - sum(historyCount(beamIx));                  
                    [~,R] = histc(rand(abs(diff),1),cumsum([0;double(transpose([dataTOPAS(:).current]))./double(sum([dataTOPAS(:).current]))]));           
                    idx = 1:length(dataTOPAS);
                    randIx = idx(R);
                    
                    newCurr = num2cell(arrayfun(@plus,double([dataTOPAS(randIx).current]),-1*sign(diff)*ones(1,abs(diff))),1);
                    [dataTOPAS(randIx).current] = newCurr{:};
                end

                historyCount(beamIx) = historyCount(beamIx) * obj.numOfRuns;
                
                
                %sort dataTOPAS according to energy
                [~,ixSorted] = sort([dataTOPAS(:).energy]);
                dataTOPAS = dataTOPAS(ixSorted);
                
                %write TOPAS data base file
                fieldSetupFileName = sprintf('beamSetup_%s_field%d.txt',obj.label,beamIx);
                fileID = fopen(fullfile(obj.workingDir,fieldSetupFileName),'w');
                obj.writeFieldHeader(fileID,beamIx);
                
                %Write modality specific info
                switch stf(beamIx).radiationMode
                    case 'protons'
                        fprintf(fileID,'s:Sim/ParticleName = "proton"\n');
                        fprintf(fileID,'u:Sim/ParticleMass = 1.0\n');
                        
                        particleA = 1;
                        particleZ = 1;
                        
                        modules = obj.modules_protons;
                        
                    case 'carbon'
                        fprintf(fileID,'s:Sim/ParticleName = "GenericIon(6,12)"\n');
                        fprintf(fileID,'u:Sim/ParticleMass = 12.0\n');
                        
                        particleA = 12;
                        particleZ = 6;
                        
                        modules = obj.modules_GenericIon;
                        
                    case 'helium'
                        fprintf(fileID,'s:Sim/ParticleName = "GenericIon(2,4)"\n');
                        fprintf(fileID,'u:Sim/ParticleMass = 4.0\n');
                        
                        particleA = 4;
                        particleZ = 2;
                        
                        modules = obj.modules_GenericIon;
                        %{
                    case 'photons' %This modality is not yet used!
                           
                        fprintf(fileID,'s:Sim/ParticleName = "gamma"\n');
                        fprintf(fileID,'u:Sim/ParticleMass = 0\n');
                        
                        particleA = 0;
                        particleZ = 0;
                        
                        obj.modules_gamma;
                        %}
                    otherwise
                        obj.matRad_cfg.dispError('Invalid radiation mode %s!',stf.radiationMode)
                end
                
                fprintf(fileID,'d:Sim/GantryAngle = %f deg\n',stf(beamIx).gantryAngle);
                fprintf(fileID,'d:Sim/CouchAngle = %f deg\n',stf(beamIx).couchAngle);
                
                fprintf(fileID,'d:Tf/TimelineStart = 0. ms\n');
                fprintf(fileID,'d:Tf/TimelineEnd = %i ms\n', 10 * cutNumOfBixel);
                fprintf(fileID,'i:Tf/NumberOfSequentialTimes = %i\n', cutNumOfBixel);
                fprintf(fileID,'dv:Tf/Beam/Spot/Times = %i ', cutNumOfBixel);
                fprintf(fileID,num2str(linspace(10,cutNumOfBixel*10,cutNumOfBixel)));
                fprintf(fileID,' ms\n');
                %fprintf(fileID,'uv:Tf/Beam/Spot/Values = %i %s\n',cutNumOfBixel,num2str(collectBixelIdx));
                
                if isfield(baseData.machine.data,'energySpectrum') && obj.useEnergySpectrum
                    
                    obj.matRad_cfg.dispInfo('Beam energy spectrum available\n');
                    energySpectrum = [baseData.machine.data(:).energySpectrum];
                    nbSpectrumPoints = length(energySpectrum(1).energy_MeVpN);
                    
                    [~,energyIx] = ismember([dataTOPAS.nominalEnergy],[baseData.machine.data.energy]);
                    
                    fprintf(fileID,'s:So/PencilBeam/BeamEnergySpectrumType = "Continuous"\n');
                    fprintf(fileID,'dv:So/PencilBeam/BeamEnergySpectrumValues = %d %s MeV\n',nbSpectrumPoints,strtrim(sprintf('Tf/Beam/EnergySpectrum/Energy/Point%03d/Value ',1:nbSpectrumPoints)));
                    fprintf(fileID,'uv:So/PencilBeam/BeamEnergySpectrumWeights = %d %s\n',nbSpectrumPoints,strtrim(sprintf('Tf/Beam/EnergySpectrum/Weight/Point%03d/Value ',1:nbSpectrumPoints)));
                    points_energy = reshape([energySpectrum(energyIx).energy_MeVpN],[],length(energyIx));
                    points_weight = reshape([energySpectrum(energyIx).weight],[],length(energyIx));
                    for spectrumPoint=1:nbSpectrumPoints
                        fprintf(fileID,'s:Tf/Beam/EnergySpectrum/Energy/Point%03d/Function = "Step"\n',spectrumPoint);
                        fprintf(fileID,'dv:Tf/Beam/EnergySpectrum/Energy/Point%03d/Times = Tf/Beam/Spot/Times ms\n',spectrumPoint);
                        fprintf(fileID,'dv:Tf/Beam/EnergySpectrum/Energy/Point%03d/Values = %d %s MeV\n',spectrumPoint,cutNumOfBixel,strtrim(sprintf('%f ',particleA*points_energy(spectrumPoint,:))));
                        fprintf(fileID,'s:Tf/Beam/EnergySpectrum/Weight/Point%03d/Function = "Step"\n',spectrumPoint);
                        fprintf(fileID,'dv:Tf/Beam/EnergySpectrum/Weight/Point%03d/Times = Tf/Beam/Spot/Times ms\n',spectrumPoint);
                        fprintf(fileID,'uv:Tf/Beam/EnergySpectrum/Weight/Point%03d/Values = %d %s\n',spectrumPoint,cutNumOfBixel,strtrim(sprintf('%f ',points_weight(spectrumPoint,:))));
                    end
                    
                end
                
                fprintf(fileID,'s:Tf/Beam/Energy/Function = "Step"\n');
                fprintf(fileID,'dv:Tf/Beam/Energy/Times = Tf/Beam/Spot/Times ms\n');
                fprintf(fileID,'dv:Tf/Beam/Energy/Values = %i ', cutNumOfBixel);
                
                fprintf(fileID,num2str(particleA*[dataTOPAS.energy])); %Transform total energy with atomic number
                fprintf(fileID,' MeV\n');
                
                
                switch obj.beamProfile
                    case 'biGaussian'
                        fprintf(fileID,'s:Tf/Beam/EnergySpread/Function = "Step"\n');
                        fprintf(fileID,'dv:Tf/Beam/EnergySpread/Times = Tf/Beam/Spot/Times ms\n');
                        fprintf(fileID,'uv:Tf/Beam/EnergySpread/Values = %i ', cutNumOfBixel);
                        fprintf(fileID,num2str([dataTOPAS.energySpread]));
                        fprintf(fileID,'\n');
                        
                        fprintf(fileID,'s:Tf/Beam/SigmaX/Function = "Step"\n');
                        fprintf(fileID,'dv:Tf/Beam/SigmaX/Times = Tf/Beam/Spot/Times ms\n');
                        fprintf(fileID,'dv:Tf/Beam/SigmaX/Values = %i ', cutNumOfBixel);
                        fprintf(fileID,num2str([dataTOPAS.spotSize]));
                        fprintf(fileID,' mm\n');
                        fprintf(fileID,'s:Tf/Beam/SigmaXPrime/Function = "Step"\n');
                        fprintf(fileID,'dv:Tf/Beam/SigmaXPrime/Times = Tf/Beam/Spot/Times ms\n');
                        fprintf(fileID,'uv:Tf/Beam/SigmaXPrime/Values = %i ', cutNumOfBixel);
                        fprintf(fileID,num2str([dataTOPAS.divergence]));
                        fprintf(fileID,'\n');
                        fprintf(fileID,'s:Tf/Beam/CorrelationX/Function = "Step"\n');
                        fprintf(fileID,'dv:Tf/Beam/CorrelationX/Times = Tf/Beam/Spot/Times ms\n');
                        fprintf(fileID,'uv:Tf/Beam/CorrelationX/Values = %i ', cutNumOfBixel);
                        fprintf(fileID,num2str([dataTOPAS.correlation]));
                        fprintf(fileID,'\n');
                        
                        fprintf(fileID,'s:Tf/Beam/SigmaY/Function = "Step"\n');
                        fprintf(fileID,'dv:Tf/Beam/SigmaY/Times = Tf/Beam/Spot/Times ms\n');
                        fprintf(fileID,'dv:Tf/Beam/SigmaY/Values = %i ', cutNumOfBixel);
                        fprintf(fileID,num2str([dataTOPAS.spotSize]));
                        fprintf(fileID,' mm\n');
                        fprintf(fileID,'s:Tf/Beam/SigmaYPrime/Function = "Step"\n');
                        fprintf(fileID,'dv:Tf/Beam/SigmaYPrime/Times = Tf/Beam/Spot/Times ms\n');
                        fprintf(fileID,'uv:Tf/Beam/SigmaYPrime/Values = %i ', cutNumOfBixel);
                        fprintf(fileID,num2str([dataTOPAS.divergence]));
                        fprintf(fileID,'\n');
                        fprintf(fileID,'s:Tf/Beam/CorrelationY/Function = "Step"\n');
                        fprintf(fileID,'dv:Tf/Beam/CorrelationY/Times = Tf/Beam/Spot/Times ms\n');
                        fprintf(fileID,'uv:Tf/Beam/CorrelationY/Values = %i ', cutNumOfBixel);
                        fprintf(fileID,num2str([dataTOPAS.correlation]));
                        fprintf(fileID,'\n');
                    case 'simple'
                        fprintf(fileID,'s:Tf/Beam/FocusFWHM/Function = "Step"\n');
                        fprintf(fileID,'dv:Tf/Beam/FocusFWHM/Times = Tf/Beam/Spot/Times ms\n');
                        fprintf(fileID,'dv:Tf/Beam/FocusFWHM/Values = %i ', cutNumOfBixel);
                        fprintf(fileID,num2str([dataTOPAS.focusFWHM]));
                        fprintf(fileID,' mm\n');
                end
                
                if obj.pencilBeamScanning
                    fprintf(fileID,'s:Tf/Beam/AngleX/Function = "Step"\n');
                    fprintf(fileID,'dv:Tf/Beam/AngleX/Times = Tf/Beam/Spot/Times ms\n');
                    fprintf(fileID,'dv:Tf/Beam/AngleX/Values = %i ', cutNumOfBixel);
                    fprintf(fileID,num2str([dataTOPAS.angleX]));
                    fprintf(fileID,' rad\n');
                    fprintf(fileID,'s:Tf/Beam/AngleY/Function = "Step"\n');
                    fprintf(fileID,'dv:Tf/Beam/AngleY/Times = Tf/Beam/Spot/Times ms\n');
                    fprintf(fileID,'dv:Tf/Beam/AngleY/Values = %i ', cutNumOfBixel);
                    fprintf(fileID,num2str([dataTOPAS.angleY]));
                    fprintf(fileID,' rad\n');
                end
                
                fprintf(fileID,'s:Tf/Beam/PosX/Function = "Step"\n');
                fprintf(fileID,'dv:Tf/Beam/PosX/Times = Tf/Beam/Spot/Times ms\n');
                fprintf(fileID,'dv:Tf/Beam/PosX/Values = %i ', cutNumOfBixel);
                fprintf(fileID,num2str([dataTOPAS.posX]));
                fprintf(fileID,' mm\n');
                fprintf(fileID,'s:Tf/Beam/PosY/Function = "Step"\n');
                fprintf(fileID,'dv:Tf/Beam/PosY/Times = Tf/Beam/Spot/Times ms\n');
                fprintf(fileID,'dv:Tf/Beam/PosY/Values = %i ', cutNumOfBixel);
                fprintf(fileID,num2str([dataTOPAS.posY]));
                fprintf(fileID,' mm\n');
                
                fprintf(fileID,'s:Tf/Beam/Current/Function = "Step"\n');
                fprintf(fileID,'dv:Tf/Beam/Current/Times = Tf/Beam/Spot/Times ms\n');
                fprintf(fileID,'iv:Tf/Beam/Current/Values = %i ', cutNumOfBixel);
                fprintf(fileID,num2str([dataTOPAS.current]));
                fprintf(fileID,'\n\n');
                
                %Range shifter in/out
                if ~isempty(raShis)
                    fprintf(fileID,'#Range Shifter States:\n');
                    for r = 1:numel(raShis)
                        fprintf(fileID,'s:Tf/Beam/%sOut/Function = "Step"\n',raShis(r).topasID);
                        fprintf(fileID,'dv:Tf/Beam/%sOut/Times = Tf/Beam/Spot/Times ms\n',raShis(r).topasID);
                        fprintf(fileID,'uv:Tf/Beam/%sOut/Values = %i ', raShis(r).topasID, cutNumOfBixel);
                        fprintf(fileID,num2str([dataTOPAS.raShiOut]));
                        fprintf(fileID,'\n\n');
                    end
                end
                
                
                % NozzleAxialDistance
                fprintf(fileID,'d:Ge/Nozzle/TransZ = -%f mm\n', nozzleToAxisDistance);
                if obj.pencilBeamScanning
                    fprintf(fileID,'d:Ge/Nozzle/RotX = Tf/Beam/AngleX/Value rad\n');
                    fprintf(fileID,'d:Ge/Nozzle/RotY = Tf/Beam/AngleY/Value rad\n');
                    fprintf(fileID,'d:Ge/Nozzle/RotZ = 0.0 rad\n');
                end
                
                %Range Shifter Definition
                for r = 1:numel(raShis)
                    obj.writeRangeShifter(fileID,raShis(r),sourceToNozzleDistance);
                end
                
                switch obj.beamProfile
                    case 'biGaussian'
                        TOPAS_beamSetup = fileread('TOPAS_beamSetup_biGaussian.txt.in');
                    case 'simple'
                        TOPAS_beamSetup = fileread('TOPAS_beamSetup_generic.txt.in');
                end
                
                fprintf(fileID,'%s\n',TOPAS_beamSetup);
                
                %translate patient according to beam isocenter
                fprintf(fileID,'d:Ge/Patient/TransX      = %f mm\n',0.5*ct.resolution.x*(ct.cubeDim(2)+1)-stf(beamIx).isoCenter(1));
                fprintf(fileID,'d:Ge/Patient/TransY      = %f mm\n',0.5*ct.resolution.y*(ct.cubeDim(1)+1)-stf(beamIx).isoCenter(2));
                fprintf(fileID,'d:Ge/Patient/TransZ      = %f mm\n',0.5*ct.resolution.z*(ct.cubeDim(3)+1)-stf(beamIx).isoCenter(3));
                fprintf(fileID,'d:Ge/Patient/RotX=0. deg\n');
                fprintf(fileID,'d:Ge/Patient/RotY=0. deg\n');
                fprintf(fileID,'d:Ge/Patient/RotZ=0. deg\n');
                
                %load topas modules depending on the particle type
                fprintf(fileID,'\n# MODULES\n');
                moduleString = cellfun(@(s) sprintf('"%s"',s),modules,'UniformOutput',false);
                fprintf(fileID,'sv:Ph/Default/Modules = %d %s\n',length(modules),strjoin(moduleString,' '));
                
                fclose(fileID);
                %write run scripts for TOPAS
                for runIx = 1:obj.numOfRuns
                    runFileName = sprintf('%s_field%d_run%d.txt',obj.label,beamIx,runIx);
                    fileID = fopen(fullfile(obj.workingDir,runFileName),'w');
                    obj.writeRunHeader(fileID,beamIx,runIx);
                    fprintf(fileID,['i:Ts/Seed = ',num2str(runIx),'\n']);
                    fprintf(fileID,'includeFile = ./%s\n',fieldSetupFileName);
                    obj.writeScorers(fileID);
                    if obj.scoreDij
                        fprintf(fileID,'s:Tf/ImageName/Function = "Step"\n');
                        % create time feature scorer and save with original rays and bixel names
                        imageName = ['sv:Tf/ImageName/Values = ',num2str(cutNumOfBixel),strcat(' "ray',string([dataTOPAS.ray]))+strcat('_bixel',string([dataTOPAS.bixel]),'"')];
                        fprintf(fileID,'%s\n',strjoin(imageName));
                        fprintf(fileID,'dv:Tf/ImageName/Times = Tf/Beam/Spot/Times ms\n');
                        fprintf(fileID,'s:Sc/Patient/Tally_DoseToMedium/SplitByTimeFeature = "ImageName"');
                    end
                    fclose(fileID);
                end
            end
            
            
            
            
            
            %Bookkeeping
            obj.MCparam.nbParticlesTotal = nParticlesTotal;
            obj.MCparam.nbHistoriesTotal = sum(historyCount);
            obj.MCparam.nbParticlesField = nBeamParticlesTotal;
            obj.MCparam.nbHistoriesField = historyCount;
        end
        
        
        function writePatient(obj,ct,pln)
            % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % matRad export CT RSP data for TOPAS simulation
            %
            % call
            %   matRad_exportCtTOPAS(ct, path, material)
            %
            % input
            %   ct:             ct cube
            %   runsPath:       path where to save the files for MC simulation
            %   basematerial:   base material to be scaled to corresponding RSP
            %
            % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            
            
            % the image cube contains the indexing of materials
            % since at the moment TOPAS does not support ushort
            % the materials should have indexes between 0 and 32767
            % therefore, the maximum length of the vector is 32768
            %answer = ''
            %while isempty(answer)
            %  prompt = 'Length of the imaging vector [2-32768]:';
            %  answer = input(prompt)
            %  if answer > 1 && answer <= 32768
            %    vecRSPlength=answer
            %  else
            %    answer = ''
            %  end
            %end
            
            medium = obj.rsp_basematerial;
            vecRSPlength=obj.rsp_vecLength;
            
            if isequal(obj.arrayOrdering,'C')
                if obj.matRad_cfg.logLevel > 2
                	disp('Exporting cube in C ordering...')
                end
                permutation = [3 1 2];
            else
                if obj.matRad_cfg.logLevel > 2
                    disp('Exporting cube in FORTRAN ordering...')
                end
                permutation = [2 1 3];
            end
            
            if isfield(pln.propMC,'materialConverter')
                cubeExport = pln.propMC.materialConverter;
            else
                cubeExport = obj.materialConversion; %'RSP'; %'HUSchneiderToWater';
            end
            checkMaterial = false;
            rspCubeMethod = obj.rsp_methodCube;
            
            
            paramFile = obj.outfilenames.patientParam;
            dataFile = obj.outfilenames.patientCube;
            
            %Bookkeeping
            obj.MCparam.imageCubeOrdering = obj.arrayOrdering;
            obj.MCparam.imageCubeConversionType = obj.materialConversion;
            obj.MCparam.imageCubeFile = obj.outfilenames.patientCube;
            obj.MCparam.imageCubeDim = ct.cubeDim;
            obj.MCparam.imageVoxelDimension = ct.resolution;
            
            
            nbVoxels = prod(ct.cubeDim);
            
            outfile = fullfile(obj.workingDir, paramFile);
            obj.matRad_cfg.dispInfo('Writing data to %s\n',outfile)
            fID = fopen(outfile,'w+');
            
            
            
            switch cubeExport
                case 'CustomWaterRSP'
                    rspCube = ct.cube{1};
                    rspCube = permute(rspCube,permutation); %  X,Y,Z ordering
                    
                    fbase = fopen(['materials/' medium '.txt'],'r');
                    while ~feof(fbase)
                        strLine = fgets(fbase); %# read line by line
                        fprintf(fID,'%s',strLine);
                    end
                    fclose(fbase);
                    
                    minRSP = min(rspCube(:));
                    maxRSP = max(rspCube(:));
                    
                    % to avoid zero density
                    if minRSP<1.e-6
                        minRSP=1.e-6;
                    end
                    
                    % case of homogenous water medium (i.e., RSP=1)
                    if (length(unique(ct.cube{1})) == 1) && (unique(ct.cube{1} == 1))
                        fprintf(fID,'s:Ge/Patient/Parent="World"\n');
                        fprintf(fID,'s:Ge/Patient/Type= "TsBox"\n');
                        fprintf(fID,'s:Ge/Patient/Material = "%s"\n',medium);
                        fprintf(fID,'d:Ge/Patient/HLX      = %f mm\n',0.5*ct.cubeDim(2)*ct.resolution.x);
                        fprintf(fID,'d:Ge/Patient/HLY      = %f mm\n',0.5*ct.cubeDim(1)*ct.resolution.y);
                        fprintf(fID,'d:Ge/Patient/HLZ      = %f mm\n',0.5*ct.cubeDim(3)*ct.resolution.z);
                        fprintf(fID,'i:Ge/Patient/XBins    = %d\n',ct.cubeDim(2));
                        fprintf(fID,'i:Ge/Patient/YBins    = %d\n',ct.cubeDim(1));
                        fprintf(fID,'i:Ge/Patient/ZBins    = %d\n',ct.cubeDim(3));
                        
                        cube = NaN;
                        % otherwise
                    else
                        
                        % to avoid issues with homogenous media
                        if minRSP+1.e-6>maxRSP
                            warning('Use only one RSP value')
                            vecRSPlength = 2;
                            minRSP = 0.5*maxRSP;
                        end
                        
                        dRSP = (maxRSP-minRSP)/(vecRSPlength-1);
                        upperRSP = maxRSP+dRSP;
                        
                        ixRSP = round((rspCube-minRSP)/dRSP)+1;
                        
                        fprintf(fID,'s:Ge/Patient/Parent="World"\n');
                        fprintf(fID,'i:Ma/%s/VariableDensityBins = %d\n',medium,vecRSPlength);
                        fprintf(fID,'u:Ma/%s/VariableDensityMin = %f\n',medium,minRSP);
                        fprintf(fID,'u:Ma/%s/VariableDensityMax = %f\n',medium,upperRSP);
                        
                        if rspCubeMethod == 1
                            fprintf(fID,'s:Ge/Patient/Type= "TsBox"\n');
                            fprintf(fID,'s:Ge/Patient/Material = "%s"\n',medium);
                            fprintf(fID,'d:Ge/Patient/HLX      = %f mm\n',0.5*ct.cubeDim(2)*ct.resolution.x);
                            fprintf(fID,'d:Ge/Patient/HLY      = %f mm\n',0.5*ct.cubeDim(1)*ct.resolution.y);
                            fprintf(fID,'d:Ge/Patient/HLZ      = %f mm\n',0.5*ct.cubeDim(3)*ct.resolution.z);
                            fprintf(fID,'i:Ge/Patient/XBins    = %d\n',ct.cubeDim(2));
                            fprintf(fID,'i:Ge/Patient/YBins    = %d\n',ct.cubeDim(1));
                            fprintf(fID,'i:Ge/Patient/ZBins    = %d\n',ct.cubeDim(3));
                            fprintf(fID,'sv:Ge/Patient/VoxelMaterials = %d\n',nbVoxels);
                            
                            voxelString = num2str(ixRSP(:)'-1,['"' medium '_VariableDensityBin_%d"\n']);
                            
                            %for ix=1:nbVoxels
                            %    fprintf(h,'"%s"\n',[ medium '_VariableDensityBin_' num2str(ixRSP(ix)-1)]);
                            %end
                            fprintf(fID,voxelString);
                            
                            if checkMaterial
                                for ix=1:nbVoxels
                                    rspMaterial{ix} = [ medium '_VariableDensityBin_' num2str(ixRSP(ix)-1)];
                                end
                                materialsUsed = unique(rspMaterial);
                                
                                %      fprintf(h,'sv:Sc/ExtractData/Material = %d ',length(materialsUsed))
                                %      for ix=1:length(materialsUsed)
                                %      fprintf(h,'"%s" ',materialsUsed{ix});
                                %      end
                                %      fprintf(h,'\n')
                            end
                            fclose(fID);
                            cube = rspCube;
                            
                        elseif rspCubeMethod == 2
                            fprintf(fID,'s:Ge/Patient/Type = "TsImageCube"\n');
                            fprintf(fID,'b:Ge/Patient/DumpImagingValues = "True"\n');
                            fprintf(fID,'s:Ge/Patient/BaseMaterial = "%s"\n',medium);
                            fprintf(fID,'i:Ge/Patient/MaterialIxMax = %d\n',vecRSPlength);
                            fprintf(fID,'s:Ge/Patient/InputDirectory = "./"\n');
                            fprintf(fID,'s:Ge/Patient/InputFile = "%s"\n',dataFile);
                            fprintf(fID,'s:Ge/Patient/ImagingtoMaterialConverter = "matrad"\n');
                            fprintf(fID,'i:Ge/Patient/NumberOfVoxelsX = %d\n',ct.cubeDim(2));
                            fprintf(fID,'i:Ge/Patient/NumberOfVoxelsY = %d\n',ct.cubeDim(1));
                            fprintf(fID,'iv:Ge/Patient/NumberOfVoxelsZ = 1 %d\n',ct.cubeDim(3));
                            fprintf(fID,'d:Ge/Patient/VoxelSizeX       = %.3f mm\n',ct.resolution.x);
                            fprintf(fID,'d:Ge/Patient/VoxelSizeY       = %.3f mm\n',ct.resolution.y);
                            fprintf(fID,'dv:Ge/Patient/VoxelSizeZ       = 1 %.3f mm\n',ct.resolution.z);
                            fprintf(fID,'s:Ge/Patient/DataType  = "SHORT"\n');
                            fclose(fID);
                            
                            % write data
                            fID = fopen(fullfile(obj.workingDir, dataFile),'w');
                            fwrite(fID,ixRSP(:)-1,'short');
                            fclose(fID);
                            cube = rspCube;
                            
                        end
                    end
                    
                case {'HUToWaterSchneider','HUToWaterSchneider_LungPhantom'}
                    huCube = int32(permute(ct.cubeHU{1},permutation));
                    
                    rspHlut = matRad_readHLUT('matRad_default.hlut');
                    densityCorrection = [];
                    for i = 1:size(rspHlut,1)-1
                        startVal = rspHlut(i,1);
                        endVal = rspHlut(i+1,1);
                        range = startVal:1:endVal-1;
                        densityCorrection(end+1:end+numel(range)) = matRad_interp1(rspHlut(:,1),rspHlut(:,2),range);
                    end
                    densityCorrection(end+1) = rspHlut(end,2); %add last missing value
                    
                    
                    %Write the Schneider Converter
                    fprintf(fID,'dv:Ge/Patient/DensityCorrection = %d %s %s\n',numel(densityCorrection),num2str(densityCorrection,'%f '),'g/cm3');
                    fprintf(fID,'iv:Ge/Patient/SchneiderHounsfieldUnitSections = 2 %d %d\n',rspHlut(1,1),rspHlut(end,1)+1);
                    fprintf(fID,'uv:Ge/Patient/SchneiderDensityOffset = 1 1\n');
                    fprintf(fID,'uv:Ge/Patient/SchneiderDensityFactor = 1 0\n');
                    fprintf(fID,'uv:Ge/Patient/SchneiderDensityFactorOffset = 1 %d\n\n',-rspHlut(1,1));
                    fprintf(fID,'i:Ge/Patient/MinImagingValue = %d\n',rspHlut(1,1));
                    
                    %fprintf(h,'iv:Ge/Patient/SchneiderHUToMaterialSections = 2 %d %d\n',rspHlut(1,1),rspHlut(end,1)+1);
                    %fprintf(h,'sv:Ge/Patient/SchneiderElements = 2 "Hydrogen" "Oxygen"\n');
                    %fprintf(h,'uv:Ge/Patient/SchneiderMaterialsWeight1 = 2 0.111894 0.888106\n');
                    
                    %At least include air?
                    if contains(cubeExport,'Lung')
                        fprintf(fID,'iv:Ge/Patient/SchneiderHUToMaterialSections = 4 %d %d %d %d\n',rspHlut(1,1),rspHlut(2,1),0,rspHlut(end,1)+1);
                        fprintf(fID,'sv:Ge/Patient/SchneiderElements = 5 "Hydrogen" "Oxygen" "Nitrogen" "Carbon" "Phosphorus"\n');
                        fprintf(fID,'uv:Ge/Patient/SchneiderMaterialsWeight1 = 5 0.0 0.23479269 0.76508170 0.00012561 0.0\n');
                        fprintf(fID,'uv:Ge/Patient/SchneiderMaterialsWeight2 = 5 0.111894 0.888106 0.0 0.0 0.0\n');
                        fprintf(fID,'uv:Ge/Patient/SchneiderMaterialsWeight3 = 5 0.10404040 0.75656566 0.03131313 0.10606061 0.00202020\n');
                        fprintf(fID,'dv:Ge/Patient/SchneiderMaterialMeanExcitationEnergy = 3 85.7 75.300000 78.0 eV\n');
                    else
                        fprintf(fID,'iv:Ge/Patient/SchneiderHUToMaterialSections = 3 %d %d %d\n',rspHlut(1,1),rspHlut(2,1),rspHlut(end,1)+1);
                        fprintf(fID,'sv:Ge/Patient/SchneiderElements = 4 "Hydrogen" "Oxygen" "Nitrogen" "Carbon"\n');
                        fprintf(fID,'uv:Ge/Patient/SchneiderMaterialsWeight1 = 4 0.0 0.23479269 0.76508170 0.00012561\n');
                        fprintf(fID,'uv:Ge/Patient/SchneiderMaterialsWeight2 = 4 0.111894 0.888106 0.0 0.0\n');
                        fprintf(fID,'dv:Ge/Patient/SchneiderMaterialMeanExcitationEnergy = 2 85.7 78.0 eV\n');
                    end
                    
                    %Write the Patient
                    fprintf(fID,'s:Ge/Patient/Parent="World"\n');
                    fprintf(fID,'s:Ge/Patient/Type = "TsImageCube"\n');
                    fprintf(fID,'b:Ge/Patient/DumpImagingValues = "True"\n');
                    %fprintf(h,'s:Ge/Patient/BaseMaterial = "%s"\n',medium);
                    %fprintf(h,'i:Ge/Patient/MaterialIxMax = %d\n',vecRSPlength);
                    fprintf(fID,'s:Ge/Patient/InputDirectory = "./"\n');
                    fprintf(fID,'s:Ge/Patient/InputFile = "%s"\n',dataFile);
                    fprintf(fID,'s:Ge/Patient/ImagingtoMaterialConverter = "Schneider"\n');
                    fprintf(fID,'i:Ge/Patient/NumberOfVoxelsX = %d\n',ct.cubeDim(2));
                    fprintf(fID,'i:Ge/Patient/NumberOfVoxelsY = %d\n',ct.cubeDim(1));
                    fprintf(fID,'iv:Ge/Patient/NumberOfVoxelsZ = 1 %d\n',ct.cubeDim(3));
                    fprintf(fID,'d:Ge/Patient/VoxelSizeX       = %.3f mm\n',ct.resolution.x);
                    fprintf(fID,'d:Ge/Patient/VoxelSizeY       = %.3f mm\n',ct.resolution.y);
                    fprintf(fID,'dv:Ge/Patient/VoxelSizeZ       = 1 %.3f mm\n',ct.resolution.z);
                    fprintf(fID,'s:Ge/Patient/DataType  = "SHORT"\n');
                    %fprintf(h,'includeFile = HUtoMaterialSchneiderWater.txt');
                    fclose(fID);
                    
                    % write data
                    fID = fopen(fullfile(obj.workingDir, dataFile),'w');
                    fwrite(fID,huCube,'short');
                    fclose(fID);
                    cube = huCube;
                case {'HUToWaterSchneider_mod','HUToWaterSchneider_custom'}
                    huCube = int32(permute(ct.cubeHU{1},permutation));
                    
                    rspHlut = matRad_readHLUT('matRad_default.hlut');
                    
                    %Write the Patient
                    fprintf(fID,'s:Ge/Patient/Parent="World"\n');
                    fprintf(fID,'s:Ge/Patient/Type = "TsImageCube"\n');
                    fprintf(fID,'b:Ge/Patient/DumpImagingValues = "True"\n');
                    fprintf(fID,'s:Ge/Patient/InputDirectory = "./"\n');
                    fprintf(fID,'s:Ge/Patient/InputFile = "%s"\n',dataFile);
                    fprintf(fID,'s:Ge/Patient/ImagingtoMaterialConverter = "Schneider"\n');
                    fprintf(fID,'i:Ge/Patient/NumberOfVoxelsX = %d\n',ct.cubeDim(2));
                    fprintf(fID,'i:Ge/Patient/NumberOfVoxelsY = %d\n',ct.cubeDim(1));
                    fprintf(fID,'iv:Ge/Patient/NumberOfVoxelsZ = 1 %d\n',ct.cubeDim(3));
                    fprintf(fID,'d:Ge/Patient/VoxelSizeX       = %.3f mm\n',ct.resolution.x);
                    fprintf(fID,'d:Ge/Patient/VoxelSizeY       = %.3f mm\n',ct.resolution.y);
                    fprintf(fID,'dv:Ge/Patient/VoxelSizeZ       = 1 %.3f mm\n',ct.resolution.z);
                    fprintf(fID,'s:Ge/Patient/DataType  = "SHORT"\n');
                    
                    %Write the Schneider Converter
                    %                     fprintf(fID,'i:Ge/Patient/MinImagingValue = %d\n',rspHlut(1,1));
                    %                     fprintf(fID,'i:Ge/Patient/MinImagingValue = %d\n',-1000);
                    
                    if strcmp(cubeExport,'HUToWaterSchneider_custom')
                        fname = fullfile(obj.thisFolder,filesep,obj.converterFolder,filesep,obj.infilenames.matConv_Schneider_custom);
                    else
                        fname = fullfile(obj.thisFolder,filesep,obj.converterFolder,filesep,obj.infilenames.matConv_Schneider_mod);
                    end
                    obj.matRad_cfg.dispInfo('Reading modulation Schneider Converter from %s\n',fname);
                    matConv_mod = fileread(fname);
                    fprintf(fID,'%s\n',matConv_mod);
                    
                    if ~(isfield(ct,'modulated') && ct.modulated)
                        ct.sampledDensities = [];
                    end
                    if strcmp(cubeExport,'HUToWaterSchneider_custom')
                        fprintf(fID,'\n');
                        fprintf(fID,'iv:Ge/Patient/SchneiderHounsfieldUnitSections = 9 -1000 -98 15 23 101 2001 2995 2996 %i \n',2996 + numel(ct.sampledDensities));
                        fprintf(fID,'iv:Ge/Patient/SchneiderHUToMaterialSections = 27 -1000 -950 -120 -83 -53 -23 7 18 80 120 200 300 400 500 600 700 800 900 1000 1100 1200 1300 1400 1500 2995 2996 %i \n',2996 + numel(ct.sampledDensities));
                        fprintf(fID,'dv:Ge/Patient/DensityCorrection = %i 9.35212 5.55269 4.14652 3.41395 2.9645 2.66061 2.44144 2.27588 2.1464 2.04239 1.95698 1.8856 1.82506 1.77307 1.72791 1.68835 1.6534 1.62229 1.59442 1.56932 1.54659 1.52591 1.50701 1.48968 1.47373 1.45899 1.44534 1.43265 1.42084 1.40981 1.39949 1.3898 1.38071 1.37214 1.36406 1.35643 1.34921 1.34237 1.33587 1.3297 1.32383 1.31824 1.31291 1.30781 1.30295 1.29829 1.29383 1.28956 1.28546 1.28152 1.1284 1.12519 1.12209 1.11912 1.11625 1.11348 1.11081 1.10823 1.10574 1.10333 1.101 1.09875 1.09657 1.09445 1.0924 1.09041 1.08848 1.08661 1.08479 1.08303 1.08131 1.07964 1.07801 1.07643 1.0749 1.0734 1.07194 1.07052 1.06913 1.06778 1.06646 1.06518 1.06392 1.0627 1.0615 1.06033 1.05919 1.05808 1.05698 1.05592 1.05487 1.05386 1.05286 1.05188 1.05092 1.04999 1.04907 1.04817 1.04728 1.04643 1.04558 1.04475 1.04394 1.04314 1.04235 1.04158 1.04084 1.04009 1.03936 1.03866 1.03796 1.03726 1.0366 1.03593 1.03527 1.03463 1.03401 1.03339 1.03277 1.03218 1.03159 1.03101 1.03044 1.02988 1.02933 1.02878 1.02825 1.02772 1.0272 1.0267 1.02619 1.0257 1.02522 1.02473 1.02426 1.02379 1.02334 1.02288 1.02243 1.022 1.02156 1.02113 1.02072 1.02029 1.01988 1.01948 1.01909 1.01869 1.0183 1.01793 1.01754 1.01717 1.0168 1.01644 1.01608 1.01572 1.01538 1.01503 1.01469 1.01436 1.01403 1.01369 1.01337 1.01306 1.01273 1.01242 1.01212 1.01181 1.0115 1.01121 1.01092 1.01062 1.01033 1.01006 1.00977 1.00949 1.00922 1.00895 1.00867 1.00842 1.00815 1.00789 1.00763 1.00738 1.00713 1.00687 1.00663 1.00639 1.00614 1.00591 1.00568 1.00544 1.0052 1.00498 1.00476 1.00453 1.00431 1.00409 1.00387 1.00366 1.00345 1.00323 1.00302 1.00282 1.00261 1.00241 1.00221 1.00201 1.00181 1.00162 1.00143 1.00123 1.00104 1.00086 1.00067 1.00049 1.0003 1.00012 0.99994 0.99977 0.99959 0.99941 0.99924 0.99907 0.9989 0.99873 0.99856 0.9984 0.99824 0.99807 0.99791 0.99775 0.99759 0.99744 0.99728 0.99713 0.99697 0.99682 0.99667 0.99652 0.99637 0.99623 0.99608 0.99594 0.99579 0.99565 0.99551 0.99537 0.99523 0.99509 0.99496 0.99482 0.99469 0.99455 0.99442 0.99429 0.99416 0.99403 0.9939 0.99378 0.99365 0.99352 0.9934 0.99328 0.99315 0.99303 0.99292 0.99279 0.99267 0.99256 0.99244 0.99232 0.99221 0.99209 0.99198 0.99187 0.99176 0.99164 0.99153 0.99143 0.99131 0.9912 0.9911 0.991 0.99088 0.99078 0.99068 0.99057 0.99047 0.99037 0.99027 0.99016 0.99006 0.98997 0.98987 0.98977 0.98967 0.98958 0.98948 0.98939 0.98929 0.98919 0.9891 0.98901 0.98892 0.98882 0.98873 0.98864 0.98855 0.98846 0.98838 0.98828 0.9882 0.98811 0.98803 0.98794 0.98785 0.98777 0.98768 0.9876 0.98752 0.98743 0.98735 0.98727 0.98719 0.9871 0.98703 0.98695 0.98687 0.98679 0.98671 0.98663 0.98655 0.98648 0.9864 0.98633 0.98625 0.98617 0.9861 0.98602 0.98595 0.98588 0.9858 0.98573 0.98566 0.98559 0.98552 0.98545 0.98538 0.9853 0.98524 0.98517 0.9851 0.98503 0.98496 0.98489 0.98482 0.98476 0.98469 0.98462 0.98456 0.98449 0.98443 0.98436 0.9843 0.98423 0.98417 0.98411 0.98404 0.98398 0.98392 0.98386 0.98379 0.98373 0.98367 0.98361 0.98355 0.98349 0.98343 0.98337 0.98331 0.98325 0.98319 0.98314 0.98308 0.98302 0.98296 0.98291 0.98285 0.98279 0.98274 0.98268 0.98263 0.98257 0.98251 0.98246 0.98241 0.98235 0.9823 0.98224 0.98219 0.98214 0.98208 0.98203 0.98198 0.98193 0.98188 0.98182 0.98177 0.98172 0.98167 0.98162 0.98157 0.98152 0.98147 0.98142 0.98137 0.98132 0.98127 0.98122 0.98118 0.98113 0.98108 0.98103 0.98098 0.98094 0.98089 0.98084 0.98079 0.98075 0.98071 0.98066 0.98061 0.98057 0.98052 0.98047 0.98043 0.98039 0.98034 0.9803 0.98025 0.98021 0.98016 0.98012 0.98008 0.98003 0.97999 0.97995 0.9799 0.97986 0.97982 0.97978 0.97974 0.9797 0.97966 0.97961 0.97957 0.97953 0.97949 0.97945 0.97941 0.97937 0.97933 0.97929 0.97925 0.97921 0.97917 0.97913 0.97909 0.97905 0.97902 0.97898 0.97894 0.9789 0.97886 0.97882 0.97878 0.97875 0.97871 0.97867 0.97864 0.9786 0.97856 0.97853 0.97849 0.97845 0.97841 0.97838 0.97835 0.97831 0.97827 0.97824 0.9782 0.97817 0.97813 0.9781 0.97806 0.97803 0.97799 0.97796 0.97793 0.97789 0.97786 0.97782 0.97779 0.97776 0.97772 0.97769 0.97766 0.97762 0.97759 0.97756 0.97753 0.97749 0.97746 0.97743 0.9774 0.97736 0.97733 0.9773 0.97727 0.97724 0.97721 0.97717 0.97714 0.97711 0.97708 0.97705 0.97702 0.97699 0.97696 0.97693 0.9769 0.97687 0.97684 0.97681 0.97678 0.97675 0.97672 0.97669 0.97666 0.97663 0.9766 0.97658 0.97654 0.97651 0.97649 0.97646 0.97643 0.9764 0.97637 0.97634 0.97632 0.97629 0.97626 0.97623 0.9762 0.97618 0.97615 0.97612 0.9761 0.97607 0.97604 0.97602 0.97599 0.97596 0.97593 0.97591 0.97588 0.97585 0.97583 0.9758 0.97577 0.97575 0.97573 0.9757 0.97567 0.97565 0.97562 0.97559 0.97557 0.97555 0.97552 0.9755 0.97547 0.97544 0.97542 0.9754 0.97537 0.97534 0.97532 0.9753 0.97527 0.97525 0.97522 0.9752 0.97517 0.97515 0.97513 0.9751 0.97508 0.97506 0.97503 0.97501 0.97499 0.97496 0.97494 0.97492 0.97489 0.97487 0.97485 0.97482 0.9748 0.97478 0.97475 0.97473 0.97471 0.97469 0.97466 0.97464 0.97462 0.9746 0.97458 0.97455 0.97453 0.97451 0.97449 0.97447 0.97444 0.97442 0.9744 0.97438 0.97436 0.97433 0.97432 0.97429 0.97427 0.97425 0.97423 0.97421 0.97419 0.97417 0.97415 0.97413 0.97411 0.97409 0.97407 0.97404 0.97402 0.974 0.97398 0.97396 0.97394 0.97392 0.9739 0.97388 0.97386 0.97384 0.97382 0.9738 0.97378 0.97376 0.97375 0.97373 0.97371 0.97369 0.97367 0.97365 0.97363 0.97361 0.97359 0.97357 0.97355 0.97353 0.97352 0.9735 0.97348 0.97346 0.97344 0.97342 0.97341 0.97338 0.97337 0.97335 0.97333 0.97331 0.97329 0.97328 0.97326 0.97324 0.97322 0.9732 0.97319 0.97317 0.97315 0.97313 0.97311 0.9731 0.97308 0.97306 0.97305 0.97303 0.97301 0.97299 0.97298 0.97296 0.97294 0.97292 0.97291 0.97289 0.97287 0.97286 0.97284 0.97282 0.97281 0.97279 0.97277 0.97276 0.97274 0.97272 0.97271 0.97269 0.97267 0.97266 0.97264 0.97262 0.97261 0.97259 0.97258 0.97256 0.97254 0.97253 0.97251 0.97249 0.97248 0.97246 0.97245 0.97243 0.97242 0.9724 0.97239 0.97237 0.97235 0.97234 0.97232 0.97231 0.97229 0.97228 0.97226 0.97225 0.97223 0.97222 0.9722 0.97218 0.97217 0.97216 0.97214 0.97213 0.97211 0.9721 0.97208 0.97207 0.97205 0.97203 0.97202 0.97201 0.97199 0.97198 0.97196 0.97195 0.97193 0.97192 0.97191 0.97189 0.97188 0.97186 0.97185 0.97183 0.97182 0.97181 0.97179 0.97178 0.97176 0.97175 0.97174 0.97172 0.97171 0.97169 0.97168 0.97167 0.97165 0.97164 0.97163 0.97161 0.9716 0.97158 0.97157 0.97156 0.97154 0.97153 0.97152 0.9715 0.97149 0.97148 0.97146 0.97145 0.97144 0.97142 0.97141 0.9714 0.97139 0.97137 0.97136 0.97135 0.97133 0.97132 0.97131 0.9713 0.97128 0.97127 0.97126 0.97124 0.97123 0.97122 0.97121 0.97119 0.97118 0.97117 0.97116 0.97114 0.97113 0.97112 0.97111 0.97109 0.97108 0.97107 0.97106 0.97105 0.97103 0.97102 0.97101 0.971 0.97098 0.97097 0.97096 0.97095 0.97094 0.97093 0.97091 0.9709 0.97089 0.97088 0.97086 0.97085 0.97084 0.97083 0.97082 0.97081 0.9708 0.97078 0.97077 0.97076 0.97075 0.97074 0.97073 0.97072 0.9707 0.97069 0.97068 0.97067 0.97066 0.97065 0.97064 0.97062 0.97061 0.9706 0.97059 0.97058 0.97057 0.97056 0.94514 0.94512 0.94511 0.9451 0.94509 0.94508 0.94507 0.94506 0.94505 0.94504 0.94503 0.94502 0.94501 0.945 0.94499 0.94498 0.94497 0.94496 0.94495 0.94494 0.94492 0.94492 0.9449 0.94489 0.94488 0.94487 0.94486 0.94485 0.94484 0.94483 0.94482 0.94481 0.94569 0.94582 0.94594 0.94607 0.94619 0.95156 0.95168 0.95181 0.95194 0.95206 0.95219 0.95231 0.95244 0.95256 0.95268 0.95281 0.95293 0.95306 0.95318 0.9533 0.95342 0.95355 0.95367 0.95379 0.95391 0.95403 0.95416 0.95428 0.9544 0.95452 0.95464 0.95476 0.95488 0.955 0.95512 0.96074 0.96086 0.96098 0.9611 0.96122 0.96134 0.96146 0.96158 0.9617 0.96181 0.96193 0.96205 0.96217 0.96229 0.9624 0.96252 0.96264 0.96275 0.96287 0.96298 0.9631 0.96322 0.96333 0.96345 0.96356 0.96367 0.96379 0.9639 0.96402 0.96413 0.96818 0.96829 0.96841 0.96852 0.96864 0.96874 0.96886 0.96897 0.96909 0.96919 0.96931 0.96943 0.96953 0.96964 0.96976 0.96986 0.96998 0.97009 0.9702 0.97031 0.97042 0.97053 0.97064 0.97075 0.97087 0.97098 0.9711 0.97123 0.97135 0.97146 0.97526 0.97538 0.97549 0.97561 0.97573 0.97584 0.97596 0.97608 0.97752 0.97848 0.97944 0.98701 0.98798 0.98895 0.98992 0.99089 0.99197 0.99181 0.99166 0.99151 0.99135 0.99119 0.99104 0.99088 0.99073 0.99057 0.99042 0.99027 0.99012 0.98997 0.98982 0.98967 0.98951 0.98936 0.98874 0.98811 0.98749 0.98687 0.98625 0.98563 0.98502 0.9844 0.98378 0.98317 0.98256 0.98195 0.98134 0.98073 0.98013 0.97952 0.97892 0.97832 0.97771 0.97711 0.97651 0.97592 0.97532 0.97472 0.97414 0.97355 0.97295 0.97236 0.97178 0.97119 0.9706 0.97002 0.96943 0.96885 0.96828 0.9677 0.96712 0.96654 0.96596 0.97464 0.97406 0.97348 0.97291 0.97233 0.97176 0.97119 0.97062 0.97005 0.96948 0.96891 0.96834 0.96777 0.96721 0.96665 0.96609 0.96553 0.96497 0.96441 0.96385 0.9633 1.00233 1.00225 1.00216 1.00208 1.002 1.00192 1.00184 1.00175 1.00167 1.00159 1.00151 1.00143 1.00134 1.00126 1.00118 1.0011 1.00102 1.00093 1.00085 1.00035 1.00027 1.00019 1.0001 1.00002 0.99994 0.99987 0.99979 0.9997 0.99962 0.99954 0.99946 0.99938 0.9993 0.99922 0.99914 0.99906 0.99899 0.9989 0.99882 0.99874 0.99867 0.99859 0.9985 0.99842 0.99835 0.99827 0.99819 0.99811 0.99803 0.99795 0.99788 0.9978 0.99771 0.99764 0.99756 0.99749 0.99741 0.99732 0.99725 0.99717 0.9971 0.99702 0.99694 0.99686 0.99678 0.99671 0.99663 0.99655 0.99647 0.9964 0.99632 0.99625 0.99616 0.99609 0.99601 0.99594 0.99586 0.99578 0.99571 0.99563 0.99556 0.99548 0.9954 0.99533 0.99525 0.99518 0.99511 0.99502 0.99495 0.99488 0.9948 0.99473 0.99465 0.99457 0.9945 0.99443 0.99435 0.99427 0.9942 1.00458 1.0045 1.00443 1.00435 1.00428 1.0042 1.00413 1.00406 1.00397 1.0039 1.00383 1.00376 1.00368 1.0036 1.00353 1.00346 1.00339 1.00331 1.00323 1.00316 1.00309 1.00301 1.00294 1.00286 1.00279 1.00272 1.00265 1.00258 1.0025 1.00243 1.00235 1.00228 1.00221 1.00213 1.00206 1.00199 1.00192 1.00185 1.00177 1.0017 1.00163 1.00156 1.00149 1.00141 1.00134 1.00127 1.0012 1.00113 1.00105 1.00098 1.00091 1.00084 1.00077 1.00069 1.00062 1.00055 1.00048 1.00041 1.00033 1.00027 1.0002 1.00013 1.00006 0.99998 0.99991 0.99984 0.99977 0.99971 0.99963 0.99956 0.99949 0.99942 0.99935 0.99928 0.99921 0.99914 0.99907 0.99901 0.99893 0.99886 0.99879 0.99873 0.99866 0.99858 0.99852 0.99845 0.99838 0.99831 0.99824 0.99817 0.9981 0.99804 0.99797 0.99789 0.99783 0.99776 0.99769 0.99763 0.99755 0.99749 1.00929 1.00922 1.00916 1.00908 1.00901 1.00895 1.00888 1.00881 1.00874 1.00867 1.00861 1.00854 1.00847 1.0084 1.00833 1.00827 1.0082 1.00813 1.00806 1.00799 1.00793 1.00786 1.0078 1.00772 1.00766 1.00759 1.00753 1.00746 1.00739 1.00732 1.00726 1.00719 1.00713 1.00705 1.00699 1.00692 1.00686 1.00679 1.00672 1.00666 1.00659 1.00653 1.00646 1.00639 1.00633 1.00626 1.0062 1.00613 1.00606 1.006 1.00593 1.00587 1.0058 1.00573 1.00567 1.00561 1.00554 1.00548 1.00541 1.00534 1.00528 1.00522 1.00515 1.00508 1.00502 1.00496 1.00489 1.00483 1.00476 1.0047 1.00463 1.00457 1.00451 1.00444 1.00438 1.00431 1.00425 1.00419 1.00412 1.00406 1.00399 1.00393 1.00387 1.0038 1.00374 1.00368 1.00361 1.00355 1.00348 1.00342 1.00336 1.0033 1.00324 1.00317 1.00311 1.00305 1.00298 1.00292 1.00285 1.00279 1.01335 1.01329 1.01323 1.01316 1.0131 1.01304 1.01298 1.01292 1.01285 1.01278 1.01272 1.01266 1.0126 1.01253 1.01247 1.01241 1.01235 1.01229 1.01222 1.01216 1.0121 1.01204 1.01198 1.01191 1.01185 1.01179 1.01173 1.01167 1.0116 1.01154 1.01148 1.01142 1.01136 1.0113 1.01124 1.01118 1.01112 1.01106 1.01099 1.01093 1.01087 1.01081 1.01075 1.01069 1.01063 1.01057 1.01051 1.01045 1.01038 1.01033 1.01027 1.01021 1.01015 1.01008 1.01003 1.00997 1.00991 1.00985 1.00978 1.00972 1.00967 1.00961 1.00955 1.00948 1.00943 1.00937 1.00931 1.00925 1.00919 1.00913 1.00907 1.00902 1.00896 1.00889 1.00883 1.00878 1.00872 1.00866 1.0086 1.00854 1.00848 1.00843 1.00837 1.0083 1.00825 1.00819 1.00813 1.00808 1.00801 1.00796 1.0079 1.00784 1.00779 1.00772 1.00767 1.00761 1.00755 1.0075 1.00743 1.00738 1.01627 1.01621 1.01615 1.01609 1.01603 1.01598 1.01592 1.01587 1.0158 1.01575 1.01569 1.01563 1.01558 1.01551 1.01546 1.0154 1.01535 1.01529 1.01523 1.01517 1.01512 1.01506 1.015 1.01494 1.01489 1.01483 1.01478 1.01472 1.01466 1.0146 1.01455 1.01449 1.01444 1.01438 1.01432 1.01427 1.01421 1.01416 1.01409 1.01404 1.01398 1.01393 1.01388 1.01381 1.01376 1.01371 1.01365 1.0136 1.01353 1.01348 1.01343 1.01337 1.01332 1.01326 1.0132 1.01315 1.0131 1.01304 1.01298 1.01293 1.01287 1.01282 1.01277 1.01271 1.01265 1.0126 1.01255 1.01249 1.01243 1.01238 1.01232 1.01227 1.01222 1.01216 1.01211 1.01205 1.012 1.01195 1.01189 1.01183 1.01178 1.01173 1.01168 1.01162 1.01156 1.01151 1.01146 1.01141 1.01135 1.01129 1.01124 1.01119 1.01114 1.01108 1.01103 1.01097 1.01092 1.01087 1.01081 1.01076 1.01979 1.01974 1.01968 1.01962 1.01957 1.01952 1.01947 1.01942 1.01936 1.0193 1.01925 1.0192 1.01915 1.01909 1.01904 1.01899 1.01894 1.01888 1.01883 1.01877 1.01872 1.01867 1.01862 1.01856 1.01851 1.01846 1.01841 1.01836 1.0183 1.01825 1.0182 1.01815 1.0181 1.01804 1.01799 1.01794 1.01789 1.01783 1.01778 1.01773 1.01768 1.01763 1.01758 1.01752 1.01747 1.01742 1.01737 1.01732 1.01726 1.01721 1.01716 1.01711 1.01706 1.017 1.01695 1.0169 1.01685 1.0168 1.01675 1.0167 1.01665 1.0166 1.01655 1.01649 1.01644 1.01639 1.01634 1.0163 1.01624 1.01619 1.01614 1.01609 1.01604 1.01599 1.01594 1.01589 1.01584 1.01579 1.01573 1.01569 1.01564 1.01559 1.01554 1.01548 1.01543 1.01539 1.01534 1.01529 1.01523 1.01518 1.01514 1.01509 1.01504 1.01499 1.01494 1.01489 1.01484 1.01479 1.01474 1.01469 1.02359 1.02354 1.02349 1.02344 1.02339 1.02334 1.02329 1.02324 1.02319 1.02314 1.02309 1.02304 1.023 1.02294 1.02289 1.02285 1.0228 1.02275 1.0227 1.02265 1.0226 1.02255 1.02251 1.02245 1.0224 1.02236 1.02231 1.02226 1.02221 1.02216 1.02211 1.02207 1.02202 1.02196 1.02192 1.02187 1.02182 1.02178 1.02172 1.02168 1.02163 1.02158 1.02154 1.02148 1.02144 1.02139 1.02134 1.0213 1.02124 1.0212 1.02115 1.0211 1.02106 1.021 1.02096 1.02091 1.02087 1.02082 1.02077 1.02072 1.02067 1.02063 1.02058 1.02053 1.02049 1.02043 1.02039 1.02035 1.02029 1.02025 1.0202 1.02016 1.02011 1.02006 1.02002 1.01996 1.01992 1.01988 1.01983 1.01978 1.01973 1.01969 1.01964 1.01959 1.01955 1.0195 1.01945 1.01941 1.01936 1.01932 1.01927 1.01922 1.01918 1.01913 1.01909 1.01903 1.01899 1.01895 1.0189 1.01886 1.02794 1.02789 1.02785 1.0278 1.02776 1.02771 1.02767 1.02762 1.02756 1.02752 1.02748 1.02744 1.02739 1.02734 1.02729 1.02725 1.0272 1.02717 1.02711 1.02706 1.02702 1.02697 1.02693 1.02688 1.02684 1.02679 1.02675 1.0267 1.02666 1.02661 1.02657 1.02652 1.02648 1.02643 1.02638 1.02635 1.0263 1.02626 1.0262 1.02616 1.02611 1.02607 1.02603 1.02598 1.02594 1.02589 1.02585 1.0258 1.02576 1.02571 1.02567 1.02563 1.02558 1.02554 1.02549 1.02545 1.02541 1.02536 1.02531 1.02526 1.02523 1.02519 1.02514 1.02509 1.02505 1.025 1.02496 1.02493 1.02487 1.02483 1.02479 1.02474 1.0247 1.02465 1.02461 1.02457 1.02453 1.02448 1.02444 1.0244 1.02435 1.02431 1.02427 1.02421 1.02417 1.02414 1.02409 1.02405 1.024 1.02396 1.02391 1.02387 1.02384 1.02379 1.02374 1.0237 1.02366 1.02361 1.02357 1.02353 1.03085 1.03081 1.03077 1.03073 1.03068 1.03064 1.0306 1.03055 1.03051 1.03046 1.03042 1.03038 1.03034 1.03029 1.03025 1.03021 1.03017 1.03012 1.03008 1.03004 1.03 1.02995 1.02992 1.02986 1.02983 1.02978 1.02975 1.0297 1.02966 1.02961 1.02957 1.02954 1.02949 1.02945 1.0294 1.02937 1.02932 1.02928 1.02923 1.02919 1.02916 1.02911 1.02908 1.02902 1.02899 1.02894 1.0289 1.02887 1.02882 1.02878 1.02873 1.0287 1.02866 1.02861 1.02857 1.02853 1.02849 1.02845 1.0284 1.02836 1.02833 1.02828 1.02825 1.02819 1.02816 1.02812 1.02807 1.02804 1.02799 1.02795 1.02791 1.02788 1.02783 1.02779 1.02775 1.0277 1.02767 1.02763 1.02758 1.02754 1.02751 1.02746 1.02742 1.02738 1.02735 1.0273 1.02726 1.02723 1.02717 1.02714 1.0271 1.02707 1.02702 1.02698 1.02694 1.0269 1.02686 1.02682 1.02678 1.02673 1.03276 1.03272 1.03268 1.03264 1.03259 1.03256 1.03252 1.03248 1.03243 1.0324 1.03235 1.03232 1.03228 1.03224 1.03219 1.03216 1.03212 1.03208 1.03204 1.032 1.03196 1.03192 1.03188 1.03183 1.0318 1.03176 1.03172 1.03169 1.03164 1.0316 1.03156 1.03153 1.03149 1.03144 1.03141 1.03136 1.03133 1.03129 1.03124 1.03121 1.03117 1.03113 1.0311 1.03105 1.03101 1.03097 1.03094 1.0309 1.03085 1.03082 1.03078 1.03074 1.03071 1.03067 1.03062 1.03059 1.03055 1.03051 1.03047 1.03043 1.03039 1.03036 1.03032 1.03027 1.03024 1.0302 1.03016 1.03013 1.03009 1.03004 1.03001 1.02998 1.02993 1.02989 1.02986 1.02981 1.02978 1.02975 1.02971 1.02966 1.02963 1.02959 1.02955 1.02951 1.02948 1.02944 1.0294 1.02936 1.02932 1.02929 1.02925 1.02921 1.02918 1.02914 1.02909 1.02906 1.02903 1.02899 1.02894 1.02891 1.03634 1.03631 1.03626 1.03623 1.03619 1.03615 1.03611 1.03608 1.03603 1.036 1.03596 1.03593 1.03589 1.03585 1.03581 1.03578 1.03574 1.03571 1.03566 1.03563 1.03559 1.03556 1.03552 1.03548 1.03544 1.03541 1.03537 1.03534 1.03529 1.03526 1.03522 1.03518 1.03515 1.03511 1.03507 1.03503 1.035 1.03496 1.03492 1.03488 1.03485 1.03482 1.03478 1.03474 1.0347 1.03467 1.03463 1.0346 1.03456 1.03452 1.03449 1.03445 1.03442 1.03438 1.03434 1.03431 1.03427 1.03423 1.0342 1.03416 1.03412 1.03408 1.03405 1.03402 1.03397 1.03394 1.03391 1.03387 1.03383 1.0338 1.03376 1.03373 1.03369 1.03365 1.03362 1.03358 1.03355 1.03352 1.03347 1.03344 1.03341 1.03337 1.03334 1.0333 1.03326 1.03323 1.03319 1.03316 1.03312 1.03308 1.03305 1.03301 1.03298 1.03294 1.0329 1.03287 1.03284 1.0328 1.03276 1.03273 1.03851 1.03848 1.03844 1.03841 1.03837 1.03833 1.0383 1.03826 1.03823 1.03819 1.03816 1.03812 1.03809 1.03805 1.03802 1.03798 1.03795 1.03792 1.03788 1.03785 1.03781 1.03778 1.03774 1.03771 1.03767 1.03764 1.0376 1.03757 1.03752 1.0375 1.03746 1.03743 1.0374 1.03735 1.03732 1.03729 1.03726 1.03722 1.03718 1.03715 1.03711 1.03709 1.03705 1.03701 1.03697 1.03694 1.03692 1.03688 1.03684 1.0368 1.03677 1.03674 1.03671 1.03667 1.03663 1.0366 1.03657 1.03653 1.0365 1.03646 1.03643 1.0364 1.03636 1.03633 1.0363 1.03626 1.03623 1.03619 1.03616 1.03613 1.03609 1.03606 1.03603 1.03599 1.03596 1.03593 1.03589 1.03586 1.03583 1.03579 1.03576 1.03573 1.03569 1.03566 1.03563 1.03559 1.03556 1.03553 1.03548 1.03546 1.03543 1.03539 1.03536 1.03532 1.03528 1.03526 1.03523 1.03519 1.03515 1.03512 1.04093 1.0409 1.04087 1.04083 1.0408 1.04076 1.04074 1.0407 1.04067 1.04063 1.0406 1.04057 1.04054 1.0405 1.04047 1.04043 1.04041 1.04037 1.04034 1.0403 1.04027 1.04024 1.04021 1.04017 1.04014 1.04011 1.04008 1.04004 1.04001 1.03998 1.03994 1.03992 1.03988 1.03985 1.03981 1.03978 1.03975 1.03972 1.03968 1.03965 1.03962 1.03959 1.03956 1.03952 1.03949 1.03945 1.03943 1.0394 1.03936 1.03933 1.03929 1.03927 1.03923 1.0392 1.03917 1.03913 1.03911 1.03907 1.03904 1.03901 1.03898 1.03895 1.03891 1.03888 1.03885 1.03882 1.03879 1.03875 1.03872 1.03869 1.03866 1.03863 1.03859 1.03856 1.03853 1.0385 1.03847 1.03844 1.0384 1.03837 1.03835 1.03831 1.03828 1.03825 1.03821 1.03819 1.03815 1.03812 1.03809 1.03805 1.03803 1.038 1.03796 1.03793 1.0379 1.03787 1.03784 1.03781 1.03777 1.03774 1.04358 1.04356 1.04352 1.04349 1.04346 1.04343 1.0434 1.04337 1.04333 1.0433 1.04327 1.04324 1.04321 1.04317 1.04315 1.04311 1.04309 1.04305 1.04302 1.04299 1.04296 1.04293 1.04291 1.04286 1.04284 1.04281 1.04278 1.04275 1.04272 1.04268 1.04266 1.04262 1.0426 1.04256 1.04253 1.0425 1.04247 1.04244 1.04241 1.04237 1.04235 1.04231 1.04229 1.04225 1.04222 1.04219 1.04217 1.04213 1.0421 1.04208 1.04204 1.04202 1.04198 1.04195 1.04192 1.04189 1.04186 1.04183 1.0418 1.04177 1.04174 1.04171 1.04168 1.04164 1.04162 1.04159 1.04156 1.04154 1.04149 1.04147 1.04144 1.04141 1.04138 1.04135 1.04132 1.04129 1.04126 1.04123 1.0412 1.04117 1.04114 1.04111 1.04108 1.04105 1.04102 1.04099 1.04096 1.04093 1.0409 1.04087 1.04084 1.04081 1.04079 1.04075 1.04072 1.04069 1.04067 1.04064 1.0406 1.04058 1.04455 1.04453 1.0445 1.04447 1.04443 1.04441 1.04438 1.04435 1.04431 1.04429 1.04426 1.04423 1.0442 1.04417 1.04414 1.04411 1.04408 1.04406 1.04402 1.044 1.04397 1.04394 1.04391 1.04388 1.04385 1.04382 1.04379 1.04377 1.04373 1.0437 1.04367 1.04365 1.04362 1.04359 1.04356 1.04353 1.0435 1.04347 1.04344 1.04341 1.04339 1.04336 1.04333 1.04329 1.04327 1.04324 1.04321 1.04319 1.04315 1.04313 1.0431 1.04307 1.04304 1.04301 1.04299 1.04295 1.04293 1.0429 1.04287 1.04284 1.04281 1.04279 1.04276 1.04273 1.0427 1.04267 1.04265 1.04261 1.04259 1.04255 1.04253 1.04251 1.04247 1.04244 1.04241 1.04239 1.04236 1.04233 1.0423 1.04228 1.04225 1.04222 1.04219 1.04216 1.04214 1.0421 1.04208 1.04205 1.04202 1.042 1.04196 1.04194 1.04192 1.04188 1.04185 1.04183 1.0418 1.04178 1.04175 1.04171 1.04169 1.04166 1.04163 1.0416 1.04158 1.04155 1.04153 1.04149 1.04146 1.04144 1.04141 1.04138 1.04136 1.04133 1.0413 1.04128 1.04125 1.04122 1.04119 1.04117 1.04113 1.04111 1.04109 1.04105 1.04103 1.041 1.04097 1.04095 1.04092 1.04089 1.04086 1.04084 1.04081 1.04077 1.04075 1.04073 1.0407 1.04067 1.04064 1.04062 1.04059 1.04056 1.04054 1.04051 1.04048 1.04046 1.04043 1.0404 1.04037 1.04035 1.04032 1.04029 1.04027 1.04024 1.04021 1.04018 1.04016 1.04014 1.04011 1.04007 1.04005 1.04003 1.04 1.03997 1.03994 1.03992 1.03989 1.03986 1.03983 1.03981 1.03979 1.03975 1.03973 1.0397 1.03968 1.03965 1.03962 1.0396 1.03957 1.03955 1.03951 1.03949 1.03947 1.03944 1.03941 1.03938 1.03936 1.03934 1.03931 1.03928 1.03925 1.03923 1.03921 1.03918 1.03915 1.03912 1.0391 1.03907 1.03904 1.03902 1.03899 1.03897 1.03894 1.03891 1.03889 1.03886 1.03883 1.03881 1.03878 1.03876 1.03873 1.0387 1.03868 1.03865 1.03863 1.0386 1.03858 1.03855 1.03852 1.0385 1.03847 1.03845 1.03842 1.03839 1.03836 1.03834 1.03832 1.03829 1.03826 1.03824 1.03821 1.03819 1.03816 1.03813 1.03811 1.03809 1.03806 1.03803 1.038 1.03798 1.03796 1.03793 1.03791 1.03788 1.03786 1.03783 1.0378 1.03778 1.03775 1.03773 1.0377 1.03768 1.03766 1.03763 1.0376 1.03757 1.03755 1.03753 1.0375 1.03747 1.03745 1.03743 1.0374 1.03737 1.03735 1.03732 1.03729 1.03727 1.03724 1.03722 1.0372 1.03717 1.03715 1.03712 1.0371 1.03707 1.03704 1.03702 1.03699 1.03697 1.03694 1.03692 1.0369 1.03687 1.03684 1.03682 1.0368 1.03678 1.03675 1.03672 1.0367 1.03667 1.03665 1.03662 1.03659 1.03657 1.03655 1.03652 1.03649 1.03647 1.03645 1.03642 1.0364 1.03637 1.03635 1.03633 1.0363 1.03628 1.03625 1.03623 1.0362 1.03618 1.03616 1.03613 1.03611 1.03608 1.03606 1.03603 1.03601 1.03598 1.03596 1.03593 1.03591 1.03589 1.03586 1.03584 1.03581 1.03579 1.03576 1.03574 1.03571 1.03569 1.03567 1.03564 1.03562 1.03559 1.03557 1.03555 1.03552 1.0355 1.03547 1.03545 1.03543 1.0354 1.03538 1.03535 1.03533 1.03531 1.03528 1.03525 1.03523 1.03521 1.03519 1.03516 1.03513 1.03511 1.03509 1.03506 1.03503 1.03501 1.03499 1.03497 1.03494 1.03491 1.03489 1.03487 1.03485 1.03482 1.0348 1.03478 1.03475 1.03473 1.03471 1.03468 1.03466 1.03463 1.03461 1.03459 1.03456 1.03453 1.03451 1.03449 1.03447 1.03445 1.03442 1.0344 1.03438 1.03435 1.03432 1.0343 1.03428 1.03426 1.03423 1.0342 1.03418 1.03416 1.03414 1.03412 1.03409 1.03407 1.03405 1.03402 1.034 1.03397 1.03395 1.03393 1.03391 1.03389 1.03386 1.03384 1.03381 1.03379 1.03377 1.03374 1.03372 1.0337 1.03368 1.03365 1.03363 1.0336 1.03358 1.03356 1.03353 1.03351 1.03349 1.03347 1.03344 1.03342 1.03339 1.03337 1.03335 1.03333 1.03331 1.03328 1.03326 1.03323 1.03321 1.03319 1.03317 1.03315 1.03312 1.0331 1.03308 1.03305 1.03303 1.03301 1.03299 1.03297 1.03294 1.03291 1.03289 1.03287 1.03285 1.03282 1.0328 1.03278 1.03276 1.03273 1.03271 1.03269 1.03267 1.03264 1.03262 1.0326 1.03258 1.03255 1.03253 1.03251 1.03248 1.03246 1.03244 1.03242 1.0324 1.03237 1.03235 1.03233 1.03231 1.03229 1.03226 1.03224 1.03222 1.0322 1.03217 1.03214 1.03212 1.03211 1.03209 1.03206 1.03203 1.03201 1.032 1.03197 1.03195 1.03192 1.0319 1.03189 1.03186 1.03184 1.03181 1.03179 1.03177 1.03181 1.03185 1.03189 1.03193 1.03197 1.03201 1.03204 1.03209 1.03212 1.03216 1.0322 1.03224 1.03227 1.03232 1.03236 1.03239 1.03244 1.03247 1.03251 1.03255 1.03259 1.03262 1.03267 1.0327 1.03274 1.03278 1.03282 1.03286 1.0329 1.03294 1.03297 1.03301 1.03305 1.03309 1.03313 1.03317 1.0332 1.03324 1.03328 1.03332 1.03336 1.0334 1.03344 1.03347 1.03351 1.03355 1.03359 1.03362 1.03366 1.0337 1.03374 1.03378 1.03382 1.03386 1.03389 1.03393 1.03397 1.03401 1.03404 1.03408 1.03412 1.03416 1.03419 1.03424 1.03428 1.03431 1.03435 1.03439 1.03443 1.03446 1.0345 1.03454 1.03458 1.03461 1.03465 1.03469 1.03473 1.03477 1.0348 1.03484 1.03488 1.03492 1.03495 1.03499 1.03503 1.03507 1.0351 1.03514 1.03518 1.03522 1.03526 1.03529 1.03533 1.03536 1.03541 1.03544 1.03548 1.03551 1.03555 1.03559 1.03563 1.03567 1.0357 1.03574 1.03578 1.03582 1.03585 1.03589 1.03592 1.03596 1.036 1.03604 1.03608 1.03611 1.03615 1.03618 1.03623 1.03626 1.0363 1.03633 1.03637 1.03641 1.03645 1.03648 1.03652 1.03656 1.03659 1.03663 1.03667 1.03671 1.03674 1.03678 1.03681 1.03685 1.03688 1.03692 1.03696 1.037 1.03704 1.03707 1.03711 1.03714 1.03718 1.03722 1.03726 1.03729 1.03733 1.03736 1.0374 1.03744 1.03747 1.03751 1.03755 1.03759 1.03762 1.03766 1.03769 1.03773 1.03776 1.0378 1.03783 1.03787 1.03791 1.03795 1.03799 1.03802 1.03806 1.03809 1.03813 1.03816 1.0382 1.03823 1.03827 1.03831 1.03835 1.03838 1.03842 1.03846 1.03849 1.03853 1.03856 1.0386 1.03863 1.03867 1.0387 1.03874 1.03878 1.03881 1.03885 1.03888 1.03892 1.03896 1.03899 1.03903 1.03907 1.0391 1.03914 1.03918 1.03921 1.03925 1.03928 1.03932 1.03935 1.03939 1.03942 1.03946 1.03949 1.03953 1.03956 1.0396 1.03964 1.03967 1.03971 1.03974 1.03978 1.03981 1.03985 1.03988 1.03992 1.03995 1.03999 1.04002 1.04006 1.0401 1.04013 1.04017 1.0402 1.04024 1.04027 1.04031 1.04034 1.04038 1.04041 1.04045 1.04049 1.04052 1.04056 1.04059 1.04063 1.04066 1.0407 1.04073 1.04077 1.0408 1.04084 1.04087 1.04091 1.04094 1.04098 1.04101 1.04105 1.04108 1.04111 1.04115 1.04118 1.04122 1.04125 1.04129 1.04133 1.04136 1.0414 1.04143 1.04147 1.0415 1.04154 1.04157 1.04161 1.04164 1.04167 1.0417 1.04174 1.04178 1.04181 1.04185 1.04188 1.04192 1.04195 1.04199 1.04202 1.04205 1.04209 1.04212 1.04215 1.04219 1.04223 1.04226 1.0423 1.04233 1.04237 1.0424 1.04243 1.04246 1.0425 1.04253 1.04257 1.04261 1.04264 1.04268 1.04271 1.04274 1.04277 1.04281 1.04284 1.04288 1.04291 1.04295 1.04298 1.04301 1.04305 1.04308 1.04312 1.04315 1.04319 1.04322 1.04325 1.04329 1.04332 1.04335 1.04339 1.04343 1.04346 1.04349 1.04352 1.04356 1.04359 1.04363 1.04366 1.0437 1.04373 1.04376 1.04379 1.04383 1.04387 1.0439 1.04393 1.04396 1.044 1.04403 1.04407 1.0441 1.04413 1.04416 1.0442 1.04423 1.04427 1.0443 1.04433 1.04437 1.0444 1.04444 1.04447 1.0445 1.04453 1.04457 1.0446 1.04464 1.04467 1.0447 1.04474 1.04477 1.04481 1.04484 1.04487 1.0449 1.04494 1.04497 1.045 1.04503 1.04507 1.04511 1.04514 1.04517 1.0452 1.04524 1.04527 1.0453 1.04533 1.04537 1.0454 1.04543 1.04547 1.0455 1.04554 1.04557 1.0456 1.04563 1.04567 1.0457 1.04573 1.04576 1.0458 1.04583 1.04586 1.0459 1.04593 1.04596 1.04599 1.04603 1.04606 1.0461 1.04612 1.04616 1.04619 1.04623 1.04625 1.04629 1.04633 1.04636 1.04639 1.04642 1.04646 1.04648 1.04652 1.04655 1.04659 1.04661 1.04665 1.04669 1.04671 1.04675 1.04678 1.04681 1.04684 1.04688 1.04691 1.04694 1.04697 1.04701 1.04704 1.04707 1.04711 1.04714 1.04717 1.0472 1.04724 1.04726 1.0473 1.04733 1.04736 1.04739 1.04743 1.04746 1.04749 1.04753 1.04756 1.04759 1.04762 1.04766 1.04768 1.04772 1.04775 1.04778 1.04781 1.04785 1.04788 1.04791 1.04795 1.04797 1.04801 1.04804 1.04807 1.0481 1.04814 1.04816 1.0482 1.04823 1.04826 1.0483 1.04833 1.04836 1.04839 1.04842 1.04845 1.04849 1.04852 1.04855 1.04858 1.04861 1.04865 1.04868 1.04871 1.04874 1.04877 1.0488 1.04884 1.04886 1.0489 1.04893 1.04896 1.04899 1.04903 1.04906 1.04909 1.04912 1.04915 1.04918 1.04921 1.04925 1.04928 1.04931 1.04934 1.04937 1.04941 1.04943 1.04947 1.0495 1.04953 1.04956 1.04959 1.04962 1.04966 1.04968 1.04972 1.04975 1.04978 1.04981 1.04984 1.04988 1.0499 1.04994 1.04997 1.05 1.05003 1.05006 1.05009 1.05012 1.05015 1.05019 1.05022 1.05025 1.05028 1.05031 1.05034 1.05037 1.05041 1.05043 1.05047 1.05049 1.05053 1.05056 1.05059 1.05062 1.05065 1.05068 1.05071 1.05075 1.05077 1.05081 1.05083 1.05087 1.0509 1.05093 1.05096 1.05099 1.05102 1.05105 1.05109 1.05111 1.05115 1.05117 1.05121 1.05123 1.05127 1.0513 1.05133 1.05136 1.05139 1.05142 1.05145 1.05148 1.05151 1.05155 1.05157 1.05161 1.05163 1.05167 1.0517 1.05173 1.05176 1.05179 1.05182 1.05185 1.05188 1.05191 1.05194 1.05197 1.052 1.05203 1.05206 1.0521 1.05212 1.05216 1.05218 1.05222 1.05224 1.05228 1.0523 1.05234 1.05236 1.0524 1.05243 1.05246 1.05249 1.05252 1.05255 1.05258 1.05261 1.05264 1.05267 1.0527 1.05273 1.05276 1.05279 1.05282 1.05285 1.05288 1.05291 1.05294 1.05297 1.053 1.05303 1.05306 1.05309 1.05312 1.05315 1.05318 1.05321 1.05324 1.05327 1.0533 1.05333 1.05336 1.05339 1.05342 1.05345 1.05347 1.05351 1.05354 1.05357 1.0536 1.05363 1.05366 1.05368 1.05372 1.05374 1.05378 1.0538 1.05383 1.05386 1.05389 1.05393 1.05395 1.05398 1.05401 1.05404 1.05407 1.0541 1.05413 1.05416 1.05419 1.05422 1.05425 1.05428 1.05431 1.05434 1.05437 1.05439 1.05443 1.05445 1.05449 1.05451 1.05454 1.05457 1.0546 1.05463 1.05466 1.05469 1.05472 1.05475 1.05478 1.05481 1.05483 1.05487 1.05489 1.05492 1.05496 1.05498 1.05501 1.05504 1.05507 1.0551 1.05513 1.05516 1.05519 1.05521 1.05525 1.05527 1.0553 1.05534 1.05536 1.05539 1.05542 1.05545 1.05548 1.05551 1.05553 1.05557 1.05559 1.05562 1.05565 1.05568 1.05571 1.05574 1.05577 1.0558 1.05583 1.05585 1.05588 1.05591 1.05594 1.05597 1.056 1.05603 1.05606 1.05609 1.05611 1.05614 1.05617 1.0562 1.05623 1.05626 1.05628 1.05632 1.05634 1.05637 1.0564 1.05643 1.05646 1.05649 1.05652 1.05654 1.05657 1.0566 1.05663 1.05666 1.05669 1.05672 1.05674 1.05678 1.0568 1.05683 1.05686 1.05689 1.05691 1.05695 1.05697 1.057 1.05703 1.05706 1.05709 1.05711 1.05715 1.05717 1.0572 1.05723 1.05726 1.05728 1.05731 1.05734 1.05737 1.05739 1.05743 1.05746 1.05748 1.05751 1.05754 1.05757 1.05759 1.05763 1.05765 1.05768 1.05771 1.05774 1.05777 1.05779 1.05782 1.05785 1.05788 1.0579 1.05793 1.05796 1.05799 1.05802 1.05805 1.05807 1.0581 1.05813 1.05816 1.05819 1.05821 1.05824 1.05827 1.0583 1.05832 1.05835 1.05838 1.05841 1.05844 1.05846 1.0585 1.05852 1.05855 1.05858 1.05861 1.05863 1.05866 1.05869 1.05872 1.05874 1.05877 1.0588 1.05883 1.05886 1.05888 1.05891 1.05894 1.05897 1.05899 1.05902 1.05905 1.05908 1.0591 1.05913 1.05916 1.05919 1.05922 1.05924 1.05927 1.05929 1.05933 1.05935 1.05938 1.0594 1.05943 1.05946 1.05949 1.05952 1.05954 1.05957 1.0596 1.05963 1.05965 1.05968 1.05971 1.05974 1.05976 1.05979 1.05982 1.05984 1.05988 1.0599 1.05993 1.05995 1.05998 1.06001 1.06004 1.06006 1.06009 1.06012 1.06015 1.06018 1.0602 1.06023 1.06025 1.06028 1.06031 1.06034 1.06036 1.06039 1.06041 1.06044 1.06047 1.0605 1.06053 1.06055 1.06058 1.06061 1.06064 1.06066 1.06069 1.06071 1.06074 1.06077 1.0608 1.06083 1.06085 1.06088 1.0609 1.06093 1.06096 1.06099 1.06101 1.06104 1.06106 1.06109 1.06112 1.06115 1.06118 1.0612 1.06123 1.06125 1.06128 1.06131 1.06133 1.06136 1.06139 1.06141 1.06144 1.06147 1.06149 1.06152 1.06155 1.06158 1.0616 1.06163 1.06165 1.06168 1.06171 1.06173 1.06176 1.06179 1.06182 1.06184 1.06187 1.06189 1.06192 1.06195 1.06197 1.062 1.06203 1.06205 1.06208 1.06211 1.06213 1.06216 1.06218 1.06221 1.06224 1.06227 1.06229 1.06232 1.06234 1.06237 1.06239 1.06242 1.06245 1.06248 1.0625 1.06253 1.06256 1.06258 1.06261 1.06263 1.06266 1.06268 1.06271 1.06274 1.06277 1.06279 1.06282 1.06285 1.06287 1.0629 1.06292 1.06295 1.06297 1.063 1.06302 1.06305 1.06308 1.06311 1.06313 1.06316 1.06319 1.06321 1.06324 1.06326 1.06329 1.06331 1.06334 1.06315 1.06295 1.06275 1.06255 1.00275 ',3996 + numel(ct.sampledDensities));
                        fprintf(fID,[repmat('%g ', 1, numel(ct.sampledDensities)),'g/cm3\n'],ct.sampledDensities');
                    end
                    
                    fclose(fID);
                    
                    % write data
                    fID = fopen(fullfile(obj.workingDir, dataFile),'w');
                    fwrite(fID,huCube,'short');
                    fclose(fID);
                    cube = huCube;
                    
                otherwise
                    obj.matRad_cfg.dispError('Material Conversion rule "%s" not implemented (yet)!\n',cubeExport);
            end
            
            obj.MCparam.imageCube = cube;
            
            
        end
        
        function writeRangeShifter(obj,fID,rangeShifter,sourceToNozzleDistance)
            
            %Hardcoded PMMA range shifter for now
            pmma_rsp = 1.165;
            rsWidth = rangeShifter.eqThickness / pmma_rsp;
            
            fprintf(fID,'s:Ge/%s/Parent   = "Nozzle"\n',rangeShifter.topasID);
            fprintf(fID,'s:Ge/%s/Type     = "TsBox"\n',rangeShifter.topasID);
            fprintf(fID,'s:Ge/%s/Material = "Lucite"\n',rangeShifter.topasID);
            fprintf(fID,'d:Ge/%s/HLX      = 250 mm\n',rangeShifter.topasID);
            fprintf(fID,'d:Ge/%s/HLY      = 250 mm\n',rangeShifter.topasID);
            fprintf(fID,'d:Ge/%s/HLZ      = %f  mm\n',rangeShifter.topasID,rsWidth/2);
            fprintf(fID,'d:Ge/%s/TransX   = 500 mm * Tf/Beam/%sOut/Value\n',rangeShifter.topasID,rangeShifter.topasID);
            fprintf(fID,'d:Ge/%s/TransY   = 0   mm\n',rangeShifter.topasID);
            fprintf(fID,'d:Ge/%s/TransZ   = %f mm\n',rangeShifter.topasID,rangeShifter.sourceRashiDistance - sourceToNozzleDistance);
            
        end
        
        function writeMCparam(obj)
            %write MCparam file with basic parameters
            
            MCparam = obj.MCparam;
            save(fullfile(obj.workingDir,'MCparam.mat'),'MCparam','-v7');
        end
        
    end
end

