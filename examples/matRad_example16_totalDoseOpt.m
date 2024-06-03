%% Example: Photon Treatment Plan
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Copyright 2017 the matRad development team. 
% 
% This file is part of the matRad project. It is subject to the license 
% terms in the LICENSE file found in the top-level directory of this 
% distribution and at https://github.com/e0404/matRad/LICENSES.txt. No part 
% of the matRad project, including this file, may be copied, modified, 
% propagated, or distributed except according to the terms contained in the 
% LICENSE file.
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
matRad_rc; %If this throws an error, run it from the parent directory first to set the paths
load('TG119.mat');


%% Treatment Plan

pln.radiationMode = 'photons';  
pln.machine       = 'Generic';

quantityOpt    = 'physicalDose';                                     
modelName      = 'none';  
pln.numOfFractions         = 30;
pln.propStf.gantryAngles   = [0:40:359];
pln.propStf.couchAngles    = zeros(1,numel(pln.propStf.gantryAngles));
pln.propStf.bixelWidth     = 5;
pln.propStf.numOfBeams      = numel(pln.propStf.gantryAngles);
pln.propStf.isoCenter       = ones(pln.propStf.numOfBeams,1) * matRad_getIsoCenter(cst,ct,0);

pln.propDoseCalc.doseGrid.resolution.x = 3; % [mm]
pln.propDoseCalc.doseGrid.resolution.y = 3; % [mm]
pln.propDoseCalc.doseGrid.resolution.z = 3; % [mm
pln.propSeq.runSequencing = 0;

pln.propOpt.runDAO        = 0;
pln.propOpt.calcTotalDose = 1;

% retrieve bio model parameters
pln.bioParam = matRad_bioModel(pln.radiationMode,quantityOpt, modelName);

% retrieve scenarios for dose calculation and optimziation
pln.multScen = matRad_multScen(ct,'nomScen');


%% Generate Beam Geometry STF
stf = matRad_generateStf(ct,cst,pln);

%% Dose Calculation
dij = matRad_calcDoseInfluence(ct,cst,stf,pln);

%% Inverse Optimization for IMRT
resultGUI = matRad_fluenceOptimization(dij,cst,pln);
%%
resultGUItmp =  matRad_calcTotalDose(resultGUI,pln.numOfFractions);
resultGUI = matRad_appendResultGUI(resultGUI,resultGUItmp,0,'TD');

