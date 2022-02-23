function pln = matRad_setPlanParameters(ct,cst,baseData,RBE_model)
pln.propStf.gantryAngles = [40,340,300];
pln.propStf.couchAngles = [0,0,0];
pln.propStf.bixelWidth = 3;
pln.propStf.numOfBeams = size(pln.propStf.gantryAngles,2);
pln.propStf.isoCenter = ones(pln.propStf.numOfBeams,1)*matRad_getIsoCenter(cst,ct);

pln.numOfFractions = 60;

s = split(baseData,'_');
pln.radiationMode   = s{1};
pln.machine = strjoin(s(2:end),'_');

if contains(RBE_model,{'WED','MCN','LEM','RBE'})
    quantityOpt                 = 'RBExD';
    pln.propOpt.bioOptimization = 'RBExD';
else
    quantityOpt                 = 'none';
    pln.propOpt.bioOptimization = 'none';
end

pln.propDoseCalc.lateralCutOff = 0.99;
pln.propDoseCalc.airOffsetCorrection = 1;
pln.propDoseCalc.calcLET = 0;

% retrieve scenarios for dose calculation and optimziation
% pln = rmfield(pln,{'bioParam','multScen','DicomInfo'});
pln.bioParam = matRad_bioModel(pln.radiationMode,quantityOpt,RBE_model);
pln.multScen = matRad_multScen(ct,'nomScen'); % optimize on the nominal scenario
pln.propDoseCalc.resolution = struct('x',1,'y',1,'z',1);

pln.propHeterogeneity.calcHetero = true;
pln.propOpt.runDAO        = 0;
pln.propOpt.runSequencing = 0;

% Monte Carlo settings
pln.propMC.numOfRuns = 5;
pln.propMC.engine = 'TOPAS';
pln.propMC.histories = 1e8;
pln.propMC.externalCalculation = true;
end