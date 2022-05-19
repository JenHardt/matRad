function resultGUI = matRad_calcCubesMC(w,dij,scenNum)
% matRad computation of all cubes for the resultGUI struct
% which is used as result container and for visualization in matRad's GUI
%
% call
%   resultGUI = matRad_calcCubesMC(w,dij)
%   resultGUI = matRad_calcCubesMC(w,dij,scenNum)
%
% input
%   w:       bixel weight vector
%   dij:     dose influence matrix
%   scenNum: optional: number of scenario to calculated (default 1)
%
% output
%   resultGUI: matRad result struct
%
% References
%   -
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Copyright 2015 the matRad development team.
%
% This file is part of the matRad project. It is subject to the license
% terms in the LICENSE file found in the top-level directory of this
% distribution and at https://github.com/e0404/matRad/LICENSES.txt. No part
% of the matRad project, including this file, may be copied, modified,
% propagated, or distributed except according to the terms contained in the
% LICENSE file.
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin < 3
    scenNum = 1;
end

resultGUI.w = w;

% get bixel - beam correspondence
if all(w == 1)
    for i = 1:dij.numOfBeams
        beamInfo(i).suffix = ['_beam', num2str(i)];
        beamInfo(i).logIx  = ([1:dij.numOfBeams]' == i);
    end
    resultGUI.w = ones(dij.numOfBeams,1);
else
    for i = 1:dij.numOfBeams
        beamInfo(i).suffix = ['_beam', num2str(i)];
        beamInfo(i).logIx  = (dij.beamNum == i);
    end
end
beamInfo(dij.numOfBeams+1).suffix = '';
beamInfo(dij.numOfBeams+1).logIx  = true(size(resultGUI.w));
%

fnames = fieldnames(dij);

% compute physical dose for all beams individually and together
for i = 1:length(beamInfo)
    resultGUI.(['physicalDose', beamInfo(i).suffix]) = reshape(full(dij.physicalDose{scenNum} * (resultGUI.w .* beamInfo(i).logIx)),dij.doseGrid.dimensions);
end
if isfield(dij,'physicalDose_std')
    for i = 1:length(beamInfo)
        resultGUI.(['physicalDose_std', beamInfo(i).suffix]) = sqrt(reshape(full(dij.physicalDose_std{scenNum}.^2 * (resultGUI.w .* beamInfo(i).logIx)),dij.doseGrid.dimensions));
    end
end
if isfield(dij,'physicalDose_batchStd')
    for i = 1:length(beamInfo)
        resultGUI.(['physicalDose_batchStd', beamInfo(i).suffix]) = sqrt(reshape(full(dij.physicalDose_batchStd{scenNum}.^2 * (resultGUI.w .* beamInfo(i).logIx)),dij.doseGrid.dimensions));
    end
end

% compute doseToWater analogously
if isfield(dij,'doseToWater')
    for i = 1:length(beamInfo)
        resultGUI.(['doseToWater', beamInfo(i).suffix]) = reshape(full(dij.doseToWater{scenNum} * (resultGUI.w .* beamInfo(i).logIx)),dij.doseGrid.dimensions);
    end
end
if isfield(dij,'doseToWater_std')
    for i = 1:length(beamInfo)
        resultGUI.(['doseToWater_std', beamInfo(i).suffix]) = sqrt(reshape(full(dij.doseToWater_std{scenNum}.^2 * (resultGUI.w .* beamInfo(i).logIx)),dij.doseGrid.dimensions));
    end
end
if isfield(dij,'doseToWater_batchStd')
    for i = 1:length(beamInfo)
        resultGUI.(['doseToWater_batchStd', beamInfo(i).suffix]) = sqrt(reshape(full(dij.doseToWater_batchStd{scenNum}.^2 * (resultGUI.w .* beamInfo(i).logIx)),dij.doseGrid.dimensions));
    end
end

% consider LET
if isfield(dij,'LET')
    for i = 1:length(beamInfo)
        resultGUI.(['LET', beamInfo(i).suffix])     = reshape(full(dij.LET{scenNum} * beamInfo(i).logIx),dij.doseGrid.dimensions);
        % Set LET zero where the dose is less or equal to 0
        resultGUI.(['LET', beamInfo(i).suffix])(resultGUI.(['physicalDose', beamInfo(i).suffix]) <= 0) = 0;
    end
end



% consider biological optimization
if any(contains(dij.MC_tallies,'alpha'))
    for i = 1:length(beamInfo)

        wBeam = (resultGUI.w .* beamInfo(i).logIx);
        ix = dij.bx(:,scenNum)~=0 & resultGUI.(['physicalDose', beamInfo(i).suffix])(:) > 0;

        if i < 4
            abFields = fnames(contains(fnames,{'alpha','beta'}));
            for j = 1:numel(abFields)
                resultGUI.([abFields{j}, beamInfo(i).suffix])       = reshape(full(dij.(abFields{j}){scenNum} * beamInfo(i).logIx),dij.doseGrid.dimensions);
                resultGUI.([abFields{j}, beamInfo(i).suffix])(~ix)	= 0;
            end
        end

        for ab = 1:length(abFields)
            model = strsplit(abFields{ab},'_');
            models{ab} = char(model(2));
        end
        models = unique(models);

        for m = 1:length(models)
        resultGUI.(['effect_' models{m} beamInfo(i).suffix])       = full(dij.(['mAlphaDose_' models{m}]){scenNum} * wBeam + (dij.(['mSqrtBetaDose_' models{m}]){scenNum} * wBeam).^2);
        resultGUI.(['effect_' models{m} beamInfo(i).suffix])       = reshape(resultGUI.(['effect_' models{m} beamInfo(i).suffix]),dij.doseGrid.dimensions);

        resultGUI.(['RBExD_' models{m} beamInfo(i).suffix])        = zeros(size(resultGUI.(['effect_' models{m} beamInfo(i).suffix])));
        resultGUI.(['RBExD_' models{m} beamInfo(i).suffix])(ix)    = (sqrt(dij.ax(ix).^2 + 4 .* dij.bx(ix) .* resultGUI.(['effect_' models{m} beamInfo(i).suffix])(ix)) - dij.ax(ix))./(2.*dij.bx(ix));

        resultGUI.(['RBE_' models{m} beamInfo(i).suffix])          = resultGUI.(['RBExD_' models{m} beamInfo(i).suffix])./resultGUI.(['physicalDose', beamInfo(i).suffix]);        
        end
    end
end

% Add non-processed MC tallies
% Note that the tallies are already computed per beam and altogether
% if isfield(dij,'MC_tallies')
%     for f = 1:numel(dij.MC_tallies)
%         tally = dij.MC_tallies{f};
%         % skip tallies processed above
%         if ~isfield(resultGUI,tally) && ~contains(tally,'std')
%             tallyCut = strsplit(tally,'_');
%             if size(tallyCut,2) > 1
%                 beamNum = str2num(cell2mat(regexp(tally,'\d','Match')));
%                 resultGUI.(tally) = reshape(full(dij.(tallyCut{1}){scenNum}(:,beamNum)),dij.doseGrid.dimensions);
%             else
%                 resultGUI.(tally) = reshape(full(dij.(tally){scenNum} * wBeam),dij.doseGrid.dimensions);
%             end
%         end
%     end
% end

% group similar fields together
resultGUI = orderfields(resultGUI);

% interpolation if dose grid does not match ct grid
if any(dij.ctGrid.dimensions~=dij.doseGrid.dimensions)
    myFields = fieldnames(resultGUI);
    for i = 1:numel(myFields)

        if numel(resultGUI.(myFields{i})) == dij.doseGrid.numOfVoxels

            % interpolate!
            resultGUI.(myFields{i}) = matRad_interp3(dij.doseGrid.x,dij.doseGrid.y',dij.doseGrid.z, ...
                resultGUI.(myFields{i}), ...
                dij.ctGrid.x,dij.ctGrid.y',dij.ctGrid.z,'linear',0);

        end

    end
end

end



