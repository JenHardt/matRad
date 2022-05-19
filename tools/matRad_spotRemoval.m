function [dij,stf] = matRad_spotRemoval(dij,w,varargin)
% matRad spot removal tool
%
% call
%   dij =           matRad_spotRemoval(dij,w)
%   [dij,stf] =     matRad_spotRemoval(dij,w,stf)
%   machine = matRad_getAlphaBetaCurves(machine,cst,modelName,overrideAB)
% Example full call for protons
%   machine = matRad_getAlphaBetaCurves(machine,pln,cst,'MCN','override')
% input
%   machine:                matRad machine file to change
%   varargin (optional):    cst:        matRad cst struct (for custom alpha/beta,
%                                       otherwise default is alpha=0.1, beta=0.05;
%                           modelName:  specify RBE modelName
%                   	    overrideAB: calculate new alpha beta even if available
%                                       and override
%
% output
%   machine:                updated machine file with alpha/beta curves
%
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Copyright 2021 the matRad development team.
%
% This file is part of the matRad project. It is subject to the license
% terms in the LICENSE file found in the top-level directory of this
% distribution and at https://github.com/e0404/matRad/LICENSES.txt. No part
% of the matRad project, including this file, may be copied, modified,
% propagated, or distributed except according to the terms contained in the
% LICENSE file.
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

matRad_cfg = MatRad_Config.instance();

% set threshold for spot removal to 3% of the mean weight.
if ~isempty(varargin)
    for i = 1:nargin-2
        if isstruct(varargin{i})
            stf = varargin{i};
        elseif isscalar(varargin{i})
            thres = varargin{i};
        end
    end
end

if ~exist('thres','var')
    thres = 0.03;
%     thres = 0.00001;
end


if exist('stf','var') && nargout > 1
    calcStf = true;
end

newSpots = w>thres*mean(w);
% newSpots = w>thres*max(w);







if ((sum(newSpots) ~= numel(w)) && sum(newSpots) ~= dij.totalNumOfBixels) && any(size(w)>1)
    dij.cutWeights = w(newSpots);
    
    dij.bixelNum = dij.bixelNum(newSpots);
    dij.rayNum = dij.rayNum(newSpots);
    dij.beamNum = dij.beamNum(newSpots);
    dij.totalNumOfBixels = sum(newSpots);
    
    dij.physicalDose{1} = dij.physicalDose{1}(:,newSpots);
    if isfield(dij,'mAlphaDose')
        dij.mAlphaDose{1} = dij.mAlphaDose{1}(:,newSpots);
        dij.mSqrtBetaDose{1} = dij.mSqrtBetaDose{1}(:,newSpots);
    end
    if isfield(dij,'mLETDose')
        dij.mLETDose{1} = dij.mLETDose{1}(:,newSpots);
    end
    [~,beamNumIdx] = unique(dij.beamNum);
    beamNumIdx = [0;beamNumIdx(2:end)-1;dij.totalNumOfBixels];
    
    for b = 1:dij.numOfBeams
        currRaysInBeam = dij.rayNum(beamNumIdx(b)+1:beamNumIdx(b+1));
        currBixelsInRay = dij.bixelNum(beamNumIdx(b)+1:beamNumIdx(b+1));
        [rayCnt,rayIdx] = unique(currRaysInBeam);
        
        if calcStf
            numOfBixelsPerRay = groupcounts(currRaysInBeam);
            cutRays = ismember([1:dij.numOfRaysPerBeam(b)]',rayCnt);
            if any(~cutRays)
                stf(b).ray = stf(b).ray(cutRays);
                stf(b).numOfRays = sum(cutRays);
            end
            for i = 1:stf(b).numOfRays
                bixelCurrRay{i} = currBixelsInRay(rayIdx(i):rayIdx(i)+numOfBixelsPerRay(i)-1);
            end
            for f = 1:stf(b).numOfRays
                stf(b).ray(f).energy = stf(b).ray(f).energy(bixelCurrRay{f});
                stf(b).ray(f).focusIx = stf(b).ray(f).focusIx(bixelCurrRay{f});
                stf(b).ray(f).rangeShifter = stf(b).ray(f).rangeShifter(bixelCurrRay{f});
            end
            stf(b).numOfBixelsPerRay = numOfBixelsPerRay';
            stf(b).totalNumOfBixels = sum(stf(b).numOfBixelsPerRay);
        end
        
        dij.numOfRaysPerBeam(b) = numel(rayCnt);
    end
    
    dij.totalNumOfRays = sum(dij.numOfRaysPerBeam);
    dij.numOfRemovedSpots = sum(~newSpots);
    matRad_cfg.dispWarning([num2str(sum(~newSpots)),'/',num2str(numel(newSpots)) ,' spots have been removed below ',num2str(100*thres),'% of the mean weight.\n'])
else
    matRad_cfg.dispWarning('no spots have been removed.')
    dij.cutWeights = w;
end
end

