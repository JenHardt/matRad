classdef matRad_spotRemovalDij
    %MATRAD_SPOTREMOVALDIJ2 Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        
        defaultRelativeThreshold
        defaultAbsoluteThreshold
        weights
        newSpots
        
    end
    
    methods
        function reset(obj)
            %Set all default properties for matRad's computations
            obj.setDefaultProperties();
        end

        function setDefaultProperties(obj)
            obj.defaultRelativeThreshold = 0.03;
            obj.defaultAbsoluteThreshold = 100;
            
        end
        
        function obj = matRad_spotRemovalDij(dij,w)

            obj.weights = w;
            obj.newSpots = obj.weights>thres*mean(obj.weights);
            
        end
        
        function obj = calcNewSpots(obj)
            
            
        end
        
        function [dij,stf] = parseOutput(obj,dij,stf)
            
            if ((sum(obj.newSpots) ~= numel(w)) && sum(obj.newSpots) ~= dij.totalNumOfBixels) && any(size(w)>1)
                dij.cutWeights = w(obj.newSpots);
                
                dij.bixelNum = dij.bixelNum(obj.newSpots);
                dij.rayNum = dij.rayNum(obj.newSpots);
                dij.beamNum = dij.beamNum(obj.newSpots);
                dij.totalNumOfBixels = sum(obj.newSpots);
                
                dij.physicalDose{1} = dij.physicalDose{1}(:,obj.newSpots);
                if isfield(dij,'mAlphaDose')
                    dij.mAlphaDose{1} = dij.mAlphaDose{1}(:,obj.newSpots);
                    dij.mSqrtBetaDose{1} = dij.mSqrtBetaDose{1}(:,obj.newSpots);
                end
                if isfield(dij,'mLETDose')
                    dij.mLETDose{1} = dij.mLETDose{1}(:,obj.newSpots);
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
                dij.numOfRemovedSpots = sum(~obj.newSpots);
                matRad_cfg.dispWarning([num2str(sum(~obj.newSpots)),'/',num2str(numel(obj.newSpots)) ,' spots have been removed below ',num2str(100*thres),'% of the mean weight.\n'])
            else
                matRad_cfg.dispWarning('no spots have been removed.')
                dij.cutWeights = w;
            end
            
        end
    end
end

