function topasCube = matRad_readTopasData(folder,dij)
% function to read out 
% for beams in x-, y- or z-direction
%
% call
%   topasCube = matRad_readTopasData(folder,dij)
%
% input
%   folder:         Path to folder where TOPAS files are in (as string)
%   dij:            dij struct (this part needs update)
%
% output
%   topasCube:      struct with all read out subfields
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Copyright 2022 the matRad development team.
%
% This file is part of the matRad project. It is subject to the license
% terms in the LICENSE file found in the top-level directory of this
% distribution and at https://github.com/e0404/matRad/LICENSES.txt. No part
% of the matRad project, including this file, may be copied, modified,
% propagated, or distributed except according to the terms contained in the
% LICENSE file.
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if exist('dij')
    calcDoseDirect = false;
else
    calcDoseDirect = true;
end

load([folder filesep 'MCparam.mat']);
if iscell(MCparam.scoreReportQuantity)
    MCparam.numOfReportQuantities = length(MCparam.scoreReportQuantity);
else
    MCparam.numOfReportQuantities = 1;
    MCparam.scoreReportQuantity = {MCparam.scoreReportQuantity};
end


%Normalize with histories and particles/weight
correctionFactor = 1e6 / double(MCparam.nbHistoriesTotal); %double(MCparam.nbParticlesTotal) / double(MCparam.nbHistoriesTotal);

cubeDim = MCparam.imageCubeDim;
% load in all data
switch MCparam.outputType
    case 'csv'
        files = dir([folder filesep 'score_matRad_plan_field1_run1_*.csv']);
        MCparam.tallies = cellfun(@(x) erase(x,{'score_matRad_plan_field1_run1_','.csv'}) ,{files(:).name} ,'UniformOutput' ,false);
    case 'binary'
        files = dir([folder filesep 'score_matRad_plan_field1_run1_*.bin']);
        MCparam.tallies = cellfun(@(x) erase(x,{'score_matRad_plan_field1_run1_','.bin'}) ,{files(:).name} ,'UniformOutput' ,false);
end

if calcDoseDirect
    for t = 1:length(MCparam.tallies)
        tname = MCparam.tallies{t};
        
        for f = 1:MCparam.nbFields
            topasMeanDiff = zeros(cubeDim(1),cubeDim(2),cubeDim(3));
            for k = 1:MCparam.nbRuns
                genFileName = sprintf('score_%s_field%d_run%d_%s',MCparam.simLabel,f,k,tname);
                switch MCparam.outputType
                    case 'csv'
                        % output as csv needs to be cleaned up
                        genFullFile = fullfile(folder,[genFileName '.csv']);
                        data{k} = matRad_readCsvData(genFullFile,cubeDim);
                    case 'binary'
                        genFullFile = fullfile(folder,[genFileName '.bin']);
                        dataRead = matRad_readBinData(genFullFile,cubeDim);
                        
                        % binheader file could be used here to determine the actually used reportQuantities
                        for i = 1:numel(dataRead)
                            data.(MCparam.scoreReportQuantity{i}){k} = dataRead{i};
                        end
                        % for example the standard deviation is not calculated for alpha/beta so a loop through all
                        % reportQuantities does not work here
                        currNumOfQuantities = numel(dataRead);
%                        currReportQuantities = ;
                    otherwise
                        error('Not implemented!');
                end
            end

            % add STD quadratically
            for i = 1:currNumOfQuantities
                if contains(MCparam.scoreReportQuantity{i},'standard_deviation','IgnoreCase',true)
                    topasSum.(MCparam.scoreReportQuantity{i}) = sqrt(double(MCparam.nbHistoriesTotal)) * sqrt(sum(cat(4,data.(MCparam.scoreReportQuantity{i}){:}).^2,4));
                else
                    topasSum.(MCparam.scoreReportQuantity{i}) = sum(cat(4,data.(MCparam.scoreReportQuantity{i}){:}),4);
                end
            end

            if contains(tname,'dose','IgnoreCase',true)
                % Calculate Standard Deviation from batches
                for k = 1:MCparam.nbRuns
                    topasMeanDiff = topasMeanDiff + (data.Sum{k} - topasSum.Sum / MCparam.nbRuns).^2;
                end
                % variance of the mean
                topasVarMean = topasMeanDiff./(MCparam.nbRuns - 1)./MCparam.nbRuns;
                % std of the MEAN!
                topasStdMean = sqrt(topasVarMean);
                % std of the SUM
                topasStdSum = topasStdMean * correctionFactor * MCparam.nbRuns;
%                 topasVarSum = topasStdSum.^2;
                
                topasCube.([tname '_batchStd_beam' num2str(f)]) = topasStdSum;
                
%                 SumVarOverFields = SumVarOverFields + topasVarSum;                
                
                for i = 1:currNumOfQuantities
                    topasSum.(MCparam.scoreReportQuantity{i}) = correctionFactor .* topasSum.(MCparam.scoreReportQuantity{i});
                end

            elseif contains(tname,'alpha','IgnoreCase',true) || contains(tname,'beta','IgnoreCase',true) || contains(tname,'RBE','IgnoreCase',true) || contains(tname,'LET','IgnoreCase',true)               
                for i = 1:currNumOfQuantities
                    topasSum.(MCparam.scoreReportQuantity{i}) = topasSum.(MCparam.scoreReportQuantity{i}) ./ MCparam.nbRuns;
                end
            end
            
            % Tally per field
            if isfield(topasSum,'Sum')
                topasCube.([tname '_beam' num2str(f)]) = topasSum.Sum;
            end
            if isfield(topasSum,'Standard_Deviation')
                topasCube.([tname '_std_beam' num2str(f)]) = topasSum.Standard_Deviation;
            end
            
        end
    end
    
else % if topas dij calculation
    for t = 1:length(MCparam.tallies)
        tname = MCparam.tallies{t};
       
        for f = 1:MCparam.nbFields
            topasSum = cell(1,dij.totalNumOfBixels);
            topasSum(1,:) = {zeros(cubeDim(1),cubeDim(2),cubeDim(3))};
            for i = 1:dij.totalNumOfBixels
                for k = 1:MCparam.nbRuns
                    
                    genFileName = sprintf('score_%s_field%d_run%d_%s-ray%i_bixel%i',MCparam.simLabel,f,k,tname,dij.rayNum(i),dij.bixelNum(i));
                    switch MCparam.outputType
                        case 'csv'
                            genFullFile = fullfile(folder,[genFileName '.csv']);
                            data = matRad_readCsvData(genFullFile,cubeDim);
                        case 'binary'
                            genFullFile = fullfile(folder,[genFileName '.bin']);
                            data = matRad_readBinData(genFullFile,cubeDim);
                        otherwise
                            error('Not implemented!');
                    end
                    topasSum{i} = topasSum{i} + data;
                end
                
                topasSum{i} = correctionFactor.*topasSum{i};
                
                % Tally per field
                topasCube.([tname '_beam' num2str(f)]){i} = topasSum{i};
            end
        end
    end
end



end
