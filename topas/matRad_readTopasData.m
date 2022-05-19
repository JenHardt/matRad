function [topasCube,std] = matRad_readTopasData(folder,dij)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
if exist('dij')
    calcDoseDirect = false;
else
    calcDoseDirect = true;
end

load([folder filesep 'MCparam.mat']);

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
        topasTally = zeros(cubeDim(1),cubeDim(2),cubeDim(3));
        SumVarOverFields = zeros(cubeDim(1),cubeDim(2),cubeDim(3));
        
        for f = 1:MCparam.nbFields
            topasMeanDiff = zeros(cubeDim(1),cubeDim(2),cubeDim(3));
            for k = 1:MCparam.nbRuns
                genFileName = sprintf('score_%s_field%d_run%d_%s',MCparam.simLabel,f,k,tname);
                switch MCparam.outputType
                    case 'csv'
                        genFullFile = fullfile(folder,[genFileName '.csv']);
                        data{k} = matRad_readCsvData(genFullFile,cubeDim);
                    case 'binary'
                        genFullFile = fullfile(folder,[genFileName '.bin']);
                        dataRead = matRad_readBinData(genFullFile,cubeDim);
                        
                        for i = 1:length(MCparam.scoreReportQuantity)
                            data.(MCparam.scoreReportQuantity{i}){k} = dataRead{i};
                        end
                        %                         end
                    otherwise
                        error('Not implemented!');
                end
%                 topasSum = topasSum + data{k};
            end
            
            for i = 1:length(MCparam.scoreReportQuantity)
                topasSum.(MCparam.scoreReportQuantity{i}) = sum(cat(4,data.(MCparam.scoreReportQuantity{i}){:}),4);
            end

            if contains(tname,'dose','IgnoreCase',true)
                for i = 1:length(MCparam.scoreReportQuantity)
                    topasSum.(MCparam.scoreReportQuantity{i}) = correctionFactor .* topasSum.(MCparam.scoreReportQuantity{i});
                end
               
                % Calculate Standard Deviation from batches
                for k = 1:MCparam.nbRuns
                    topasMeanDiff = topasMeanDiff + (data.Sum{k} - topasSum.Sum ./ MCparam.nbRuns).^2;
                end
                % variance of the mean
                topasVarMean = topasMeanDiff./(MCparam.nbRuns - 1)./MCparam.nbRuns;
                % std of the MEAN!
                topasStdMean = sqrt(topasVarMean);
                % std of the SUM
                topasStdSum = topasStdMean * MCparam.nbRuns;
                topasVarSum = topasStdSum.^2;
                
                topasCube.([tname '_batchStd_beam' num2str(f)]) = topasStdSum;
                
                SumVarOverFields = SumVarOverFields + topasVarSum;
            elseif contains(tname,'alpha','IgnoreCase',true) || contains(tname,'beta','IgnoreCase',true) || contains(tname,'RBE','IgnoreCase',true) || contains(tname,'LET','IgnoreCase',true)               
                for i = 1:length(MCparam.scoreReportQuantity)
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
