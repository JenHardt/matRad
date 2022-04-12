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
MCparam.tallies = unique(MCparam.tallies);

if calcDoseDirect
    for t = 1:length(MCparam.tallies)
        tname = MCparam.tallies{t};
        topasTally = zeros(cubeDim(1),cubeDim(2),cubeDim(3));
        SumVarOverFields = zeros(cubeDim(1),cubeDim(2),cubeDim(3));
        
        for f = 1:MCparam.nbFields
            topasSum = zeros(cubeDim(1),cubeDim(2),cubeDim(3));
            topasMeanDiff = zeros(cubeDim(1),cubeDim(2),cubeDim(3));
            for k = 1:MCparam.nbRuns
                genFileName = sprintf('score_%s_field%d_run%d_%s',MCparam.simLabel,f,k,tname);
                switch MCparam.outputType
                    case 'csv'
                        genFullFile = fullfile(folder,[genFileName '.csv']);
                        data{k} = matRad_readCsvData(genFullFile,cubeDim);
                    case 'binary'
                        genFullFile = fullfile(folder,[genFileName '.bin']);
                        %                         if iscell(MCparam.scoreReportQuantity)
                        %                             [data{k},topasStd] = matRad_readBinData(genFullFile,cubeDim,numel(MCparam.scoreReportQuantity));
                        %                         else
                        data{k} = matRad_readBinData(genFullFile,cubeDim);
                        %                         end
                    otherwise
                        error('Not implemented!');
                end
                topasSum = topasSum + data{k};
            end
            
            if contains(tname,'dose','IgnoreCase',true)
                topasSum = correctionFactor.*topasSum;
                
                % Calculate Standard Deviation from batches
                for k = 1:MCparam.nbRuns
                    topasMeanDiff = topasMeanDiff + (data{k} - topasSum./MCparam.nbRuns).^2;
                end
                % variance of the mean
                topasVarMean = topasMeanDiff./(MCparam.nbRuns - 1)./MCparam.nbRuns;
                % std of the MEAN!
                topasStdMean = sqrt(topasVarMean);
                % std of the SUM
                topasStdSum = topasStdMean * MCparam.nbRuns;
                topasVarSum = topasStdSum.^2;
                
                topasCube.([tname '_std_beam' num2str(f)]) = topasStdSum;
                
                SumVarOverFields = SumVarOverFields + topasVarSum;
            elseif contains(tname,'alpha','IgnoreCase',true) || contains(tname,'beta','IgnoreCase',true) || contains(tname,'RBE','IgnoreCase',true) || contains(tname,'LET','IgnoreCase',true)
                topasSum = topasSum./ MCparam.nbRuns;
            end
            
            % Accumulate over the fields
            topasTally = topasTally + topasSum;
            % Tally per field
            topasCube.([tname '_beam' num2str(f)]) = topasSum;
            
        end
%         if contains(tname,'dose','IgnoreCase',true)
%             topasCube.([tname '_std']) = sqrt(SumVarOverFields);
%         end
%         if ~(contains(tname,'alpha','IgnoreCase',true) || contains(tname,'beta','IgnoreCase',true) || contains(tname,'RBE','IgnoreCase',true))
%             topasCube.(tname) = topasTally;
%         end
    end
    
else % if topas dij calculation
    for t = 1:length(MCparam.tallies)
        tname = MCparam.tallies{t};
        topasTally = zeros(cubeDim(1),cubeDim(2),cubeDim(3));
        
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
