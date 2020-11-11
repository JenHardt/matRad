function [topasCube] = matRad_readTopasData(folder,dij)
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

if calcDoseDirect 
    for t = 1:length(MCparam.tallies)
        tname = MCparam.tallies{t};
        topasTally = zeros(cubeDim(1),cubeDim(2),cubeDim(3));
        
        for f = 1:MCparam.nbFields
            topasSum = zeros(cubeDim(1),cubeDim(2),cubeDim(3));
            for k = 1:MCparam.nbRuns
                genFileName = sprintf('score_%s_field%d_run%d_%s',MCparam.simLabel,f,k,tname);
                switch MCparam.outputType
                    case 'csv'
                        genFullFile = fullfile(folder,[genFileName '.csv']);
                        data = readCsvData(genFullFile,cubeDim);
                    case 'binary'
                        genFullFile = fullfile(folder,[genFileName '.bin']);
                        data = readBinData(genFullFile,cubeDim);
                    otherwise
                        error('Not implemented!');
                end
                topasSum = topasSum + data;
            end
            
            topasSum = correctionFactor.*topasSum;
            
            % Tally per field
            topasCube.([tname '_beam' num2str(f)]) = topasSum;
            
            % Accumulate over the fields
            topasTally = topasTally + topasSum;
        end
        
        topasCube.(tname) = topasTally;
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
                            data = readCsvData(genFullFile,cubeDim);
                        case 'binary'
                            genFullFile = fullfile(folder,[genFileName '.bin']);
                            data = readBinData(genFullFile,cubeDim);
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

function data = readCsvData(csvFile,cubeDim)
    data = zeros(cubeDim(2),cubeDim(1),cubeDim(3));
    fID = fopen(csvFile,'r');
    dataCsv = textscan(fID,'%d %d %d %f','Delimiter',',','CommentStyle','#','CollectOutput',true);
    fclose(fID);
    ix = sub2ind([cubeDim(1) cubeDim(2) cubeDim(3)],dataCsv{1}(:,2)+1,dataCsv{1}(:,1)+1,dataCsv{1}(:,3)+1);
    data(ix) = dataCsv{2};
end

function data = readBinData(binFile,cubeDim)
    fID = fopen(binFile);
    data = fread(fID,prod(cubeDim),'double');
    fclose(fID);
    data = reshape(data,cubeDim(2),cubeDim(1),cubeDim(3));
    data = permute(data,[2 1 3]);
end

