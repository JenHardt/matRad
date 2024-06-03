function result = matRad_calcTotalDose(resultGUI,numOfFractions)

    numBeams = fieldnames(resultGUI);
    numBeams = numBeams(cellfun(@(x) contains(x,'beam'),numBeams));
    numBeams = cellfun(@(x) regexp(x,'beam(\d+)', 'tokens'),numBeams);
    numBeams = unique(cellfun(@(x) str2double(x),numBeams));
    numBeams = max(numBeams);
    
    for i = 1:numBeams
        beamInfo(i).suffix = ['_beam', num2str(i)];    
    end
    beamInfo(numBeams+1).suffix = '';
    
    % physical Dose
    if isfield(resultGUI,'physicalDose')
        for i = 1:numBeams+1 
            result.(['physicalDose', beamInfo(i).suffix ]) = numOfFractions.*resultGUI.(['physicalDose', beamInfo(i).suffix]);
        end
    end
    % effect
    if isfield(resultGUI, 'alphaDoseCube')
        for i = 1:numBeams+1
            result.(['alphaDoseCube', beamInfo(i).suffix ]) = numOfFractions.*resultGUI.(['alphaDoseCube', beamInfo(i).suffix]);
            result.(['SqrtBetaDoseCube', beamInfo(i).suffix ]) = numOfFractions.*resultGUI.(['SqrtBetaDoseCube', beamInfo(i).suffix]);
            result.(['effect', beamInfo(i).suffix ]) = numOfFractions.*resultGUI.(['effect', beamInfo(i).suffix]);
        end
    end
    % RBExD
    if isfield(resultGUI,'RBExD')
        for i = 1:numBeams+1 
            result.(['RBExD', beamInfo(i).suffix ]) = numOfFractions.*resultGUI.(['RBExD', beamInfo(i).suffix]);
        end
    end
    
    end