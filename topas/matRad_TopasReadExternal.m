function resultGUI = matRad_TopasReadExternal(folder)

matRad_cfg =  MatRad_Config.instance();

if contains(folder,'*')% && all(modulation ~= false)
    folder = dir(folder);
    for i = 1:length(folder)
        if any(isstrprop(folder(i).name(end-4:end),'digit'))
            folders{i} = [folder(i).folder filesep folder(i).name];
        end
    end
else
    folders{1} = folder;
end
folders = folders(~cellfun('isempty',folders));

for f = 1:length(folders)
    try
        currFolder = folders{f};
        topasConfig = MatRad_TopasConfig();
        topasCubes = matRad_readTopasData(currFolder);

        load([folder filesep 'MCparam.mat']);
        MCparam.tallies = unique(MCparam.tallies);

        % check if all fields are filled.
        if any(structfun(@isempty, topasCubes)) || any(structfun(@(x) all(x(:)==0), topasCubes))
            matRad_cfg.dispWarning('Field in topasCubes resulted in either empty or all zeros.')
        end

        loadedVars = load([currFolder filesep 'dij.mat'],'dij');
        dij = loadedVars.dij;
        loadedVars = load([currFolder filesep 'weights.mat'],'w');
        w = loadedVars.w;

        ctScen = 1;
        fnames = fieldnames(topasCubes);
        dij.MC_tallies = fnames;

        if ~any(contains(fnames,'alpha','IgnoreCase',true))
            for c = 1:numel(fnames)
                dij.(fnames{c}){ctScen,1} = sum(w(:,ctScen))*reshape(topasCubes.(fnames{c}),[],1);
            end
        else
            for d = 1:dij.numOfBeams
                doseFields = MCparam.tallies(contains(MCparam.tallies,'dose','IgnoreCase',true));
                for j = 1:numel(doseFields)
                    dij.(doseFields{j}){ctScen,1}(:,d)             = sum(w)*reshape(topasCubes.([doseFields{j} '_beam',num2str(d)]),[],1);
                    dij.([doseFields{j} '_std']){ctScen,1}(:,d)          = sum(w)*reshape(topasCubes.([doseFields{j} '_std_beam',num2str(d)]),[],1);
                end

                abFields = MCparam.tallies(contains(MCparam.tallies,{'alpha','beta','LET'},'IgnoreCase',true));
                for j = 1:numel(abFields)
                    dij.(abFields{j}){ctScen,1}(:,d)        = reshape(topasCubes.([abFields{j} '_beam',num2str(d)]),[],1);
                end

                % find first available model for RBE evaluation
                for j = 1:numel(topasConfig.scorerRBEmodelOrderForEvaluation)
                    if any(contains(fnames,topasConfig.scorerRBEmodelOrderForEvaluation{j}))
                        model = topasConfig.scorerRBEmodelOrderForEvaluation{j};
                        break
                    end
                end

                % Either dose to water or dose to medium (physical Dose) used to calculate alpha and sqrt(beta) doses
                dij.mAlphaDose{ctScen,1}(:,d)           = dij.(['alpha_' model]){ctScen,1}(:,d) .* dij.physicalDose{ctScen,1}(:,d);             
                dij.mSqrtBetaDose{ctScen,1}(:,d)        = sqrt(dij.(['beta_' model]){ctScen,1}(:,d)) .* dij.physicalDose{ctScen,1}(:,d);
%                 dij.mAlphaDose{ctScen,1}(:,d)           = dij.(['alpha_' model]){ctScen,1}(:,d) .* dij.doseToWater{ctScen,1}(:,d);             
%                 dij.mSqrtBetaDose{ctScen,1}(:,d)        = sqrt(dij.(['beta_' model]){ctScen,1}(:,d)) .* dij.physicalDose{ctScen,1}(:,d);
            end
        end

        if length(folders) > 1
            outDose    = matRad_calcCubes(ones(dij.numOfBeams,1),dij,1);
            if ~exist('resultGUI')
                for i = 1:dij.numOfBeams
                    beamInfo(i).suffix = ['_beam', num2str(i)];
                end
                beamInfo(dij.numOfBeams+1).suffix = '';
                for i = 1:length(beamInfo)
                    resultGUI.(['physicalDose', beamInfo(i).suffix]) = zeros(dij.ctGrid.dimensions);
                    resultGUI.(['RBExD', beamInfo(i).suffix]) = zeros(dij.ctGrid.dimensions);

                    resultGUI.(['alpha', beamInfo(i).suffix]) = {};
                    resultGUI.(['beta', beamInfo(i).suffix]) = {};
                    resultGUI.(['RBE', beamInfo(i).suffix]) = {};
                    resultGUI.(['effect', beamInfo(i).suffix]) = {};
                end
            end

            for i = 1:length(beamInfo)
                resultGUI.(['physicalDose', beamInfo(i).suffix]) = resultGUI.(['physicalDose', beamInfo(i).suffix]) + outDose.(['physicalDose', beamInfo(i).suffix])/length(folders);
                if isfield(outDose,'alpha')
                    resultGUI.(['RBExD', beamInfo(i).suffix]) = resultGUI.(['RBExD', beamInfo(i).suffix]) + outDose.(['RBExD', beamInfo(i).suffix])/length(folders);
                    resultGUI.(['alpha', beamInfo(i).suffix]){f} = outDose.(['alpha', beamInfo(i).suffix]);
                    resultGUI.(['beta', beamInfo(i).suffix]){f} = outDose.(['beta', beamInfo(i).suffix]);
                    resultGUI.(['RBE', beamInfo(i).suffix]){f} = outDose.(['RBE', beamInfo(i).suffix]);
                    resultGUI.(['effect', beamInfo(i).suffix]){f} = outDose.(['effect', beamInfo(i).suffix]);
                end
                resultGUI.samples = f;
            end
        else
            resultGUI    = matRad_calcCubesMC(ones(dij.numOfBeams,1),dij,1);
        end
    catch ME
        subFolder = strsplit(currFolder,'\');
        subFolder = subFolder{end};
        matRad_cfg.dispError(['Error in folder "',subFolder,'": ',strjoin(strsplit(ME.message,'\'),'/')]);
    end
end

end