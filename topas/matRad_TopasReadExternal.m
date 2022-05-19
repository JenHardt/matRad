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

        d = dir([currFolder filesep 'score_matRad_plan_field1_run1_*']);
        [~,~,outputType] = fileparts(d(1).name);
        files = dir([currFolder filesep 'score_matRad_plan_field1_run1_*' outputType]);
        tallies = cellfun(@(x) erase(x,{'score_matRad_plan_field1_run1_', outputType}) ,{files(:).name} ,'UniformOutput' ,false);

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

        for d = 1:dij.numOfBeams
                doseFields = tallies(contains(tallies,'dose','IgnoreCase',true));
                for j = 1:numel(doseFields)
                    dij.(doseFields{j}){ctScen,1}(:,d)             = sum(w)*reshape(topasCubes.([doseFields{j} '_beam',num2str(d)]),[],1);
                    dij.([doseFields{j} '_std']){ctScen,1}(:,d)          = sum(w)*reshape(topasCubes.([doseFields{j} '_std_beam',num2str(d)]),[],1);
                    if any(contains(fieldnames(topasCubes),'batchStd'))
                        dij.([doseFields{j} '_batchStd']){ctScen,1}(:,d)          = sum(w)*reshape(topasCubes.([doseFields{j} '_batchStd_beam',num2str(d)]),[],1);
                    end
                end
        end
        
        if any(contains(fnames,'alpha','IgnoreCase',true))
            for d = 1:dij.numOfBeams
                abFields = tallies(contains(tallies,{'alpha','beta','LET'},'IgnoreCase',true));
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
                dij.(['mAlphaDose_' model]){ctScen,1}(:,d)           = dij.(['alpha_' model]){ctScen,1}(:,d) .* dij.physicalDose{ctScen,1}(:,d);
                dij.(['mSqrtBetaDose_' model]){ctScen,1}(:,d)        = sqrt(dij.(['beta_' model]){ctScen,1}(:,d)) .* dij.physicalDose{ctScen,1}(:,d);
                %                 dij.mAlphaDose{ctScen,1}(:,d)           = dij.(['alpha_' model]){ctScen,1}(:,d) .* dij.doseToWater{ctScen,1}(:,d);
                %                 dij.mSqrtBetaDose{ctScen,1}(:,d)        = sqrt(dij.(['beta_' model]){ctScen,1}(:,d)) .* dij.physicalDose{ctScen,1}(:,d);
            end
        end

        if length(folders) > 1
            outDose    = matRad_calcCubesMC(ones(dij.numOfBeams,1),dij,1);
            %             if any(contains(fieldnames(outDose),'alpha'))
            %                 names = fieldnames(outDose);
            %                 alphaFields = names((contains(names,{'alpha','beta','effect','RBE'}) + ~contains(names,{'RBExD'})) > 1);
            %             end
            if ~exist('resultGUI')
                for i = 1:dij.numOfBeams
                    beamInfo(i).suffix = ['_beam', num2str(i)];
                end
                beamInfo(dij.numOfBeams+1).suffix = '';
                for i = 1:length(beamInfo)
                    resultGUI.(['physicalDose', beamInfo(i).suffix]) = zeros(dij.ctGrid.dimensions);
                    if any(contains(fieldnames(outDose),'alpha'))
                        resultGUI.(['RBExD_' model beamInfo(i).suffix]) = zeros(dij.ctGrid.dimensions);
                    end
                    %
                    %                     for k = 1:length(alphaFields)
                    %                         resultGUI.(alphaFields{k}) = cell(1,length(folders));
                    %                     end
                    resultGUI = orderfields(resultGUI);
                end
            end

            for i = 1:length(beamInfo)
                resultGUI.(['physicalDose', beamInfo(i).suffix]) = resultGUI.(['physicalDose', beamInfo(i).suffix]) + outDose.(['physicalDose', beamInfo(i).suffix])/length(folders);
            end
            if any(contains(fieldnames(outDose),'alpha'))
                for i = 1:length(beamInfo)
                    resultGUI.(['RBExD_' model beamInfo(i).suffix]) = resultGUI.(['RBExD_' model beamInfo(i).suffix]) + outDose.(['RBExD_' model beamInfo(i).suffix])/length(folders);
                end
            end
            resultGUI.samples = f;
        else
            resultGUI    = matRad_calcCubesMC(ones(dij.numOfBeams,1),dij,1);
        end
    catch ME
        subFolder = strsplit(currFolder,'\');
        subFolder = subFolder{end};
        matRad_cfg.dispError(['Error in line ' num2str(ME.stack(end).line) '\nError in folder "',subFolder,'": ',strjoin(strsplit(ME.message,'\'),'/')]);
    end
end

end