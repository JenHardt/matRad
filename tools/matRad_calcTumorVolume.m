function output = matRad_calcTumorVolume(ct,cst,field)
% volumePTV = Volume of the PTV segmentation in [mL] or [cc]
matRad_cfg = MatRad_Config.instance();
if nargin < 3
    field = 'PTV';
end
idx = cst(contains(cst(:,2),field),[2,4]);

if ~isempty(idx)
    voxelVolume = ct.resolution.x * ct.resolution.y * ct.resolution.z;
    volumePTV = cellfun(@numel,[idx{:,2}]') * voxelVolume / 1000;

    output = [idx(:,1),num2cell(volumePTV)];

    tabHeader = {'Structure','Volume [cc]'};
    outTable = table(output(:,1),output(:,2),'VariableNames',tabHeader);
else
    matRad_cfg.dispWarning(['Field ',field,' not found in cst.']);
    output = [];
    outTable = [];
end
end

