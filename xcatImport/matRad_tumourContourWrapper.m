function [structures,ct] = matRad_tumourContourWrapper(structures,ct,importOptions)

origNumStruct = size(structures,2);
offset = origNumStruct+1;


%% GTV
radius = ct.tumourMotion.radius;
tumourMotion = false;
indGTV = matRad_tumourContour(ct,radius,tumourMotion);

structures(offset).structName = 'GTV';
structures(offset).indices = indGTV;
offset = offset+1;

%% CTV
radius = ct.tumourMotion.radius+importOptions.margins.CTV;
tumourMotion = false;
indCTV = matRad_tumourContour(ct,radius,tumourMotion);

structures(offset).structName = 'CTV';
structures(offset).indices = indCTV;
offset = offset+1;

%% ITV
radius = ct.tumourMotion.radius+importOptions.margins.CTV;
tumourMotion = true;
indITV = matRad_tumourContour(ct,radius,tumourMotion);

structures(offset).structName = 'ITV';
structures(offset).indices = indITV;
offset = offset+1;

%% PTV_CTV
radius = ct.tumourMotion.radius+importOptions.margins.CTV+importOptions.margins.PTV;
tumourMotion = false;
indPTV_CTV = matRad_tumourContour(ct,radius,tumourMotion);

structures(offset).structName = 'PTV_CTV';
structures(offset).indices = indPTV_CTV;
offset = offset+1;

%% PTV_ITV
radius = ct.tumourMotion.radius+importOptions.margins.CTV+importOptions.margins.PTV;
tumourMotion = true;
indPTV_ITV = matRad_tumourContour(ct,radius,tumourMotion);

structures(offset).structName = 'PTV_ITV';
structures(offset).indices = indPTV_ITV;
offset = offset+1;

%% delete/rename motionVec
if isfield(ct,'motionVecX_new') && ~importOptions.keepAllFrames
    ct.motionVecX   = ct.motionVecX_new;
    ct.motionVecY   = ct.motionVecY_new;
    ct.motionVecZ   = ct.motionVecZ_new;
end
ct = rmfield(ct,{'motionVecX_new' 'motionVecY_new' 'motionVecZ_new'});

newNumStruct = size(structures,2);

deleteInd = false(newNumStruct,1);
for i = 1:origNumStruct
    if ~isempty(regexpi(structures(i).structName,'tv')) || ...
       ~isempty(regexpi(structures(i).structName,'target')) || ...
       ~isempty(regexpi(structures(i).structName,'gtv')) || ...
       ~isempty(regexpi(structures(i).structName,'ctv')) || ...
       ~isempty(regexpi(structures(i).structName,'ptv')) || ...
       ~isempty(regexpi(structures(i).structName,'boost')) || ...
       ~isempty(regexpi(structures(i).structName,'tumor'))
        
        deleteInd(i) = true;
    end
end

structures(deleteInd) = [];

end

