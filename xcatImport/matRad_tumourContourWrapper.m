function [structures,ct] = matRad_tumourContourWrapper(structures,ct,importOptions)

offset = size(structures,2)+1;

% should include something to delete any old contours with CTV, etc., in
% the name
% also get rid of the _new suffixes

%% GTV
radius = ct.tumourMotion.radius;
tumourMotion = false;
indGTV = matRad_tumourContour(ct,radius,tumourMotion);

structures(offset).structName = 'GTV_new';
structures(offset).indices = indGTV;
offset = offset+1;

%% CTV
radius = ct.tumourMotion.radius+importOptions.margins.CTV;
tumourMotion = false;
indCTV = matRad_tumourContour(ct,radius,tumourMotion);

structures(offset).structName = 'CTV_new';
structures(offset).indices = indCTV;
offset = offset+1;

%% ITV 
radius = ct.tumourMotion.radius+importOptions.margins.CTV;
tumourMotion = true;
indITV = matRad_tumourContour(ct,radius,tumourMotion);

structures(offset).structName = 'ITV_new';
structures(offset).indices = indITV;
offset = offset+1;

%% PTV_CTV
radius = ct.tumourMotion.radius+importOptions.margins.CTV+importOptions.margins.PTV;
tumourMotion = false;
indPTV_CTV = matRad_tumourContour(ct,radius,tumourMotion);

structures(offset).structName = 'PTV_CTV_new';
structures(offset).indices = indPTV_CTV;
offset = offset+1;

%% PTV_ITV
radius = ct.tumourMotion.radius+importOptions.margins.CTV+importOptions.margins.PTV;
tumourMotion = true;
indPTV_ITV = matRad_tumourContour(ct,radius,tumourMotion);

structures(offset).structName = 'PTV_ITV_new';
structures(offset).indices = indPTV_ITV;
offset = offset+1;

%% delete/rename motionVec
ct.motionVecX   = ct.motionVecX_new;
ct.motionVecY   = ct.motionVecY_new;
ct.motionVecZ   = ct.motionVecZ_new;
ct = rmfield(ct,{'motionVecX_new' 'motionVecY_new' 'motionVecZ_new'});

end

