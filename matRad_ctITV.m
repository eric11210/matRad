function ct_ITV = matRad_ctITV(ct,pln)

%% calculate probability for each phase

% prepare model
motionModel = matRad_prepModelForOpt(pln.propOpt.prop4D);

% determine probabilities
probPhase = accumarray(motionModel.indices.subPhase2PosPhase,motionModel.Pi_deltaTSample);

%% construct ITV ct

ct_fieldNames = fieldnames(ct);

for i = 1:numel(ct_fieldNames)
    
    if ~(strcmp(ct_fieldNames{i},'cube') || strcmp(ct_fieldNames{i},'cubeHU'))
        
        ct_ITV.(ct_fieldNames{i}) = ct.(ct_fieldNames{i});
    end
end

% prep ct_ITV cubes
ct_ITV.numOfCtScen  = 1;
ct_ITV.cube{1}      = zeros(ct_ITV.cubeDim);
ct_ITV.cubeHU{1}    = zeros(ct_ITV.cubeDim);

% take mean of ct cubes
for phase = 1:ct.numOfCtScen
    
    ct_ITV.cube{1}      = ct_ITV.cube{1} + probPhase(phase).*ct.cube{phase};
    ct_ITV.cubeHU{1}    = ct_ITV.cubeHU{1} + probPhase(phase).*ct.cubeHU{phase};
end

end