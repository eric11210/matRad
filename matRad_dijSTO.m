function [dij_STO,trajectory] = matRad_dijSTO(dij,pln,stf)

%% select most probable trajectory

% prepare model
motionModel = matRad_prepModelForOpt(pln.propOpt.prop4D);

% load machine properties
fileName = pln.propOpt.VMAToptions.machineConstraintFile;
try
    load(fileName,'machine');
catch
    error(['Could not find the following machine file: ' fileName ]);
end

% determine transition times
tTrans = zeros(pln.propStf.numOfBeams,1);

for i = 1:pln.propStf.numOfBeams
    
    tTrans(i) = stf(1).propVMAT.doseAngleBordersDiff./machine.constraints.gantryRotationSpeed(2);
end

% determine index map, set maxProb to true
indexMap    = (1:numel(motionModel.initProb))';
maxProb     = true;

% sample the most probable trajectory
[lSimulated,tSimulated] = matRad_runMarkovChain_Q(motionModel,tTrans,indexMap,maxProb);

% convert subphase to phase
pSimulated = motionModel.indices.subPhase2PosPhase(lSimulated);

% store trajectory information
trajectory.tSimulated = tSimulated;
trajectory.pSimulated = pSimulated;

%% construct effective dij

dij_fieldNames = fieldnames(dij);

for i = 1:numel(dij_fieldNames)
    
    if ~strcmp(dij_fieldNames{i},'physicalDose')
        
        dij_STO.(dij_fieldNames{i}) = dij.(dij_fieldNames{i});
    end
end

dij_STO.physicalDose{1} = spalloc(dij.numOfVoxels,dij.totalNumOfBixels,1);
dij.numOfScenarios      = 1;

offset = 0;
for i = 1:pln.propStf.numOfBeams
    
    dij_STO.physicalDose{1}(:,offset+(1:stf(i).totalNumOfBixels)) = dij.physicalDose{pSimulated(i)}(:,offset+(1:stf(i).totalNumOfBixels));
    
    offset = offset+stf(i).totalNumOfBixels;
end

end